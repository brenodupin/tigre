# -*- coding: utf-8 -*-
"""Module to solve overlaps in GFF3 files."""

import concurrent.futures as cf
import queue
from pathlib import Path
from typing import Callable

import pandas as pd

from . import gff3_utils, log_setup


def multiple_regions_solver(
    log: log_setup.TempLogger,
    df: pd.DataFrame,
    an: str,
) -> pd.DataFrame:
    log.warning(f"[{an}] Multiple regions found:")
    for row in df[df["type"] == "region"].itertuples():
        log.warning(f"[{an}] s: {row.start} e: {row.end} | {row.attributes}")

    region_mask = df["type"] == "region"
    if region_mask.sum() > 1:
        df = df[~region_mask | (df.index == 0)].reset_index(drop=True)
    log.warning(
        f"[{an}] Chosen s: {df.at[0, 'start']} e: {df.at[0, 'end']} | {df.at[0, 'attributes']}"
    )
    return df


def bigger_than_region_solver(
    log: log_setup.TempLogger,
    df: pd.DataFrame,
    an: str,
) -> pd.DataFrame:
    region_end = df.at[0, "end"]
    log.warning(f"[{an}] Genes with 'end' bigger than region: {region_end}")
    for row in df[df["end"] > region_end].itertuples():
        log.warning(
            f"[{an}] s: {row.start} e: {row.end} | t: {row.type} | {row.attributes}"
        )

    new_rows = []
    for index, row in df[df["end"] > region_end].iterrows():
        row_type = row.type
        log.debug(f"[{an}] Processing type {row_type}: {row.attributes}")

        if row.end // region_end > 1:
            log.error(f"[{an}] Gene is bigger than region! g:{row.end} r:{region_end}")

        # Handle overlap
        overlap = row.end - region_end
        df.at[index, "end"] = region_end
        df.at[index, "type"] = row_type + "_region"

        # Create new row
        new_row = row.copy()
        new_row.start = 1
        new_row.end = overlap
        new_row.type = row_type + "_region"
        new_rows.append(new_row)

        log.debug(
            f"[{an}] s: {new_row.start} | e: {new_row.end} | {new_row.attributes}"
        )

    return pd.concat([df, pd.DataFrame(new_rows)])


def overlaps_chooser(
    log: log_setup.TempLogger,
    overlaps_in: list[pd.Series | pd.DataFrame],
    an: str,
) -> str:
    overlaps = pd.DataFrame(overlaps_in).sort_values("start", ignore_index=True)

    log.debug(f"[{an}] Overlaps found: {overlaps.shape[0]}")
    left = overlaps[overlaps.start == overlaps.start.min()]
    if left.shape[0] > 1:
        log.warning(
            f"[{an}] More than one element with the lowest start: {overlaps.start.min()}"
        )
        left = left[left.end == left.end.max()]
        if left.shape[0] > 1 and left.end.max() == overlaps.end.max():
            log.warning(
                f"[{an}] Genes with same start and end, choosing the first one: {left.iat[0, 8]}"
            )
            return f"{name_replace_bool(left.iat[0, 8], is_left=True)}{name_replace_bool(left.iat[0, 8], is_left=False)}"

    right = overlaps[overlaps.end == overlaps.end.max()]
    if right.shape[0] > 1:
        log.warning(
            f"[{an}] More than one element with the highest end: {overlaps.end.max()}"
        )
        right = right[right.start == right.start.min()]
    return f"{name_replace_bool(left.iat[0, 8], is_left=True)}{name_replace_bool(right.iat[0, 8], is_left=False)}"


def name_replace_bool(
    s: str,
    is_left: bool,
) -> str:
    if is_left:
        return s.replace("name=", "name_left=").replace("source=", "source_left=")
    return s.replace("name=", "name_right=").replace("source=", "source_right=")


def create_header(
    an: str,
    region_att: str,
    region_end: int,
) -> str:
    header = f"##gff-version 3\n#!gff-spec-version 1.26\n#!processor [PLACE_H0LDER_NAME]\n##sequence-region {an} 1 {region_end}\n"
    taxon_start = region_att.find("taxon:") + 6
    header += f"##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={region_att[taxon_start : region_att.find(';', taxon_start)]}\n"
    return header


def overlap_solver(
    log: log_setup.TempLogger,
    df: pd.DataFrame,
    an: str,
) -> pd.DataFrame:
    df_clean = df.iloc[0:1].copy()  # copia region, index 0

    next_index = 1  # itera sobre as linhas do dataframe, começando da segunda
    for index, row in df.loc[1:].iterrows():
        if index is not next_index:
            continue

        end_region = row.end
        next_index += 1
        overlaps = []

        # peek next rows to check for overlaps
        while (next_index < len(df)) and (end_region >= df.at[next_index, "start"]):
            # overlap found, append to list and update next_index
            # so the while will check the next row
            overlaps.append(df.loc[next_index])

            if end_region < df.at[next_index, "end"]:
                end_region = df.at[next_index, "end"]

            next_index += 1

        if overlaps:
            overlaps.insert(0, row.copy())  # insert current row at index 0
            log.debug(
                f"[{an}] overlap found, region: {row.start} - {end_region} | length: {end_region - row.start + 1}"
            )
            for r in overlaps:
                log.trace(f"[{an}] s: {r.start} | e: {r.end} | {r.attributes}")

            # modify the current row with overlap_region info
            row.end = end_region
            row.type = "overlapping_feature_set"
            # row.strand = '+'
            row.attributes = overlaps_chooser(log, overlaps, an)
            log.debug(f"[{an}] Att: {row.attributes}")

        # append the row to the df_clean, be it modified (and therefor overlap_region) or not
        df_clean = pd.concat([df_clean, row.to_frame().T], ignore_index=True)

    return df_clean


def clean_attr(
    log: log_setup.TempLogger,
    row: pd.Series,
) -> str:
    return f"name={row.gene_id};source={row.seqid}|{row.type}|{row.start}|{row.end}|{row.strand}|{row.gene_id};"


def exec_single_an(
    log: log_setup.TempLogger,
    an: str,
    gff_in: Path,
    gff_out: Path,
    clean_func: Callable[[log_setup.TempLogger, pd.Series], str],
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    try:
        log.info(f"[{an}] -- Start --")

        log.trace(f"[{an}] Loading GFF3 file: {gff_in}")
        df = gff3_utils.load_gff3(
            gff_in,
            usecols=gff3_utils.GFF3_COLUMNS,
            query_string='type == ["gene", "region"]',
        )
        df = gff3_utils.filter_orfs(df)
        header = create_header(df.at[0, "attributes"], an, df.at[0, "end"])
        df["gene_id"] = df["attributes"].str.extract(gff3_utils._RE_ID, expand=False)  # type: ignore[call-overload]

        if (df["type"] == "region").sum() > 1:
            log.trace(f"[{an}] More than one region found, solving...")
            df = multiple_regions_solver(log, df, an)

        log.trace(f"[{an}] Cleaning attributes...")
        df.loc[1:, "attributes"] = df.loc[1:].apply(clean_func, axis=1, log=log)

        if df["end"].idxmax() != 0:
            log.trace(f"[{an}] Region end is not the maximum, fixing...")
            df = bigger_than_region_solver(log, df, an)

        log.trace(f"[{an}] Sorting DataFrame by start and type...")
        df.at[df[df["type"] == "region"].index[0], "type"] = "_"
        df = df.sort_values(by=["start", "type"], ascending=True, ignore_index=True)
        df.at[0, "type"] = "region"

        log.trace(f"[{an}] Overlaps solver...")
        df = overlap_solver(log, df, an)
        df = df.drop(columns=["gene_id"], errors="ignore")

        # add header to the file
        with open(gff_out, "w") as f:
            f.write(header)
            df.to_csv(f, sep="\t", header=False, index=False)

        log.info(f"[{an}] -- End --")
        return True, an, log.get_records()
    except Exception as e:
        log.error(f"[{an}] Error: {str(e)}")
        log.info(f"[{an}] -- End --")
        return False, an, log.get_records()


def exec_single_an_with_pool(
    pool: queue.Queue[log_setup.TempLogger],
    an: str,
    gff_in: Path,
    gff_out: Path,
    clean_func: Callable[[log_setup.TempLogger, pd.Series], str],
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Execute a single AN with a thread pool."""
    log = pool.get()
    try:
        return exec_single_an(log, an, gff_in, gff_out, clean_func)
    finally:
        log.clear()
        pool.put(log)


def orchestration(
    log: log_setup.GDTLogger,
    tsv_path: Path,
    gff_ext: str = ".gff3",
    gff_suffix: str = "",
    out_suffix: str = "_clean",
    an_column: str = "AN",
    workers: int = 0,
    clean_func: Callable[[log_setup.TempLogger, pd.Series], str] = clean_attr,
) -> None:
    """Orchestrates the execution of the ans."""
    tsv = pd.read_csv(tsv_path, sep="\t")

    gff_in_builder = gff3_utils.PathBuilder(gff_ext).use_folder_builder(
        tsv_path.parent,
        gff_suffix,
    )
    gff_out_builder = gff3_utils.PathBuilder(gff_ext).use_folder_builder(
        tsv_path.parent,
        out_suffix,
    )
    gff3_utils.check_file_in_tsv(
        log,
        tsv,
        gff_in_builder,
        an_column,
    )

    logger_pool: queue.Queue[log_setup.TempLogger] = queue.Queue()
    for _ in range(workers + 1):
        logger_pool.put(log.spawn_buffer())

    log.info(f"Starting processing {tsv.shape[0]} ANs with {workers} workers...")
    with cf.ThreadPoolExecutor(max_workers=workers) as executor:
        tasks = [
            executor.submit(
                exec_single_an_with_pool,
                logger_pool,
                an,
                gff_in_builder.build(an),
                gff_out_builder.build(an),
                clean_func,
            )
            for an in tsv[an_column]
        ]

        for future in cf.as_completed(tasks):
            success, an, records = future.result()
            if not success:
                log.error(f"[{an}] Failed to process.")

            for record in records:
                log.log(record[0], record[1])

        log.info("All processed.")
