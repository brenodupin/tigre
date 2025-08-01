# -*- coding: utf-8 -*-
"""Module to solve overlaps in GFF3 files."""

import concurrent.futures as cf
from pathlib import Path
from typing import Callable

import pandas as pd

from . import gff3_utils, log_setup

species_url = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="


def multiple_regions_solver(
    log: log_setup.TempLogger,
    df: pd.DataFrame,
    an: str,
) -> pd.DataFrame:
    log.warning("Multiple regions found:")
    for row in df[df["type"] == "region"].itertuples():
        log.warning(f" s: {row.start} e: {row.end} | {row.attributes}")

    region_mask = df["type"] == "region"
    if region_mask.sum() > 1:
        df = df[~region_mask | (df.index == 0)].reset_index(drop=True)
    log.warning(
        f" Chosen s: {df.at[0, 'start']} e: {df.at[0, 'end']} | {df.at[0, 'attributes']}"
    )
    return df


def bigger_than_region_solver(
    log: log_setup.TempLogger,
    df: pd.DataFrame,
    an: str,
) -> pd.DataFrame:
    region_end = df.at[0, "end"]
    log.warning(f"Genes with 'end' bigger than region: {region_end}")
    for row in df[df["end"] > region_end].itertuples():
        log.warning(f" s: {row.start} e: {row.end} | t: {row.type} | {row.attributes}")

    new_rows = []
    for index, row in df[df["end"] > region_end].iterrows():
        row_type = row.type
        log.debug(f"Processing type {row_type}: {row.attributes}")

        if row.end // region_end > 1:
            log.error(f" Gene is bigger than region! g:{row.end} r:{region_end}")

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

        log.debug(f" s: {new_row.start} | e: {new_row.end} | {new_row.attributes}")

    return pd.concat([df, pd.DataFrame(new_rows)])


def overlaps_chooser(
    log: log_setup.TempLogger,
    overlaps_in: list[pd.Series | pd.DataFrame],
    an: str,
) -> str:
    overlaps = pd.DataFrame(overlaps_in).sort_values("start", ignore_index=True)

    left = overlaps[overlaps.start == overlaps.start.min()]
    if left.shape[0] > 1:
        log.trace(
            f" More than one element with the lowest start: {overlaps.start.min()}"
        )
        left = left[left.end == left.end.max()]
        if left.shape[0] > 1 and left.end.max() == overlaps.end.max():
            log.trace(
                f" Genes with same start and end, choosing the first one: {left.iat[0, 8]}"
            )
            return f"{name_replace(left.iat[0, 8], is_left=True)}{name_replace(left.iat[0, 8])}"

    right = overlaps[overlaps.end == overlaps.end.max()]
    if right.shape[0] > 1:
        log.trace(f" More than one element with the highest end: {overlaps.end.max()}")
        right = right[right.start == right.start.min()]
    return (
        f"{name_replace(left.iat[0, 8], is_left=True)}{name_replace(right.iat[0, 8])}"
    )


def name_replace(
    string: str,
    is_left: bool = False,
) -> str:
    if is_left:
        return string.replace("name=", "name_left=").replace("source=", "source_left=")
    return string.replace("name=", "name_right=").replace("source=", "source_right=")


def create_header(
    an: str,
    attributes: str,
    end: int,
) -> str:
    header = f"##gff-version 3\n#!gff-spec-version 1.26\n#!processor [PLACE_H0LDER_NAME]\n##sequence-region {an} 1 {end}\n"
    taxon_match = gff3_utils._RE_region_taxon.search(attributes)
    header += f"##species {species_url + taxon_match.group(1) if taxon_match else 'taxon_id_not_found'}\n"
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
                f"overlap found, region: {row.start} - {end_region} | length: {end_region - row.start + 1}"
            )
            for r in overlaps:
                log.trace(f" {r.start} | {r.end} | {r.type} | {r.attributes}")

            # modify the current row with overlap_region info
            row.end = end_region
            row.type = "overlapping_feature_set"
            # row.strand = '+'
            row.attributes = overlaps_chooser(log, overlaps, an)
            log.debug(f"Result: {row.attributes}")

        # append the row to the df_clean, be it modified (and therefor overlap_region) or not
        df_clean = pd.concat([df_clean, row.to_frame().T], ignore_index=True)

    return df_clean


def clean_attr(
    row: pd.Series,
    log: log_setup.TempLogger,
) -> str:
    return f"name={row.gene_id};source={row.seqid}|{row.type}|{row.start}|{row.end}|{row.strand}|{row.gene_id};"


def clean_an(
    log: log_setup.TempLogger,
    gff_in: Path,
    gff_out: Path,
    clean_func: Callable[[pd.Series, log_setup.TempLogger], str],
    query_string: str,
    keep_orfs: bool,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    try:
        log.info(f"[{gff_in.name}] -- Start --")

        log.trace(f" Loading GFF3 file: {gff_in}")
        df = gff3_utils.load_gff3(
            gff_in,
            usecols=gff3_utils.GFF3_COLUMNS,
            query_string=query_string,
        )
        if not keep_orfs:
            df = gff3_utils.filter_orfs(df)
        an = df.at[0, "seqid"]
        header = create_header(an, df.at[0, "attributes"], df.at[0, "end"])
        df["gene_id"] = df["attributes"].str.extract(gff3_utils._RE_ID, expand=False)  # type: ignore[call-overload]

        if (df["type"] == "region").sum() > 1:
            df = multiple_regions_solver(log, df, an)

        log.trace(f" Cleaning function: {clean_func}.")
        df.loc[1:, "attributes"] = df.loc[1:].apply(clean_func, axis=1, log=log)

        if df["end"].idxmax() != 0:
            df = bigger_than_region_solver(log, df, an)

        log.trace(" Sorting DataFrame by start and type...")
        df.at[df[df["type"] == "region"].index[0], "type"] = "_"
        df = df.sort_values(by=["start", "type"], ascending=True, ignore_index=True)
        df.at[0, "type"] = "region"

        log.trace(" Overlaps solver...")
        df = overlap_solver(log, df, an)
        df = df.drop(columns=["gene_id"], errors="ignore")

        # add header to the file
        with open(gff_out, "w", encoding="utf-8") as file_handler:
            file_handler.write(header)
            df.to_csv(file_handler, sep="\t", header=False, index=False)

        log.info(f"[{an}] -- End --")
        return True, an, log.get_records()

    except Exception as e:
        log.error(f"[{an}] Error: {e}")
        log.info(f"[{an}] -- End --")
        return False, an, log.get_records()


def clean_multiple(
    log: log_setup.GDTLogger,
    tsv_path: Path,
    workers: int,
    gff_in_ext: str = ".gff3",
    gff_in_suffix: str = "",
    gff_out_ext: str = ".gff3",
    gff_out_suffix: str = "_clean",
    an_column: str = "AN",
    clean_func: Callable[[pd.Series, log_setup.TempLogger], str] = clean_attr,
    query_string: str = gff3_utils.QS_GENE_TRNA_RRNA_REGION,
    keep_orfs: bool = False,
    overwrite: bool = False,
) -> None:
    """Orchestrates the execution of the ans."""
    tsv = pd.read_csv(tsv_path, sep="\t")

    gff_in_builder = gff3_utils.PathBuilder(gff_in_ext).use_folder_builder(
        tsv_path.parent,
        gff_in_suffix,
    )
    gff_out_builder = gff3_utils.PathBuilder(gff_out_ext).use_folder_builder(
        tsv_path.parent,
        gff_out_suffix,
    )
    gff3_utils.check_tsv(
        log, tsv, gff_in_builder, gff_out_builder, overwrite, an_column
    )

    log.info(f"Starting processing {tsv.shape[0]} ANs with {workers} workers...")
    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        tasks = [
            executor.submit(
                clean_an,
                log.spawn_buffer(),
                gff_in_builder.build(an),
                gff_out_builder.build(an),
                clean_func,
                query_string,
                keep_orfs,
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
