# -*- coding: utf-8 -*-

import concurrent.futures as cf
from pathlib import Path
from typing import Any

import pandas as pd

from . import gff3_utils, log_setup


def multiple_regions_solver(
    df: pd.DataFrame, an: str, log: log_setup.GDTLogger
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
    df: pd.DataFrame, an: str, log: log_setup.GDTLogger
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
    overlaps_in: list[pd.Series | pd.DataFrame], an: str, log: log_setup.GDTLogger
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


def name_replace_bool(s: str, is_left: bool) -> str:
    if is_left:
        return s.replace("name=", "name_left=").replace("source=", "source_left=")
    return s.replace("name=", "name_right=").replace("source=", "source_right=")


def create_header(region_att: str, an: str, region_end: int) -> str:
    header = f"##gff-version 3\n#!gff-spec-version 1.26\n#!processor [PLACE_H0LDER_NAME]\n##sequence-region {an} 1 {region_end}\n"
    taxon_start = region_att.find("taxon:") + 6
    header += f"##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={region_att[taxon_start : region_att.find(';', taxon_start)]}\n"
    return header


def overlap_solver(df: pd.DataFrame, an: str, log: log_setup.GDTLogger) -> pd.DataFrame:
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
            row.attributes = overlaps_chooser(overlaps, an, log)
            log.debug(f"[{an}] Att: {row.attributes}")

        # append the row to the df_clean, be it modified (and therefor overlap_region) or not
        df_clean = pd.concat([df_clean, row.to_frame().T], ignore_index=True)

    return df_clean


def clean_attr(
    row: pd.Series, an: str, gene_dict: Any, log: log_setup.GDTLogger
) -> str:
    gene_id = row.attributes.split(";", 1)[0].split("ID=", 1)[1]
    try:
        gene_label = gene_dict[gene_id].label
    except KeyError:
        #log.error(f"[{an}] No gene_dict for {gene_id}")
        return f"name=NO_GENE_DICT_{gene_id};source={row.seqid}|{row.type}|{row.start}|{row.end}|{row.strand}|{gene_id};"

    return f"name={gene_label};source={row.seqid}|{row.type}|{row.start}|{row.end}|{row.strand}|{gene_id};"


def exec_single_an(
    an: str, an_path: Path, clean_path: Path, log: log_setup.GDTLogger
) -> tuple[bool, str]:
    try:
        log.info(f"[{an}] -- Start --")

        log.trace(f"[{an}] Loading GFF3 file: {an_path}")
        df = gff3_utils.load_gff3(
            an_path,
            usecols=gff3_utils.GFF3_COLUMNS,
            query_string='type == ["gene", "region"]',
        )
        df = gff3_utils.filter_orfs(df)
        header = create_header(df.at[0, "attributes"], an, df.at[0, "end"])

        if (df["type"] == "region").sum() > 1:
            log.trace(f"[{an}] More than one region found, solving...")
            df = multiple_regions_solver(df, an, log)

        log.trace(f"[{an}] Cleaning attributes...")
        df.loc[1:, "attributes"] = df.loc[1:].apply(clean_attr, axis=1, an=an, gene_dict={}, log=log)

        if df["end"].idxmax() != 0:
            log.trace(f"[{an}] Region end is not the maximum, fixing...")
            df = bigger_than_region_solver(df, an, log)

        log.trace(f"[{an}] Sorting DataFrame by start and type...")
        df.at[df[df["type"] == "region"].index[0], "type"] = "_"
        df = df.sort_values(by=["start", "type"], ascending=True, ignore_index=True)
        df.at[0, "type"] = "region"

        log.trace(f"[{an}] Overlaps solver...")
        df = overlap_solver(df, an, log)

        # add header to the file
        with open(clean_path, "w") as f:
            f.write(header)
            df.to_csv(f, sep="\t", header=False, index=False)

        log.info(f"[{an}] -- End --")
        return True, an
    except Exception as e:
        log.error(f"[{an}] Error: {str(e)}")
        log.info(f"[{an}] -- End --")
        return False, an


def orchestration(tsv_path: Path, log: log_setup.GDTLogger) -> None:
    """Orchestrates the execution of the ans."""
    num_workers = 38

    tsv = pd.read_csv(tsv_path, sep="\t")
    folder = tsv_path.parent

    an_builder = gff3_utils.GFFPathBuilder().use_folder_builder(folder)
    clean_builder = gff3_utils.GFFPathBuilder().use_folder_builder(
        folder, suffix="_clean"
    )

    print(f"Processing {len(tsv)} ans with {num_workers} workers")
    with cf.ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(
                exec_single_an,
                an,
                an_builder.build(an),
                clean_builder.build(an),
                log,
            )
            for an in tsv["AN"]
        ]
    cf.wait(futures)

    for future in futures:
        r = future.result()
        if not r[0]:
            print(f"Error on {r[1]}")


def clean_log(folder: str) -> None:
    time_l, msg_lvl_l, msg_l = [], [], []
    with open(f"{folder}_overlaps_raw_final.log", "r") as f:
        for line in f.readlines():
            time, msg_lvl, msg = line.split(" - ", 2)
            time_l.append(time)
            msg_lvl_l.append(msg_lvl)
            msg_l.append(msg.strip())

    df_log = pd.DataFrame({"time": time_l, "msg_lvl": msg_lvl_l, "msg": msg_l})
    df_log[["an", "msg"]] = df_log["msg"].str.split("] ", n=1, expand=True)
    df_log["an"] = df_log["an"].str.replace("[", "")
    df_log["time"] = pd.to_datetime(df_log["time"], format="%Y-%m-%d %H:%M:%S,%f")
    df_log = df_log[df_log["msg"] != "-- End --"]
    df_log = df_log[["an", "time", "msg_lvl", "msg"]]

    with open(f"{folder}_overlaps_clean_final.log", "w+") as f:

        def format_row(row: pd.Series) -> str:
            if row.msg == "-- Start --":
                return f"{row.time} - {row.msg_lvl} - [{row.an}] {row.msg}"
            return f"{row.time} - {row.msg_lvl} - {row.msg}"

        df_log["formatted"] = df_log.apply(format_row, axis=1)
        df_log["extra_newline"] = (
            df_log["msg"].str.startswith("Att: ").map({True: "\n", False: ""})
        )

        for an, group in df_log.groupby("an", sort=False):
            group_text = "\n".join(group["formatted"] + group["extra_newline"]) + "\n\n"
            f.write(group_text)
