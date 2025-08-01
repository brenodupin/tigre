# -*- coding: utf-8 -*-
"""Module to extract intergenic regions from GFF3 files."""

import concurrent.futures as cf
import re
from pathlib import Path
from typing import Union

import pandas as pd

from . import gff3_utils, log_setup

START_IDX = 3
END_IDX = 4

NAME_LEFT_IDX = 9  # df['name_left']
SOURCE_LEFT_IDX = 10  # df['source_left']

NAME_RIGHT_IDX = 11  # df['name_right']
SOURCE_RIGHT_IDX = 12  # df['source_right']

_RE_name = re.compile(r"name=([^;]+)")
_RE_source = re.compile(r"source=([^;]+)")

_RE_name_left = re.compile(r"name_left=([^;]+)")
_RE_source_left = re.compile(r"source_left=([^;]+)")

_RE_name_right = re.compile(r"name_right=([^;]+)")
_RE_source_right = re.compile(r"source_right=([^;]+)")


def extract_intergenic_regions(
    log: log_setup.TempLogger,
    gff_in: Path,
    gff_out: Path,
    add_region: bool = False,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Extract intergenic regions from a GFF3 file and save to a new file.

    Args:
        log (log_setup.GDTLogger): Logger instance for logging.
        gff_in (Path): Input GFF3 file path.
        gff_out (Path): Output GFF3 file path.

    """
    try:
        log.debug(f"Extracting intergenic regions from {gff_in} to {gff_out}...")

        header = []
        with open(gff_in, "r") as gff_file:
            for line in gff_file:
                if not line.startswith("#"):
                    break

                header.append(line.strip())

            df = gff3_utils.load_gff3(
                gff_in,
                usecols=gff3_utils.GFF3_COLUMNS,
            )

        # pop region line out of df
        region = df.iloc[0]
        df = df.iloc[1:].reset_index(drop=True)

        region_size = region.iat[END_IDX]
        seqid = region.at["seqid"]
        source = region.at["source"]
        circular = "is_circular=true" in region.at["attributes"].lower()

        log.debug(
            f"AN: {seqid} | Source: {source} | Size: {region_size} | Circular: {circular}"
        )

        df = _split_attributes(log, df)

        boundary = _solve_boundaries(
            log,
            df,
            region_size,
            seqid,
            source,
            circular,
        )

        # this will alter df! by adding +1 to end and -1 to start
        df_ig = _create_intergenic(
            log,
            df,
            seqid,
            source,
        )

        if df_ig.empty and not boundary:
            log.error(
                f"No intergenic regions found. Did {gff_in} contain any features?"
            )
            return False, seqid, log.get_records()

        if boundary:
            df_ig = (
                pd.concat([df_ig, pd.DataFrame(boundary)])
                .sort_values(by=["start", "end"])
                .reset_index(drop=True)
            )

        start = 2 if "intergenic_region_merged" in df_ig["type"].values else 1
        df_ig["attributes"] = [
            f"ID={i};{attr}" for i, attr in enumerate(df_ig["attributes"], start)
        ]

        log.trace(f"Writing intergenic regions to {gff_out}.")
        with open(gff_out, "w") as f:
            f.write("\n".join(header) + "\n")
            if add_region:
                log.trace("Adding region line to output.")
                f.write(
                    f"{seqid}\t{source}\tregion\t1\t{region_size}\t.\t+\t.\t"
                    f"{region.attributes}\n"
                )
            df_ig.to_csv(f, sep="\t", index=False, header=False)

        return True, seqid, log.get_records()
    except Exception as e:
        log.error(f"Error extracting intergenic regions: {e}")
        return False, "", log.get_records()


def _solve_boundaries(
    log: log_setup.TempLogger,
    df: pd.DataFrame,
    region_size: int,
    seqid: str,
    source: str,
    circular: bool,
) -> list[dict[str, Union[str, int]]]:
    # Legend:
    # [---] feature (gene, trna, rrna)
    # [== igr ==] intergenic region
    # || genome boundary, i.e. position where region_size rolls over to 1

    first_start = df.iat[0, START_IDX]  # of a feature
    last_end = df.iat[-1, END_IDX]  # of a feature
    log.debug(f"First start: {first_start}, Last end: {last_end}")

    # ...---][== igr ==]||[== igr ==][---...
    # if the first start is greater than 1 and the last end is less than
    # the region size, we have an intergenic region that spans the whole
    # region, so we create an intergenic region from first_start to
    # last_end. This only happens if the region is circular, otherwise
    # we skip the merge and let the other checks handle the intergenic regions.
    if circular and first_start > 1 and last_end < region_size:
        return [
            {
                "seqid": seqid,
                "source": source,
                "type": "intergenic_region_merged",
                "start": last_end + 1,
                "end": region_size + first_start - 1,
                "score": ".",
                "strand": "+",
                "phase": ".",
                "attributes": f"name_up={df.iat[-1, NAME_RIGHT_IDX]};"
                f"source_up={df.iat[-1, SOURCE_RIGHT_IDX]};"
                f"name_dw={df.iat[0, NAME_LEFT_IDX]};"
                f"source_dw={df.iat[0, SOURCE_LEFT_IDX]};",
            }
        ]

    regions: list[dict[str, Union[str, int]]] = []
    # ...||[== igr ==][---...
    # if the first start is greater than 1, we have an intergenic region
    # before the first feature, so we create an intergenic region from 1 to
    # first_start - 1.
    if first_start > 1:
        regions.append(
            {
                "seqid": seqid,
                "source": source,
                "type": "intergenic_region",
                "start": 1,
                "end": first_start - 1,
                "score": ".",
                "strand": "+",
                "phase": ".",
                "attributes": "name_up="
                f"{df.iat[-1, NAME_RIGHT_IDX] if circular else 'region_start'};"
                "source_up="
                f"{df.iat[-1, SOURCE_RIGHT_IDX] if circular else 'region_start'};"
                f"name_dw={df.iat[0, NAME_LEFT_IDX]};"
                f"source_dw={df.iat[0, SOURCE_LEFT_IDX]};",
            }
        )
    # [---...][== igr ==]||...
    # if the last end is less than the region size, we have an intergenic
    # region after the last feature, so we create an intergenic region from
    # last_end + 1 to region_size.
    if last_end < region_size:
        regions.append(
            {
                "seqid": seqid,
                "source": source,
                "type": "intergenic_region",
                "start": last_end + 1,
                "end": region_size,
                "score": ".",
                "strand": "+",
                "phase": ".",
                "attributes": f"name_up={df.iat[-1, NAME_RIGHT_IDX]};"
                f"source_up={df.iat[-1, SOURCE_RIGHT_IDX]};"
                "name_dw="
                f"{df.iat[0, NAME_LEFT_IDX] if circular else 'region_end'};"
                "source_dw="
                f"{df.iat[0, SOURCE_LEFT_IDX] if circular else 'region_end'};",
            }
        )
    return regions


def _create_intergenic(
    log: log_setup.TempLogger,
    df: pd.DataFrame,
    seqid: str,
    source: str,
) -> pd.DataFrame:

    # we change the original df, so as to not have to copy it
    df["start"] = df["start"] - 1
    df["end"] = df["end"] + 1

    # we create a shifted version with only the data we'll use
    shifted = df[["end", "name_right", "source_right"]].shift(1)
    valid_gaps = shifted["end"] <= df["start"]

    if not valid_gaps.any():
        return pd.DataFrame()

    attributes = [
        f"name_up={nu};source_up={su};name_dw={nd};source_dw={sd};"
        for nu, su, nd, sd in zip(
            shifted["name_right"][valid_gaps],
            shifted["source_right"][valid_gaps],
            df["name_left"][valid_gaps],
            df["source_left"][valid_gaps],
        )
    ]

    return pd.DataFrame(
        {
            "seqid": seqid,
            "source": source,
            "type": "intergenic_region",
            "start": shifted["end"][valid_gaps].values.astype("int64"),
            "end": df["start"][valid_gaps].values,
            "score": ".",
            "strand": "+",
            "phase": ".",
            "attributes": attributes,
        }
    )


def _split_attributes(
    log: log_setup.TempLogger,
    df: pd.DataFrame,
) -> pd.DataFrame:
    df["name_left"] = (
        df["attributes"]
        .str.extract(_RE_name, expand=False)  # type: ignore[call-overload]
        .fillna(
            df["attributes"].str.extract(_RE_name_left, expand=False)  # type: ignore[call-overload]
        )
    )
    df["source_left"] = (
        df["attributes"]
        .str.extract(_RE_source, expand=False)  # type: ignore[call-overload]
        .fillna(
            df["attributes"].str.extract(_RE_source_left, expand=False)  # type: ignore[call-overload]
        )
    )

    df["name_right"] = (
        df["attributes"]
        .str.extract(_RE_name, expand=False)  # type: ignore[call-overload]
        .fillna(
            df["attributes"].str.extract(_RE_name_right, expand=False)  # type: ignore[call-overload]
        )
    )
    df["source_right"] = (
        df["attributes"]
        .str.extract(_RE_source, expand=False)  # type: ignore[call-overload]
        .fillna(
            df["attributes"].str.extract(_RE_source_right, expand=False)  # type: ignore[call-overload]
        )
    )
    return df


def extract_multiple(
    log: log_setup.GDTLogger,
    tsv_path: Path,
    an_column: str = "AN",
    workers: int = 0,
    gff_in_ext: str = ".gff3",
    gff_in_suffix: str = "_clean",
    gff_out_ext: str = ".gff3",
    gff_out_suffix: str = "_intergenic",
    add_region: bool = False,
    overwrite: bool = False,
) -> None:
    tsv = pd.read_csv(tsv_path, sep="\t")

    gff_in_builder = gff3_utils.PathBuilder(gff_in_ext).use_folder_builder(
        tsv_path.parent, gff_in_suffix
    )
    gff_out_builder = gff3_utils.PathBuilder(gff_out_ext).use_folder_builder(
        tsv_path.parent, gff_out_suffix
    )

    gff3_utils.check_tsv(
        log, tsv, gff_in_builder, gff_out_builder, overwrite, an_column
    )

    log.info(f"Starting processing {tsv.shape[0]} ANs with {workers} workers...")
    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        tasks = [
            executor.submit(
                extract_intergenic_regions,
                log.spawn_buffer(),
                gff_in_builder.build(an),
                gff_out_builder.build(an),
                add_region,
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
