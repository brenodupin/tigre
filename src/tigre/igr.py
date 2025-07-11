# -*- coding: utf-8 -*-
"""Module to extract intergenic regions from GFF3 files."""

import re
from pathlib import Path
from typing import Union

import pandas as pd

from . import gff3_utils, log_setup

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
    log: log_setup.GDTLogger,
    gff_in: Path,
    gff_out: Path,
    keep_orfs: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Extract intergenic regions from a GFF3 file and save to a new file.

    Args:
        log (log_setup.GDTLogger): Logger instance for logging.
        gff_in (Path): Input GFF3 file path.
        gff_out (Path): Output GFF3 file path.
        keep_orfs (bool): Whether to keep ORF entries in the output.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: DataFrames containing the
                                           original GFF3 data and the
                                           intergenic regions.

    """
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

    if not keep_orfs:
        df = gff3_utils.filter_orfs(df)

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
        log.warning("No intergenic regions found. Returning empty DataFrame.")
        return df, pd.DataFrame()

    if boundary:
        df_ig = (
            pd.concat([df_ig, pd.DataFrame(boundary)], ignore_index=True)
            .sort_values(by=["start", "end"])
            .reset_index(drop=True)
        )

    start = 2 if "intergenic_region_merged" in df_ig["type"].values else 1
    df_ig["attributes"] = [
        f"ID={i};{attr}" for i, attr in enumerate(df_ig["attributes"], start)
    ]

    with open(gff_out, "w") as f:
        f.write("\n".join(header) + "\n")
        df_ig.to_csv(f, sep="\t", index=False, header=False)

    return df, df_ig


def _solve_boundaries(
    log: log_setup.GDTLogger,
    df: pd.DataFrame,
    region_size: int,
    seqid: str,
    source: str,
    circular: bool,
) -> list[dict[str, Union[str, int]]]:

    first_start = df.at[0, "start"]
    last_end = df.iat[-1, END_IDX]
    log.debug(f"First start: {first_start}, Last end: {last_end}")

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
                f"{df.iat[-1, NAME_LEFT_IDX] if circular else 'region_start'};"
                "source_up="
                f"{df.iat[-1, SOURCE_LEFT_IDX] if circular else 'region_start'};"
                f"name_dw={df.iat[0, NAME_LEFT_IDX]};"
                f"source_dw={df.iat[0, SOURCE_LEFT_IDX]};",
            }
        )
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
    log: log_setup.GDTLogger,
    df: pd.DataFrame,
    seqid: str,
    source: str,
) -> pd.DataFrame:

    df["start"] = df["start"] - 1
    df["end"] = df["end"] + 1

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


def _split_attributes(log: log_setup.GDTLogger, df: pd.DataFrame) -> pd.DataFrame:
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
