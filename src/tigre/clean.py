# -*- coding: utf-8 -*-
"""Module to solve overlaps in GFF3 files."""

import concurrent.futures as cf
import re
import traceback
from pathlib import Path
from typing import Callable, TypeAlias, cast

import pandas as pd

from . import gff3_utils, log_setup

species_url = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="

HEADER_BASE: str = (
    "##gff-version 3\n" "#!gff-spec-version 1.26\n" "#!processor TIGRE clean.py\n"
)

_RE_name = re.compile(r"name=([^;]+);")
_RE_source = re.compile(r"source=([^;]+);")

_RE_name_up = re.compile(r"name_up=([^;]+);")
_RE_source_up = re.compile(r"source_up=([^;]+);")

_RE_name_dw = re.compile(r"name_dw=([^;]+);")
_RE_source_dw = re.compile(r"source_dw=([^;]+);")

_RE_name_left = re.compile(r"name_left=([^;]+);")
_RE_source_left = re.compile(r"source_left=([^;]+);")

_RE_name_right = re.compile(r"name_right=([^;]+);")
_RE_source_right = re.compile(r"source_right=([^;]+);")

START_IDX = 3
END_IDX = 4

NAME_LEFT_IDX = 10  # df['name_left']
SOURCE_LEFT_IDX = 11  # df['source_left']

NAME_RIGHT_IDX = 12  # df['name_right']
SOURCE_RIGHT_IDX = 13  # df['source_right']


def _get_column_indices(df: pd.DataFrame) -> None:
    """Update global indices for name and source columns."""
    global START_IDX, END_IDX, NAME_LEFT_IDX
    global NAME_RIGHT_IDX, SOURCE_LEFT_IDX, SOURCE_RIGHT_IDX

    assert df.columns.is_unique, "DataFrame columns must be unique"

    # cast to int to satisfy mypy, since is_unique guarantees the
    # column exists and is unique
    START_IDX = cast(int, df.columns.get_loc("start"))
    END_IDX = cast(int, df.columns.get_loc("end"))
    NAME_LEFT_IDX = cast(int, df.columns.get_loc("name_left"))
    NAME_RIGHT_IDX = cast(int, df.columns.get_loc("name_right"))
    SOURCE_LEFT_IDX = cast(int, df.columns.get_loc("source_left"))
    SOURCE_RIGHT_IDX = cast(int, df.columns.get_loc("source_right"))


def bigger_than_region_solver(
    log: log_setup.TempLogger,
    df: pd.DataFrame,
    region_end: int,
) -> pd.DataFrame:
    """Handle features with 'end' bigger than the region."""
    log.warning(f"Feature with 'end' bigger than region: {region_end}")
    for t in df[df["end"] > region_end].itertuples():
        log.warning(
            f" s: {t.start} e: {t.end} | t: {t.type} | nl: {t.name_left} |"
            f" nr: {t.name_right}"
        )

    new_rows = []
    for index, row in df[df["end"] > region_end].iterrows():
        log.debug(
            f"Processing {row.type}: nl: {row.name_left} sr: {row.source_left} |"
            f" nr: {row.name_right} sl: {row.source_right}"
        )

        if row.end // region_end > 1:
            log.error(f" Feature is bigger than region! g:{row.end} r:{region_end}")

        # Handle overlap row
        overlap = row.end - region_end
        df.at[index, "end"] = region_end
        df.at[index, "type"] = row.type + "_fragment"

        # Create new row
        new_row = row.copy()
        new_row.start = 1
        new_row.end = overlap
        new_row.type = row.type + "_fragment"
        new_rows.append(new_row)

        log.debug(
            f" 1 s: {row.start} | e: {region_end} | nl: {row.name_left} |"
            f" nr: {row.name_right}"
        )
        log.debug(
            f" 2 s: {new_row.start} | e: {new_row.end} | nl: {new_row.name_left} |"
            f" nr: {new_row.name_right}"
        )

    return pd.concat([df, pd.DataFrame(new_rows)])


def overlaps_chooser(
    log: log_setup.TempLogger,
    overlaps_in: list[pd.Series | pd.DataFrame],
) -> str:
    """Choose the best overlap from a list of overlaps."""
    overlaps = pd.DataFrame(overlaps_in).sort_values("start", ignore_index=True)

    left = overlaps[overlaps.start == overlaps.start.min()]
    if left.shape[0] > 1:
        log.trace(
            f" More than one feature with the lowest start: {overlaps.start.min()}"
        )
        left = left[left.end == left.end.max()]
        if left.shape[0] > 1 and left.end.max() == overlaps.end.max():
            log.trace(
                " Features with coinciding coordinates, choosing the first: "
                f"nl: {left.iat[0, NAME_LEFT_IDX]} nr: {left.iat[0, NAME_RIGHT_IDX]} | "
                f"att: {left.iat[0, 8]}"
            )

            return (
                f"name_left={left.iat[0, NAME_LEFT_IDX]};"
                f"source_left={left.iat[0, SOURCE_LEFT_IDX]};"
                f"name_right={left.iat[0, NAME_RIGHT_IDX]};"
                f"source_right={left.iat[0, SOURCE_RIGHT_IDX]};"
            )

    right = overlaps[overlaps.end == overlaps.end.max()]
    if right.shape[0] > 1:
        log.trace(f" More than one feature with the highest end: {overlaps.end.max()}")
        right = right[right.start == right.start.min()]

    return (
        f"name_left={left.iat[0, NAME_LEFT_IDX]};"
        f"source_left={left.iat[0, SOURCE_LEFT_IDX]};"
        f"name_right={right.iat[0, NAME_RIGHT_IDX]};"
        f"source_right={right.iat[0, SOURCE_RIGHT_IDX]};"
    )


def create_header(
    an: str,
    attributes: str,
    end: int,
) -> str:
    """Create a GFF3 header with the given attributes."""
    taxon_match = gff3_utils._RE_region_taxon.search(attributes)
    species = species_url + taxon_match.group(1) if taxon_match else "taxon_not_found"
    return HEADER_BASE + f"##sequence-region {an} 1 {end}\n##species {species}\n"


def overlap_solver(
    log: log_setup.TempLogger,
    df: pd.DataFrame,
    df_region: pd.DataFrame,
) -> pd.DataFrame:
    """Solve overlaps in the GFF3 file."""
    next_index = 0
    rows = []
    for index, row in df.iterrows():
        # since both index and next_index are int, we cant use 'is' here
        # due to python int caching values between -5 and 256 (implementation detail)
        # dont ask how i know that
        if index != next_index:
            continue

        end_region = row.end
        next_index += 1
        overlaps = []

        # peek next rows to check for overlaps
        while (next_index < len(df)) and (end_region >= df.iat[next_index, START_IDX]):
            # overlap found, append to list and update next_index
            # so the while will check the next row
            overlaps.append(df.loc[next_index])

            if end_region < df.iat[next_index, END_IDX]:
                end_region = df.iat[next_index, END_IDX]

            next_index += 1

        if overlaps:
            overlaps.insert(0, row.copy())  # insert current row at index 0
            log.debug(
                f"overlap found, region: {row.start} - {end_region} | "
                f"length: {end_region - row.start + 1}"
            )
            for r in overlaps:
                log.trace(
                    f" {r.start} | {r.end} | {r.type} | nl: {r.name_left} | "
                    f" nr: {r.name_right} | sl: {r.source_left} | sr: {r.source_right}"
                )

            # modify the current row with overlap results
            row.end = end_region
            row.type = "overlapping_feature_set"
            row.strand = "+"
            row.attributes = overlaps_chooser(log, overlaps)
            log.trace(f"Result: {row.attributes}")

        # append the current row, be it modified (and therefore overlap) or not
        rows.append(row)

    df_region = pd.concat([df_region, pd.DataFrame(rows)], ignore_index=True)

    return df_region


_NamesType: TypeAlias = "Callable[[pd.DataFrame, log_setup.TempLogger], pd.Series[str]]"


def get_names(subset: pd.DataFrame, log: log_setup.TempLogger) -> "pd.Series[str]":
    """Get names from a subset of the DataFrame."""
    return subset["gene_id"]


def format_df(
    df: pd.DataFrame, log: log_setup.TempLogger, names_func: _NamesType
) -> pd.DataFrame:
    """Gather and format the DataFrame to ensure correct data types."""
    attrs = df["attributes"].str

    name_up = attrs.extract(_RE_name_up, expand=False).astype("string")  # type: ignore[call-overload]
    name_left = attrs.extract(_RE_name_left, expand=False).astype("string")  # type: ignore[call-overload]
    name_right = attrs.extract(_RE_name_right, expand=False).astype("string")  # type: ignore[call-overload]
    name_dw = attrs.extract(_RE_name_dw, expand=False).astype("string")  # type: ignore[call-overload]
    name_general = attrs.extract(_RE_name, expand=False).astype("string")  # type: ignore[call-overload]

    source_up = attrs.extract(_RE_source_up, expand=False).astype("string")  # type: ignore[call-overload]
    source_left = attrs.extract(_RE_source_left, expand=False).astype("string")  # type: ignore[call-overload]
    source_right = attrs.extract(_RE_source_right, expand=False).astype("string")  # type: ignore[call-overload]
    source_dw = attrs.extract(_RE_source_dw, expand=False).astype("string")  # type: ignore[call-overload]
    source_general = attrs.extract(_RE_source, expand=False).astype("string")  # type: ignore[call-overload]

    df["name_left"] = name_up.fillna(name_left).fillna(name_general)
    df["source_left"] = source_up.fillna(source_left).fillna(source_general)

    df["name_right"] = name_dw.fillna(name_right).fillna(name_general)
    df["source_right"] = source_dw.fillna(source_right).fillna(source_general)

    # og_row here means original rows, direct from the gff3 file source
    # meaning they dont have name_left, source_left, name_right or source_right
    # or any other formatting that we do in subsequent steps
    og_row = (
        df["name_left"].isna()
        | df["source_left"].isna()
        | df["name_right"].isna()
        | df["source_right"].isna()
    )

    if og_row.any():
        subset = df[og_row]

        source_fallback = pd.Series(
            [
                f"{seqid}|{typ}|{start}|{end}|{strand}|{gene_id}"
                for seqid, typ, start, end, strand, gene_id in zip(
                    subset["seqid"],
                    subset["type"],
                    subset["start"].astype(str),
                    subset["end"].astype(str),
                    subset["strand"],
                    subset["gene_id"],
                )
            ],
            index=subset.index,
        )

        gene_labels = names_func(subset, log)

        df.loc[og_row, "name_left"] = subset["name_left"].fillna(gene_labels)
        df.loc[og_row, "name_right"] = subset["name_right"].fillna(gene_labels)
        df.loc[og_row, "source_left"] = subset["source_left"].fillna(source_fallback)
        df.loc[og_row, "source_right"] = subset["source_right"].fillna(source_fallback)

        df.loc[og_row, "attributes"] = [
            f"name={name};source={source};"
            for name, source in zip(
                df.loc[og_row, "name_left"],
                df.loc[og_row, "source_left"],
            )
        ]

    return df


def clean_an(
    log: log_setup.TempLogger,
    gff_in: Path,
    gff_out: Path,
    names_func: _NamesType,
    query_string: str,
    keep_orfs: bool,
    ext_filter: bool = False,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Clean a GFF3 file by removing unnecessary features and attributes.

    Args:
        log: Logger instance
        gff_in: Path to the input GFF3 file
        gff_out: Path to the output GFF3 file
        names_func: Function to clean attributes in each row
        query_string: Query string to filter the GFF3 file
        keep_orfs: Whether to keep ORFs in the GFF3 file
        ext_filter: Whether to use extended filtering for ORFs

    """
    try:
        log.debug(f"[{gff_in.name}] -- Start --")

        log.trace(f" Loading GFF3 file: {gff_in}")
        df = gff3_utils.load_gff3(
            gff_in,
            usecols=gff3_utils.GFF3_COLUMNS,
            query_string=query_string,
        )
        if not keep_orfs:
            df = gff3_utils.filter_orfs(df, extended=ext_filter)
        an = df.at[0, "seqid"]
        header = create_header(an, df.at[0, "attributes"], df.at[0, "end"])
        df["gene_id"] = df["attributes"].str.extract(gff3_utils._RE_ID, expand=False)  # type: ignore[call-overload]

        if (df["type"] == "region").sum() > 1:
            log.warning("Multiple regions found:")
            for row in df[df["type"] == "region"].itertuples():
                log.warning(f" s: {row.start} e: {row.end} | {row.attributes}")
            log.warning(" Always choosing the first region and removing the rest.")

        df_region = df.iloc[0:1].copy()
        df = df[df["type"] != "region"]
        df = format_df(df, log, names_func)

        if df["end"].max() > df_region.at[0, "end"]:
            df = bigger_than_region_solver(log, df, df_region.at[0, "end"])

        df = df.sort_values(by=["start", "type"], ascending=True, ignore_index=True)

        log.trace(" Overlaps solver...")
        _get_column_indices(df)
        df = overlap_solver(log, df, df_region)

        df = df.drop(
            columns=[
                "gene_id",
                "name_left",
                "source_left",
                "name_right",
                "source_right",
            ],
            errors="ignore",
        )

        # add header to the file
        with open(gff_out, "w", encoding="utf-8") as file_handler:
            file_handler.write(header)
            df.to_csv(file_handler, sep="\t", header=False, index=False)

        log.debug(f"[{an}] -- End --")
        return True, an, log.get_records()

    except Exception:
        if "an" not in locals():
            error_msg = traceback.format_exc()
            log.error(f"Error in {gff_in}:\n{error_msg}")
            return False, "Not defined", log.get_records()

        error_msg = traceback.format_exc()
        log.error(f"[{an}] Error in {gff_in}:\n{error_msg}")
        log.debug(f"[{an}] -- End --")
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
    query_string: str = gff3_utils.QS_GENE_TRNA_RRNA_REGION,
    keep_orfs: bool = False,
    overwrite: bool = False,
    ext_filter: bool = False,
) -> None:
    """Orchestrates the execution of `clean_an` across multiple GFF3 files.

    Args:
        log: Logger instance
        tsv_path: Path to the TSV file containing ANs
        workers: Number of worker processes to use
        gff_in_ext: Extension for input GFF3 files
        gff_in_suffix: Suffix for input GFF3 files
        gff_out_ext: Extension for output GFF3 files
        gff_out_suffix: Suffix for output GFF3 files
        an_column: Column name in the TSV file that contains ANs
        clean_func: Function to clean attributes in each row
        query_string: Query string to filter the GFF3 file
        keep_orfs: Whether to keep ORFs in the GFF3 file
        overwrite: Whether to overwrite existing output files
        ext_filter: Whether to use extended filtering for ORFs

    """
    tsv = pd.read_csv(tsv_path, sep="\t")

    gff_in_builder = gff3_utils.PathBuilder(gff_in_ext).use_folder_builder(
        tsv_path.parent,
        gff_in_suffix,
    )
    gff_out_builder = gff3_utils.PathBuilder(gff_out_ext).use_folder_builder(
        tsv_path.parent,
        gff_out_suffix,
    )

    gff3_utils.check_files(
        log,
        tsv,
        gff_in_builder,
        an_column,
        should_exist=True,
    )

    if not overwrite:
        gff3_utils.check_files(
            log,
            tsv,
            gff_out_builder,
            an_column,
            should_exist=False,
        )

    log.info(f"Starting processing {tsv.shape[0]} ANs with {workers} workers...")
    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        tasks = [
            executor.submit(
                clean_an,
                log.spawn_buffer(),
                gff_in_builder.build(an),
                gff_out_builder.build(an),
                get_names,
                query_string,
                keep_orfs,
                ext_filter,
            )
            for an in tsv[an_column]
        ]
        log.trace("Tasks submitted, waiting for completion...")

        for future in cf.as_completed(tasks):
            success, an, records = future.result()
            if not success:
                log.error(f"[{an}] Failed to process.")

            for record in records:
                log.log(record[0], record[1])

        log.info("All processed.")
