# -*- coding: utf-8 -*-
"""Module to solve overlaps in GFF3 files."""

import concurrent.futures as cf
import re
import traceback
from pathlib import Path
from typing import Callable, TypeAlias

import polars as pl

from . import gff3_utils, log_setup

species_url = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="

HEADER_BASE: str = (
    "##gff-version 3\n" "#!gff-spec-version 1.26\n" "#!processor TIGRE clean.py\n"
)

_RS_name = r"name=([^;]+);"
_RS_source = r"source=([^;]+);"

_RS_name_up = r"name_up=([^;]+);"
_RS_source_up = r"source_up=([^;]+);"

_RS_name_dw = r"name_dw=([^;]+);"
_RS_source_dw = r"source_dw=([^;]+);"

_RS_name_left = r"name_left=([^;]+);"
_RS_source_left = r"source_left=([^;]+);"

_RS_name_right = r"name_right=([^;]+);"
_RS_source_right = r"source_right=([^;]+);"

# backwards compatibility
_RE_name = re.compile(_RS_name)
_RE_source = re.compile(_RS_source)

_RE_name_up = re.compile(_RS_name_up)
_RE_source_up = re.compile(_RS_source_up)

_RE_name_dw = re.compile(_RS_name_dw)
_RE_source_dw = re.compile(_RS_source_dw)

_RE_name_left = re.compile(_RS_name_left)
_RE_source_left = re.compile(_RS_source_left)

_RE_name_right = re.compile(_RS_name_right)
_RE_source_right = re.compile(_RS_source_right)


def adjust_coords(
    log: log_setup.TempLogger,
    df: pl.DataFrame,
    region_end: int,
) -> pl.DataFrame:
    """Adjust features with 'end' and 'start' bigger than the region."""
    mask = (df["end"] > region_end) & (df["start"] > region_end)

    if not mask.any():
        return df

    log.warning(f"Feature with 'end' and 'start' bigger than region: {region_end}")

    for row in df.filter(mask).iter_rows(named=True):
        log.warning(
            f" s: {row['start']} e: {row['end']} | t: {row['type']} | "
            f"nl: {row['name_left']} | nr: {row['name_right']}"
        )
        log.debug(f" Adjusted s: {row['start'] - region_end} e: {row['end'] - region_end}")

    return df.with_columns(
        pl.when(mask)
        .then(pl.col("start") - region_end)
        .otherwise(pl.col("start"))
        .alias("start"),
        pl.when(mask)
        .then(pl.col("end") - region_end)
        .otherwise(pl.col("end"))
        .alias("end"),
    ).sort(["start", "type"], maintain_order=True)


def bigger_than_region_solver(
    log: log_setup.TempLogger,
    df: pl.DataFrame,
    region_end: int,
) -> pl.DataFrame:
    """Handle features with 'end' bigger than the region."""
    mask = (df["end"] > region_end) & (df["start"] <= region_end)
    subset = df.filter(mask)

    if subset.is_empty():
        return df

    log.warning(f"Feature with 'end' bigger than region: {region_end}")

    for row in subset.iter_rows(named=True):
        log.warning(
            f" s: {row['start']} e: {row['end']} | t: {row['type']} | "
            f"nl: {row['name_left']} | nr: {row['name_right']}"
        )
        if row["end"] // region_end > 1:
            log.error(f" Feature is bigger than region! g:{row['end']} r:{region_end}")
        log.debug(
            f"Processing {row['type']}: nl: {row['name_left']} sr: {row['source_left']} |"
            f" nr: {row['name_right']} sl: {row['source_right']}"
        )
        log.debug(
            f" 1 s: {row['start']} | e: {region_end} | nl: {row['name_left']} |"
            f" nr: {row['name_right']}"
        )
        log.debug(
            f" 2 s: 1 | e: {row['end'] - region_end} | nl: {row['name_left']} |"
            f" nr: {row['name_right']}"
        )

    overlap = subset.with_columns(
        pl.lit(1, dtype=df.schema["start"]).alias("start"),
        (pl.col("end") - region_end).alias("end"),
        (pl.col("type") + "_fragment").alias("type"),
    )

    df = df.with_columns(
        pl.col("end").clip(upper_bound=region_end),
        pl.when(mask)
        .then(pl.col("type") + "_fragment")
        .otherwise(pl.col("type"))
        .alias("type"),
    )

    return pl.concat([df, overlap])


def overlaps_chooser(
    log: log_setup.TempLogger,
    overlaps: pl.DataFrame,
) -> str:
    """Choose the best overlap from a list of overlaps."""
    overlaps = overlaps.sort("start", maintain_order=True)

    min_start = overlaps["start"].min()
    max_end = overlaps["end"].max()
    left = overlaps.filter(pl.col("start") == min_start)

    if len(left) > 1:
        log.trace(f" More than one feature with the lowest start: {min_start!s}")
        left_max_end = left["end"].max()
        left = left.filter(pl.col("end") == left_max_end)
        if len(left) > 1 and left_max_end == max_end:
            row = left.row(0, named=True)
            log.trace(
                " Features with coinciding coordinates, choosing the first: "
                f"nl: {row['name_left']} nr: {row['name_right']} | "
                f"att: {row['attributes']}"
            )
            return (
                f"name_left={row['name_left']};"
                f"source_left={row['source_left']};"
                f"name_right={row['name_right']};"
                f"source_right={row['source_right']};"
            )

    right = overlaps.filter(pl.col("end") == max_end)
    if len(right) > 1:
        log.trace(f" More than one feature with the highest end: {max_end!s}")
        right = right.filter(pl.col("start") == right["start"].min())

    left_row = left.row(0, named=True)
    right_row = right.row(0, named=True)
    return (
        f"name_left={left_row['name_left']};"
        f"source_left={left_row['source_left']};"
        f"name_right={right_row['name_right']};"
        f"source_right={right_row['source_right']};"
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
    df: pl.DataFrame,
) -> pl.DataFrame:
    """Solve overlaps in the GFF3 file."""
    df = df.with_columns(
        (pl.col("start") > pl.col("end").shift(1).cum_max())
        .fill_null(False)
        .cast(pl.Int32)
        .cum_sum()
        .alias("_group")
    )

    group_sizes = df["_group"].value_counts().filter(pl.col("count") == 1)["_group"]
    solo = df.filter(pl.col("_group").is_in(group_sizes))
    overlapping = df.filter(~pl.col("_group").is_in(group_sizes))

    if overlapping.is_empty():
        return solo.drop("_group")

    merged_rows = []
    for group_df in overlapping.partition_by("_group", maintain_order=True):
        group_df = group_df.sort("start", maintain_order=True)
        end_region = group_df["end"].max()
        row = group_df.row(0, named=True)

        log.debug(
            f"overlap found, region: {row['start']} - {end_region!s} | "
            f"length: {end_region - row['start'] + 1}"
        )
        for r in group_df.iter_rows(named=True):
            log.trace(
                f" {r['start']} | {r['end']} | {r['type']} | nl: {r['name_left']} |"
                f" nr: {r['name_right']} | sl: {r['source_left']} | sr: {r['source_right']}"
            )

        attributes = overlaps_chooser(log, group_df)
        log.trace(f"Result: {attributes}")
        merged_rows.append(
            {
                **row,
                "end": end_region,
                "type": "overlapping_feature_set",
                "strand": "+",
                "attributes": attributes,
            }
        )

    merged = pl.DataFrame(merged_rows, schema=df.schema)
    return (
        pl.concat([solo, merged])
        .drop("_group")
        .sort(["start", "type"], maintain_order=True)
    )


_NamesType: TypeAlias = (
    "Callable[[list[str], log_setup.TempLogger], tuple[list[str], log_setup.TempLogger]]"
)


def get_names(
    gene_ids: list[str],
    log: log_setup.TempLogger,
) -> "tuple[list[str], log_setup.TempLogger]":
    """Get names from a list of gene IDs."""
    return gene_ids, log


def format_df(
    df: pl.DataFrame, log: log_setup.TempLogger, names_func: _NamesType
) -> pl.DataFrame:
    """Gather and format the DataFrame to ensure correct data types."""
    df = (
        df.with_row_index("_idx")
        .with_columns(
            pl.col("attributes").str.extract(_RS_name_up, 1).alias("_nm_up"),
            pl.col("attributes").str.extract(_RS_name_left, 1).alias("_nm_left"),
            pl.col("attributes").str.extract(_RS_name_dw, 1).alias("_nm_dw"),
            pl.col("attributes").str.extract(_RS_name_right, 1).alias("_nm_right"),
            pl.col("attributes").str.extract(_RS_name, 1).alias("_nm_general"),
            pl.col("attributes").str.extract(_RS_source_up, 1).alias("_src_up"),
            pl.col("attributes").str.extract(_RS_source_left, 1).alias("_src_left"),
            pl.col("attributes").str.extract(_RS_source_dw, 1).alias("_src_dw"),
            pl.col("attributes").str.extract(_RS_source_right, 1).alias("_src_right"),
            pl.col("attributes").str.extract(_RS_source, 1).alias("_src_general"),
        )
        .with_columns(
            pl.coalesce(["_nm_up", "_nm_left", "_nm_general"]).alias("name_left"),
            pl.coalesce(["_src_up", "_src_left", "_src_general"]).alias("source_left"),
            pl.coalesce(["_nm_dw", "_nm_right", "_nm_general"]).alias("name_right"),
            pl.coalesce(["_src_dw", "_src_right", "_src_general"]).alias("source_right"),
        )
        .drop(
            [
                "_nm_up",
                "_nm_left",
                "_nm_dw",
                "_nm_right",
                "_nm_general",
                "_src_up",
                "_src_left",
                "_src_dw",
                "_src_right",
                "_src_general",
            ]
        )
    )

    og_row_mask = (
        df["name_left"].is_null()
        | df["source_left"].is_null()
        | df["name_right"].is_null()
        | df["source_right"].is_null()
    )

    if og_row_mask.any():
        subset = df.filter(og_row_mask)
        # think about improving this part, as the input grows, this will become a bottleneck
        gene_labels, log = names_func(subset["gene_id"].to_list(), log)

        updated_subset = (
            subset.with_columns(
                pl.Series("_gene_labels", gene_labels),
                pl.concat_str(
                    ["seqid", "type", "start", "end", "strand", "gene_id"],
                    separator="|",
                ).alias("_source_fallback"),
            )
            .with_columns(
                pl.coalesce(["name_left", "_gene_labels"]).alias("name_left"),
                pl.coalesce(["name_right", "_gene_labels"]).alias("name_right"),
                pl.coalesce(["source_left", "_source_fallback"]).alias("source_left"),
                pl.coalesce(["source_right", "_source_fallback"]).alias("source_right"),
            )
            .with_columns(
                pl.concat_str(
                    [
                        pl.lit("name="),
                        pl.col("name_left"),
                        pl.lit(";source="),
                        pl.col("source_left"),
                        pl.lit(";"),
                    ]
                ).alias("attributes")
            )
            .drop(["_gene_labels", "_source_fallback"])
        )

        df = df.update(updated_subset, on="_idx")

    return df.drop("_idx")


def clean_an(
    log: log_setup.TempLogger,
    gff_in: Path,
    gff_out: Path,
    names_func: _NamesType,
    query_expr: "pl.Expr",
    keep_orfs: bool,
    ext_filter: bool = False,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Clean a GFF3 file by removing unnecessary features and attributes.

    Args:
        log: Logger instance
        gff_in: Path to the input GFF3 file
        gff_out: Path to the output GFF3 file
        names_func: Function to clean attributes in each row
        query_expr: Query string to filter the GFF3 file
        keep_orfs: Whether to keep ORFs in the GFF3 file
        ext_filter: Whether to use extended filtering for ORFs

    """
    try:
        log.debug(f"[{gff_in.name}] -- Start --")

        log.trace(f" Loading GFF3 file: {gff_in}")
        df = gff3_utils.load_gff3(
            gff_in,
            usecols=gff3_utils.GFF3_COLUMNS,
            query_expr=query_expr,
            return_polars=True,
        )

        region_df = df.head(1)
        df = df.slice(1)

        if not keep_orfs:
            df = gff3_utils.filter_orfs(df, extended=ext_filter, return_polars=True)

        an = region_df.item(0, "seqid")
        region_end = region_df.item(0, "end")
        header = create_header(an, region_df.item(0, "attributes"), region_end)

        df = df.filter(pl.col("type") != "region")
        df = df.with_columns(
            pl.col("attributes").str.extract(gff3_utils._RS_ID, 1).alias("gene_id")
        )

        df = format_df(df, log, names_func)
        df = adjust_coords(log, df, region_end)
        df = bigger_than_region_solver(log, df, region_end)
        df = df.sort(["start", "type"], maintain_order=True)

        log.trace(" Overlaps solver...")
        df = overlap_solver(log, df)

        df = df.drop(
            [
                col
                for col in [
                    "gene_id",
                    "name_left",
                    "source_left",
                    "name_right",
                    "source_right",
                ]
                if col in df.columns
            ]
        )
        df = pl.concat([region_df, df])

        # add header to the file
        with open(gff_out, "w", encoding="utf-8") as file_handler:
            file_handler.write(header)
            df.write_csv(file_handler, separator="\t", include_header=False)

        log.debug(f"[{an}] -- End --")
        return True, an, log.get_records()

    except Exception:
        an_error: str = an if "an" in locals() else gff_in.name
        error_msg = traceback.format_exc()
        log.error(f"Error in {an_error}:\n{error_msg}")
        return False, an_error, log.get_records()


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
    gff3_utils._ensure_spawn(log)  # Ensure spawn method is set for multiprocessing
    tsv = pl.read_csv(tsv_path, separator="\t")
    try:
        query_expr = gff3_utils._parse_query_to_polars(query_string)
        log.info(f"Using query expression: {query_expr}")
    except Exception as e:
        log.error(f"Error parsing query string '{query_string}': {e}")
        raise

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
                query_expr,
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
