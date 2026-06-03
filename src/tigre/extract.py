# -*- coding: utf-8 -*-
"""Module to extract intergenic regions from GFF3 files."""

import concurrent.futures as cf
from pathlib import Path

import polars as pl

from . import clean, gff3_utils, log_setup


def update_header_annotator(header: list[str]) -> list[str]:
    """Update the GFF3 header with the annotator information.

    Args:
        header (list[str]): List of header lines from the GFF3 file.
        annotator (str): Name of the annotator to be added.

    Returns:
        list[str]: Updated list of header lines.

    """
    new_header = []
    for line in header:
        if line == "#!processor TIGRE clean.py":
            new_header.append("#!processor TIGRE igr.py")
        else:
            new_header.append(line)

    return new_header


def extract_intergenic_regions(
    log: log_setup.TempLogger,
    gff_in: Path,
    gff_out: Path,
    add_region: bool = False,
    feature_type: str = "intergenic_region",
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Extract intergenic regions from a GFF3 file and save to a new file.

    Args:
        log (log_setup.GDTLogger): Logger instance for logging.
        gff_in (Path): Input GFF3 file path.
        gff_out (Path): Output GFF3 file path.
        add_region (bool): Whether to add a region line to the output file.
                           Defaults is False.
        feature_type (str): Type of the intergenic region to be used in the output.
                            Defaults to "intergenic_region". If the feature spans the
                            genome boundary in circular genomes, it will be appended
                            with "_merged".


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
                gff_in, usecols=gff3_utils.GFF3_COLUMNS, return_polars=True
            )
        header = update_header_annotator(header)

        # pop region line out of df
        region_df = df.head(1)
        seqid = region_df.item(0, "seqid")
        df = df.slice(1)

        if df.is_empty():
            log.warning(f"No features found in {gff_in}.")
            write_gff3(
                log,
                pl.DataFrame(),
                region_df,
                gff_out,
                header,
                add_region,
            )
            return True, seqid, log.get_records()

        region_size = region_df.item(0, "end")
        source = region_df.item(0, "source")
        circular = "is_circular=true" in region_df.item(0, "attributes").lower()

        log.debug(
            f"AN: {seqid} | Source: {source} | Size: {region_size} | "
            f"Circular: {circular} | Feature type: {feature_type}"
        )

        df = _split_attributes(log, df)

        boundary = _solve_boundaries(
            log,
            df,
            region_size,
            seqid,
            source,
            circular,
            feature_type,
        )

        # this will alter df! by adding +1 to end and -1 to start
        df_ig = _create_intergenic(
            log,
            df,
            seqid,
            source,
            feature_type,
        )

        if df_ig.is_empty() and not boundary:
            log.warning(
                f"No intergenic regions found. Did {gff_in.name} contain any features?"
            )
            write_gff3(
                log,
                pl.DataFrame(),
                region_df,
                gff_out,
                header,
                add_region,
            )
            return True, seqid, log.get_records()

        if boundary:
            df_ig = pl.concat([df_ig, pl.DataFrame(boundary)]).sort(
                ["start", "end"], descending=[False, True]
            )

        start = 2 if df_ig["type"].is_in([f"{feature_type}_merged"]).any() else 1

        df_ig = df_ig.with_columns(
            pl.Series(
                [
                    f"ID={seqid}_{feature_type}_{i};name_up={name_up};source_up={source_up};name_dw={name_dw};source_dw={source_dw};"
                    for i, (name_up, source_up, name_dw, source_dw) in enumerate(
                        zip(
                            df_ig["name_up"],
                            df_ig["source_up"],
                            df_ig["name_dw"],
                            df_ig["source_dw"],
                        ),
                        start=start,
                    )
                ]
            ).alias("attributes")
        ).drop(["name_up", "source_up", "name_dw", "source_dw"])

        write_gff3(
            log,
            df_ig,
            region_df,
            gff_out,
            header,
            add_region,
        )
        return True, seqid, log.get_records()

    except Exception as e:
        an_error = seqid if "seqid" in locals() else gff_in.name
        log.error(f"Error in {an_error}: {e}")
        return False, an_error, log.get_records()


def write_gff3(
    log: log_setup.TempLogger,
    df: pl.DataFrame,
    region: pl.DataFrame,
    gff_out: Path,
    header: list[str],
    add_region: bool = False,
) -> None:
    """Write a DataFrame to a GFF3 file with optional region line."""
    log.trace(f"Writing intergenic regions to {gff_out}.")
    with open(gff_out, "w") as f:
        f.write("\n".join(header) + "\n")
        if add_region:
            df = pl.concat([region, df])
        df.write_csv(f, separator="\t", include_header=False)


def _solve_boundaries(
    log: log_setup.TempLogger,
    df: pl.DataFrame,
    region_size: int,
    seqid: str,
    source: str,
    circular: bool,
    feature_type: str,
) -> list[dict[str, str | int]]:
    # Legend:
    # [---] feature (gene, trna, rrna)
    # [== igr ==] intergenic region
    # || genome boundary, i.e. position where region_size rolls over to 1

    first_start = df.item(0, "start")  # of a feature
    last_end = df.item(-1, "end")  # of a feature
    # log.debug(f"First start: {first_start}, Last end: {last_end}")

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
                "type": f"{feature_type}_merged",
                "start": last_end + 1,
                "end": region_size + first_start - 1,
                "score": ".",
                "strand": "+",
                "phase": ".",
                "name_up": df.item(-1, "name_right"),
                "source_up": df.item(-1, "source_right"),
                "name_dw": df.item(0, "name_left"),
                "source_dw": df.item(0, "source_left"),
            }
        ]

    regions: list[dict[str, str | int]] = []
    # ...||[== igr ==][---...
    # if the first start is greater than 1, we have an intergenic region
    # before the first feature, so we create an intergenic region from 1 to
    # first_start - 1.
    if first_start > 1:
        regions.append(
            {
                "seqid": seqid,
                "source": source,
                "type": feature_type,
                "start": 1,
                "end": first_start - 1,
                "score": ".",
                "strand": "+",
                "phase": ".",
                "name_up": df.item(-1, "name_right") if circular else "region_start",
                "source_up": df.item(-1, "source_right") if circular else "region_start",
                "name_dw": df.item(0, "name_left"),
                "source_dw": df.item(0, "source_left"),
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
                "type": feature_type,
                "start": last_end + 1,
                "end": region_size,
                "score": ".",
                "strand": "+",
                "phase": ".",
                "name_up": df.item(-1, "name_right"),
                "source_up": df.item(-1, "source_right"),
                "name_dw": df.item(0, "name_left") if circular else "region_end",
                "source_dw": df.item(0, "source_left") if circular else "region_end",
            }
        )
    return regions


def _create_intergenic(
    log: log_setup.TempLogger,
    df: pl.DataFrame,
    seqid: str,
    source: str,
    feature_type: str,
) -> pl.DataFrame:

    df = df.with_columns(
        [
            (pl.col("start") - 1).alias("start"),
            (pl.col("end") + 1).alias("end"),
        ]
    )

    shifted = df.select(["end", "name_right", "source_right"]).shift(1)
    valid_gaps = shifted["end"] <= df["start"]

    if not valid_gaps.any():
        return pl.DataFrame()

    shifted = shifted.filter(valid_gaps)
    df = df.filter(valid_gaps)

    return pl.DataFrame(
        {
            "seqid": seqid,
            "source": source,
            "type": feature_type,
            "start": shifted["end"].cast(pl.Int64),
            "end": df["start"],
            "score": ".",
            "strand": "+",
            "phase": ".",
            "name_up": shifted["name_right"],
            "source_up": shifted["source_right"],
            "name_dw": df["name_left"],
            "source_dw": df["source_left"],
        }
    )


def _split_attributes(
    log: log_setup.TempLogger,
    df: pl.DataFrame,
) -> pl.DataFrame:
    """Extract name and source attributes with fallback to general values."""
    return (
        df.with_columns(
            [
                # Group 1 extracts the capture group ([^;]+)
                pl.col("attributes").str.extract(clean._RS_name, 1).alias("name_general"),
                pl.col("attributes")
                .str.extract(clean._RS_source, 1)
                .alias("source_general"),
                pl.col("attributes")
                .str.extract(clean._RS_name_left, 1)
                .alias("_name_left"),
                pl.col("attributes")
                .str.extract(clean._RS_name_right, 1)
                .alias("_name_right"),
                pl.col("attributes")
                .str.extract(clean._RS_source_left, 1)
                .alias("_source_left"),
                pl.col("attributes")
                .str.extract(clean._RS_source_right, 1)
                .alias("_source_right"),
            ]
        )
        .with_columns(
            [
                pl.coalesce(["_name_left", "name_general"]).alias("name_left"),
                pl.coalesce(["_name_right", "name_general"]).alias("name_right"),
                pl.coalesce(["_source_left", "source_general"]).alias("source_left"),
                pl.coalesce(["_source_right", "source_general"]).alias("source_right"),
            ]
        )
        .drop(
            [
                "name_general",
                "source_general",
                "_name_left",
                "_name_right",
                "_source_left",
                "_source_right",
            ]
        )
    )


def extract_multiple(
    log: log_setup.GDTLogger,
    tsv_path: Path,
    workers: int,
    an_column: str = "AN",
    gff_in_ext: str = ".gff3",
    gff_in_suffix: str = "_clean",
    gff_out_ext: str = ".gff3",
    gff_out_suffix: str = "_intergenic",
    add_region: bool = False,
    overwrite: bool = False,
    feature_type: str = "intergenic_region",
) -> None:
    """Extract intergenic regions from multiple GFF3 files listed in a TSV file.

    Args:
        log (log_setup.GDTLogger): Logger instance for logging.
        tsv_path (Path): Path to the TSV file containing ANs.
        workers (int): Number of worker threads/processes to use.
        an_column (str): Column name in the TSV file that contains the ANs.
                         Defaults to "AN".
        gff_in_ext (str): Extension for input GFF3 files. Defaults to ".gff3".
        gff_in_suffix (str): Suffix for input GFF3 files. Defaults to "_clean".
        gff_out_ext (str): Extension for output GFF3 files. Defaults to ".gff3".
        gff_out_suffix (str): Suffix for output GFF3 files. Defaults to "_intergenic".
        add_region (bool): Whether to add a region line to the output file.
                           Defaults is False.
        overwrite (bool): Whether to overwrite existing output files. Defaults is False.
        feature_type (str): Type of the intergenic region to be used in the output.
                            Defaults to "intergenic_region". If the feature spans the
                            genome boundary in circular genomes, it will be appended
                            with "_merged".

    """
    gff3_utils._ensure_spawn(log)
    tsv = pl.read_csv(tsv_path, separator="\t")

    gff_in_builder = gff3_utils.PathBuilder(gff_in_ext).use_folder_builder(
        tsv_path.parent, gff_in_suffix
    )
    gff_out_builder = gff3_utils.PathBuilder(gff_out_ext).use_folder_builder(
        tsv_path.parent, gff_out_suffix
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
                extract_intergenic_regions,
                log.spawn_buffer(),
                gff_in_builder.build(an),
                gff_out_builder.build(an),
                add_region,
                feature_type,
            )
            for an in tsv[an_column]
        ]
        log.info(f"Created {len(tasks)} tasks, starting log collection...")
        for future in cf.as_completed(tasks):
            success, an, records = future.result()
            if not success:
                log.error(f"[{an}] Failed to process.")

            for record in records:
                log.log(record[0], record[1])

        log.info("All processed.")
