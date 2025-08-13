# -*- coding: utf-8 -*-
"""Module to solve overlaps in GFF3 files."""

import concurrent.futures as cf
import sys
import traceback
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from . import clean, clean_gdt_server, gff3_utils, log_setup

if TYPE_CHECKING:
    import gdt  # type: ignore[import-not-found]

    from .cli import _Args


def clean_attr_gdt(
    row: "pd.Series",
    gdict: "gdt.GeneDictionary",
    log: log_setup.TempLogger,
) -> str:
    """Clean attributes using GDT."""
    try:
        label = gdict[row.gene_id].label
    except KeyError:
        log.error(
            f"Gene ID {row.gene_id} not found in GDT dictionary. Using default format."
        )
        label = row.gene_id
    return (
        f"name={label};source={row.seqid}|{row.type}|{row.start}|{row.end}|"
        f"{row.strand}|{row.gene_id};"
    )


def clean_an_gdt(
    log: log_setup.TempLogger,
    gff_in: Path,
    gff_out: Path,
    gdict: "gdt.GeneDict",
    query_string: str,
    keep_orfs: bool,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Clean a GFF3 file by removing unnecessary features and attributes with a gdict.

    Args:
        log: Logger instance
        gff_in: Path to the input GFF3 file
        gff_out: Path to the output GFF3 file
        gdict: GDT GeneDict instance
        query_string: Query string to filter the GFF3 file
        keep_orfs: Whether to keep ORFs in the GFF3 file

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
            df = gff3_utils.filter_orfs(df)
        an = df.at[0, "seqid"]
        header = clean.create_header(an, df.at[0, "attributes"], df.at[0, "end"])
        df["gene_id"] = df["attributes"].str.extract(gff3_utils._RE_ID, expand=False)  # type: ignore[call-overload]

        if (df["type"] == "region").sum() > 1:
            df = clean.multiple_regions_solver(log, df, an)

        df.loc[1:, "attributes"] = df.loc[1:].apply(
            clean_attr_gdt, log=log, gdict=gdict, axis=1
        )

        if df["end"].idxmax() != 0:
            df = clean.bigger_than_region_solver(log, df, an)

        log.trace(" Sorting DataFrame by start and type...")
        df.at[df[df["type"] == "region"].index[0], "type"] = "_"
        df = df.sort_values(by=["start", "type"], ascending=True, ignore_index=True)
        df.at[0, "type"] = "region"

        log.trace(" Overlaps solver...")
        df = clean.overlap_solver(log, df, an)
        df = df.drop(columns=["gene_id"], errors="ignore")

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


def clean_multiple_gdt(
    log: log_setup.GDTLogger,
    tsv_path: Path,
    workers: int,
    gdict: "gdt.GeneDict",
    gff_in_ext: str = ".gff3",
    gff_in_suffix: str = "",
    gff_out_ext: str = ".gff3",
    gff_out_suffix: str = "_clean",
    an_column: str = "AN",
    query_string: str = gff3_utils.QS_GENE_TRNA_RRNA_REGION,
    keep_orfs: bool = False,
    overwrite: bool = False,
) -> None:
    """Orchestrates the execution of `clean_an` across multiple GFF3 files.

    Args:
        log: Logger instance
        tsv_path: Path to the TSV file containing ANs
        workers: Number of worker processes to use
        gdict: GDT GeneDict instance
        gff_in_ext: Extension for input GFF3 files
        gff_in_suffix: Suffix for input GFF3 files
        gff_out_ext: Extension for output GFF3 files
        gff_out_suffix: Suffix for output GFF3 files
        an_column: Column name in the TSV file that contains ANs
        query_string: Query string to filter the GFF3 file
        keep_orfs: Whether to keep ORFs in the GFF3 file
        overwrite: Whether to overwrite existing output files

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
    with cf.ThreadPoolExecutor(max_workers=workers) as executor:
        tasks = [
            executor.submit(
                clean_an_gdt,
                log.spawn_buffer(),
                gff_in_builder.build(an),
                gff_out_builder.build(an),
                gdict,
                query_string,
                keep_orfs,
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


def clean_gdt_multiple(
    log: log_setup.GDTLogger,
    args: "_Args",
    workers_process: int,
    workers_threading: int,
) -> None:
    """Handle `clean` command with an GDT gdict.

    Args:
        log: Logger instance
        args: Parsed command line arguments
        workers_process: Number of worker processes to use
        workers_threading: Number of worker threads to use

    """
    gdict = load_gdt(log, args.gdt)

    server_mode = False

    if len(gdict) > 200000:
        server_mode = True

    if (not args.no_server) and (args.server or server_mode):
        log.info("Running in server mode.")
        clean_gdt_server.clean_multiple_gdt_server(
            log,
            args.tsv,
            workers_process,
            gdict,
            args.gff_in_ext,
            args.gff_in_suffix,
            args.gff_out_ext,
            args.gff_out_suffix,
            args.an_column,
            args.query_string,
            args.keep_orfs,
            args.overwrite,
        )

    else:
        log.info("Running in normal mode.")
        clean_multiple_gdt(
            log,
            args.tsv,
            workers_threading,
            gdict,
            args.gff_in_ext,
            args.gff_in_suffix,
            args.gff_out_ext,
            args.gff_out_suffix,
            args.an_column,
            args.query_string,
            args.keep_orfs,
            args.overwrite,
        )


def load_gdt(
    log: log_setup.GDTLogger,
    gdict_path: Path,
) -> "gdt.GeneDict":
    """Load GDT Library and read the gdict file."""
    try:
        import gdt
    except ImportError:
        raise SystemExit(
            "GDT package not found. Please install it to use the --gdict option."
        )

    log.debug(f"GDT package version: {gdt.__version__}")
    gdt_path: Path = Path(gdict_path).resolve()
    if not gdt_path.is_file():
        log.error(f"GDT .gdict file not found: {gdt_path}")
        sys.exit(1)

    gdict = gdt.read_gdict(gdt_path, lazy_info=False)
    log.info(f"GDT dictionary loaded from {gdt_path}")
    gdt.log_info(log, gdict)
    return gdict
