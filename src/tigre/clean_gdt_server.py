# -*- coding: utf-8 -*-
"""Module to solve overlaps in GFF3 files."""

import concurrent.futures as cf
import sys
import traceback
from pathlib import Path
from typing import TYPE_CHECKING, cast

import multiprocessing as mp
from queue import Empty

import pandas as pd

from . import clean, gff3_utils, log_setup

if TYPE_CHECKING:
    import argparse

    import gdt  # type: ignore[import-not-found]


def format_attributes(
    df: "pd.DataFrame",
    an_index: int,
    request_queue: "mp.Queue[tuple[int, list[str]]]",
    response_queue: "mp.Queue[list[str]]",
    log: log_setup.TempLogger,
) -> "pd.DataFrame":
    """Format attributes in the DataFrame using GDT."""
    gene_ids_to_query = df.loc[1:, "gene_id"].dropna().tolist()

    try:
        # Send batch request
        request_queue.put((an_index, gene_ids_to_query))
        labels_list = response_queue.get()

        # Check if we got a valid response with matching length
        if len(labels_list) != len(gene_ids_to_query):
            log.error(f"Invalid response from dictionary server. Using default format.")
            log.trace(
                f"Expected {len(gene_ids_to_query)} labels, got {len(labels_list)}."
            )
            log.trace(f"Gene IDs: {gene_ids_to_query}, Labels: {labels_list}")
            labels_list = gene_ids_to_query  # fallback to gene_ids

    except Exception as e:
        log.error(f"Error querying gene IDs: {e}. Using default format.")
        labels_list = gene_ids_to_query  # fallback to gene_ids

    # Build attributes string using labels_list directly
    df.loc[1:, "attributes"] = (
        "name="
        + pd.Series(labels_list, index=df.loc[1:].index)
        + ";source="
        + df.loc[1:, "seqid"]
        + "|"
        + df.loc[1:, "type"]
        + "|"
        + df.loc[1:, "start"].astype(str)
        + "|"
        + df.loc[1:, "end"].astype(str)
        + "|"
        + df.loc[1:, "strand"]
        + "|"
        + df.loc[1:, "gene_id"]
        + ";"
    )

    return df


def clean_an_gdt_server(
    an_index: int,
    log: log_setup.TempLogger,
    gff_in: Path,
    gff_out: Path,
    request_queue: "mp.Queue[tuple[int, list[str]]]",
    response_queue: "mp.Queue[list[str]]",
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

        df = format_attributes(df, an_index, request_queue, response_queue, log)

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

    gff3_utils.check_files(log, tsv, gff_in_builder, an_column, should_exist=True)

    if not overwrite:
        gff3_utils.check_files(log, tsv, gff_out_builder, an_column, should_exist=False)

    # create server for dictionary
    try:
        log.info("Starting dictionary server...")
        server_process, request_queue, response_queues = create_dict_server(
            gdict=gdict,
            num_queues=len(tsv),
        )
        log.info(
            f"Creating tasks for processing {tsv.shape[0]} ANs with {workers} workers..."
        )
        with cf.ProcessPoolExecutor(max_workers=workers) as executor:
            tasks = []
            for i, an in enumerate(tsv[an_column]):
                task = executor.submit(
                    clean_an_gdt_server,
                    i,
                    log.spawn_buffer(),
                    gff_in_builder.build(an),
                    gff_out_builder.build(an),
                    request_queue,
                    response_queues[i],
                    query_string,
                    keep_orfs,
                )
                tasks.append(task)

            log.info("Tasks submitted, waiting for completion...")

            for future in cf.as_completed(tasks):
                success, an, records = future.result()
                if not success:
                    log.error(f"[{an}] Failed to process.")

                for record in records:
                    log.log(record[0], record[1])

            log.info("All processed.")
    finally:
        # Stop the dictionary server gracefully
        stop_dict_server(server_process, request_queue)
        log.info("Dictionary server stopped.")


def solve_gdt_call_server(
    log: log_setup.GDTLogger, args: "argparse.Namespace", workers: int
) -> tuple[bool, str, list[log_setup._RawMsg]] | None:
    """Handle `clean` command with an GDT gdict.

    Args:
        log: Logger instance
        args: Parsed command line arguments
        workers: Number of worker processes to use

    """
    try:
        import gdt
    except ImportError:
        raise SystemExit(
            "GDT package not found. Please install it to use the --gdict option."
        )

    log.debug(f"GDT package version: {gdt.__version__}")
    gdt_path: Path = Path(args.gdt).resolve()
    if not gdt_path.is_file():
        log.error(f"GDT .gdict file not found: {gdt_path}")
        sys.exit(1)

    gdict = gdt.read_gdict(gdt_path, lazy_info=False)
    log.info(f"GDT dictionary loaded from {gdt_path}")
    gdt.log_info(log, gdict)

    clean_multiple_gdt(
        log,
        args.tsv,
        workers,
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

    return None


def run_dictionary_server(
    gdict: "gdt.GeneDict",
    request_queue: "mp.Queue[tuple[int, list[str]]]",  # Changed: now expects list[str]
    response_queues: "dict[int, mp.Queue[list[str]]]",  # Changed: now sends list[str] in order
) -> None:
    """
    Dictionary server function that runs in a separate process.

    Args:
        gdict: The dictionary to serve
        request_queue: Queue to receive (worker_id, gene_id_list) requests
        response_queues: Dict of worker_id -> Queue for sending ordered responses
    """
    query_count = 0
    an_count = 0

    while True:
        try:
            # Get request: (worker_id, list_of_gene_ids)
            an_index, gene_id_list = request_queue.get(timeout=1.0)

            # Check for shutdown signal
            if an_index == -1 and gene_id_list == ["shutdown"]:
                print(
                    f"Dictionary server shutting down after {query_count} batch queries, from {an_count} ANs."
                )
                break

            # NEW: Process gene_ids in order and return ordered results
            result_list = []
            for gene_id in gene_id_list:
                result = gdict.get(gene_id, "not_found")
                if result != "not_found":
                    result_list.append(result.label)
                else:
                    result_list.append("not_found")
                query_count += 1
            an_count += 1

            # Send ordered response back to specific AN's queue
            if an_index in response_queues:
                response_queues[an_index].put(result_list)
            else:
                print(f"Warning: Unknown an_index {an_index}")

        except Empty:
            # Timeout - continue loop to check for shutdown
            continue
        except Exception as e:
            print(f"Dictionary server error: {e}")
            continue


def create_dict_server(
    gdict: "gdt.GeneDict",
    num_queues: int,
) -> tuple[
    mp.Process,
    "mp.Queue[tuple[int, list[str]]]",
    dict[int, "mp.Queue[list[str]]"],
]:
    """
    Create and start a dictionary server process.

    Args:
        gdict: GDT GeneDict instance to serve
        num_queues: Number of response queues to create

    Returns:
        Tuple containing:
        - server_process: The dictionary server process
        - request_queue: Queue for sending (queue_id, gene_id) requests
        - response_queues: Dict mapping queue_id to response queues
    """
    # Use Manager for queues that need to be shared across ProcessPoolExecutor
    manager = mp.Manager()
    request_queue = cast("mp.Queue[tuple[int, list[str]]]", manager.Queue())
    response_queues = cast(
        dict[int, "mp.Queue[list[str]]"],
        {i: manager.Queue() for i in range(num_queues)},
    )

    # Start server process
    server_process = mp.Process(
        target=run_dictionary_server, args=(gdict, request_queue, response_queues)
    )
    server_process.start()

    return server_process, request_queue, response_queues


def stop_dict_server(
    server_process: mp.Process, request_queue: "mp.Queue[tuple[int, list[str]]]"
) -> None:
    """
    Stop the dictionary server gracefully.

    Args:
        server_process: The server process to stop
        request_queue: Request queue to send shutdown signal
    """
    if server_process and server_process.is_alive():
        # Send shutdown signal
        request_queue.put((-1, ["shutdown"]))
        server_process.join(timeout=5)

        if server_process.is_alive():
            print("Server didn't shutdown gracefully, terminating...")
            server_process.terminate()
            server_process.join()

        print("Dictionary server stopped")
