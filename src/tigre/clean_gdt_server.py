# -*- coding: utf-8 -*-
"""Module to solve overlaps in GFF3 files."""

import concurrent.futures as cf
import multiprocessing as mp
import sys
import traceback
from pathlib import Path
from queue import Empty
from typing import TYPE_CHECKING, TypeAlias, cast

import pandas as pd

from . import clean, gff3_utils, log_setup

if TYPE_CHECKING:
    import argparse

    import gdt  # type: ignore[import-not-found]

_ReqQueue: TypeAlias = "mp.Queue[tuple[list[str], mp.connection.Connection]]"


def format_attributes(
    df: "pd.DataFrame",
    queue: _ReqQueue,
    log: log_setup.TempLogger,
) -> "pd.DataFrame":
    """Format attributes in the DataFrame using GDT."""
    gene_ids_to_query = df.loc[1:, "gene_id"].dropna().tolist()

    try:
        # Create pipe for response and send request
        worker_conn, server_conn = mp.Pipe()
        queue.put((gene_ids_to_query, server_conn))
        labels_list = worker_conn.recv()  # Get response through pipe
        worker_conn.close()  # Close our end of the pipe

        # Check if we got a valid response with matching length
        if len(labels_list) != len(gene_ids_to_query):
            log.error("Invalid response from dictionary server. Using default format.")
            log.trace(
                f"Expected {len(gene_ids_to_query)} labels, got {len(labels_list)}."
            )
            log.trace(f"Gene IDs: {gene_ids_to_query}, Labels: {labels_list}")
            labels_list = gene_ids_to_query  # fallback to gene_ids

    except Exception as e:
        log.error(f"Error querying gene IDs: {e}. Using default format.")
        labels_list = gene_ids_to_query  # fallback to gene_ids

    # Log missing gene_ids and replace "not_found" with original gene_id
    for i, (gene_id, label) in enumerate(zip(gene_ids_to_query, labels_list)):
        if label == "not_found":
            log.error(
                f"Gene ID {gene_id} not found in GDT dictionary. Using default format."
            )
            labels_list[i] = gene_id  # replace "not_found" with original gene_id

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
    log: log_setup.TempLogger,
    gff_in: Path,
    gff_out: Path,
    req_queue: _ReqQueue,
    query_string: str,
    keep_orfs: bool,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Clean a GFF3 file by removing unnecessary features and attributes with a gdict.

    Args:
        an_index: Index of the AN being processed
        log: Logger instance
        gff_in: Path to the input GFF3 file
        gff_out: Path to the output GFF3 file
        req_queue: Request queue for the dictionary server
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

        df = format_attributes(df, req_queue, log)

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
        server_process, req_queue = create_dict_server(gdict)
        log.info(f"Submitiing tasks for {tsv.shape[0]} ANs with {workers} workers...")
        with cf.ProcessPoolExecutor(max_workers=workers) as executor:
            tasks = []
            for an in tsv[an_column]:
                task = executor.submit(
                    clean_an_gdt_server,
                    log.spawn_buffer(),
                    gff_in_builder.build(an),
                    gff_out_builder.build(an),
                    req_queue,
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
        stop_dict_server(server_process, req_queue)
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


def run_dictionary_server(gdict: "gdt.GeneDict", req_queue: _ReqQueue) -> None:
    """Run the dictionary server using a shared queue for communication."""
    query_count = 0

    while True:
        try:
            # Get request: (gene_ids_list, response_pipe)
            gene_ids, worker_conn = req_queue.get(timeout=1.0)

            if gene_ids == ["shutdown"]:
                worker_conn.close()
                break

            # Process batch lookup
            results = []
            for gene_id in gene_ids:
                result = gdict.get(gene_id, "not_found")
                results.append(result.label if result != "not_found" else "not_found")
                query_count += 1

            # Send response back through pipe and close it
            try:
                worker_conn.send(results)
                worker_conn.close()
            except Exception as e:
                print(f"Error sending response: {e}")

        except Empty:
            # Timeout - continue loop to check for shutdown
            continue
        except Exception as e:
            print(f"Dictionary server error: {e}")
            continue

    print(f"Dictionary server shutting down after {query_count} queries")


def create_dict_server(gdict: "gdt.GeneDict") -> tuple[mp.Process, _ReqQueue]:
    """Create dictionary server with a shared queue."""
    # Single shared queue for all requests
    manager = mp.Manager()
    request_queue = cast(_ReqQueue, manager.Queue())

    server_process = mp.Process(
        target=run_dictionary_server, args=(gdict, request_queue)
    )
    server_process.start()

    return server_process, request_queue


def stop_dict_server(server_process: mp.Process, req_queue: _ReqQueue) -> None:
    """Stop the dictionary server gracefully.

    Args:
        server_process: The server process to stop
        req_queue: Request queue to send shutdown signal

    """
    if server_process and server_process.is_alive():
        try:
            # Send shutdown signal
            req_queue.put((["shutdown"], mp.Pipe()[1]))
        except Exception as e:
            print(f"Error during shutdown signal: {e}")

        server_process.join(timeout=5)

        if server_process.is_alive():
            print("Server didn't shutdown gracefully, terminating...")
            server_process.terminate()
            server_process.join()

        print("Dictionary server stopped")
