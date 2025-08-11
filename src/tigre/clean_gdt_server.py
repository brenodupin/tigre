# -*- coding: utf-8 -*-
"""Module to solve overlaps in GFF3 files."""

import concurrent.futures as cf
import multiprocessing as mp
import traceback
from pathlib import Path
from queue import Empty
from typing import TYPE_CHECKING, TypeAlias, cast

import pandas as pd

from . import clean, gff3_utils, log_setup

if TYPE_CHECKING:

    import gdt  # type: ignore[import-not-found]

_ReqQueue: TypeAlias = "mp.queues.Queue[tuple[list[str], mp.connection.Connection]]"


def format_attributes(
    df: "pd.DataFrame",
    queue: _ReqQueue,
    log: log_setup.TempLogger,
) -> "pd.DataFrame":
    """Format attributes in the DataFrame using GDT."""
    gene_query = df.loc[1:, "gene_id"].tolist()
    labels_list: list[str] = gene_query  # default format

    worker_conn: mp.connection.Connection | None = None
    try:
        worker_conn, server_conn = mp.Pipe()
        queue.put(
            (gene_query, server_conn)
        )  # Create pipe for response and send request
        received_labels = worker_conn.recv()
        worker_conn.close()

        if isinstance(received_labels, list) and len(received_labels) == len(
            gene_query
        ):
            log.trace(f"Got {len(received_labels)} labels from dictionary server.")
            labels_list = received_labels

        else:
            log.error("Invalid response from dictionary server. Using default format.")
            log.trace(f"Expected {len(gene_query)} labels, got {len(received_labels)}.")
            log.trace(f"Gene IDs: {gene_query}, Labels: {labels_list}")

    except Exception as e:
        log.error(f"Error querying gene IDs: {e}. Using default format.")
        if worker_conn:
            try:
                worker_conn.close()
            except Exception:
                pass  # Already closed or connection error

    for i, (gene_id, label) in enumerate(zip(gene_query, labels_list)):
        if label == "not_found":
            log.error(
                f"Gene ID {gene_id} not found in GDT dictionary. Using default format."
            )
            labels_list[i] = gene_id  # replace "not_found" with original gene_id

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
        error_msg = traceback.format_exc()
        if "an" not in locals():
            log.error(f"Error in {gff_in}:\n{error_msg}")
            return False, "Not defined", log.get_records()

        log.error(f"[{an}] Error:\n{error_msg}")
        log.debug(f"[{an}] -- End --")
        return False, an, log.get_records()


def clean_multiple_gdt_server(
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
        workers: Number of worker processes to use (including the Dictionary Server)
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
    # subtracting the server process from the count, but it should never be 0
    workers = max(1, workers - 1)

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

    try:
        log.info("Starting dictionary server...")
        server_process, req_queue = create_dict_server(gdict, log.spawn_buffer())
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
        stop_dict_server(server_process, req_queue, log)


def run_dictionary_server(
    gdict: "gdt.GeneDict",
    req_queue: _ReqQueue,
    log: log_setup.TempLogger,
) -> None:
    """Run the dictionary server using a shared queue for communication."""
    query_count = 0
    batch_count = 0

    while True:
        try:
            # Get request: (gene_ids_list, response_pipe)
            gene_ids, worker_conn = req_queue.get(timeout=1.0)

            if gene_ids == ["shutdown"]:
                # Send log records back through the shutdown pipe
                try:
                    log.info(
                        f"Dictionary server resolved {query_count} queries in "
                        f"{batch_count} batches."
                    )
                    worker_conn.send(log.get_records())

                except Exception as e:
                    print(f"Error sending log records: {e}")

                finally:
                    worker_conn.close()
                    break

            # Process batch lookup
            results: list[str] = []
            for gene_id in gene_ids:
                result = gdict.get(gene_id, None)
                results.append(result.label if result else "not_found")
                query_count += 1
            batch_count += 1

            # log.trace(f"Queue got: {gene_ids}")
            # log.trace(f"Pipe sent: {results}")

            try:
                worker_conn.send(results)
                worker_conn.close()
            except Exception as e:
                log.error(f"Error sending response to {worker_conn}: {e}")

        except Empty:
            # Timeout - continue loop to check for shutdown
            continue
        except Exception as e:
            log.error(f"Dictionary server error: {e}")
            continue


def create_dict_server(
    gdict: "gdt.GeneDict",
    server_log: log_setup.TempLogger,
) -> tuple[mp.Process, _ReqQueue]:
    """Create dictionary server with a shared queue."""
    manager = mp.Manager()
    request_queue = cast(_ReqQueue, manager.Queue())  # shared queue for all requests

    server_process = mp.Process(
        target=run_dictionary_server,
        args=(gdict, request_queue, server_log),
    )
    server_process.start()

    return server_process, request_queue


def stop_dict_server(
    server_process: mp.Process,
    req_queue: _ReqQueue,
    log: "log_setup.GDTLogger",
) -> None:
    """Stop the dictionary server gracefully and collect its logs."""
    if server_process and server_process.is_alive():
        try:

            worker_conn, server_conn = mp.Pipe()
            req_queue.put((["shutdown"], server_conn))

            # Receive log records through the shutdown pipe
            records = worker_conn.recv()
            worker_conn.close()

            for record in records:
                log.log(record[0], record[1])

        except Exception as e:
            log.error(f"Error during shutdown or log collection: {e}")

        server_process.join(timeout=5)

        if server_process.is_alive():
            log.warning("Server didn't shutdown gracefully, terminating...")
            server_process.terminate()
            server_process.join()

    log.info("Dictionary server stopped")
