# -*- coding: utf-8 -*-
"""Module to solve overlaps in GFF3 files."""

import concurrent.futures as cf
import multiprocessing as mp
import sys
from functools import partial
from pathlib import Path
from queue import Empty
from typing import TYPE_CHECKING, TypeAlias, cast

import pandas as pd

from . import clean, gff3_utils, log_setup

if TYPE_CHECKING:
    import gdt  # type: ignore[import-not-found]


_ReqQueue: TypeAlias = "mp.queues.Queue[tuple[list[str], mp.connection.Connection]]"


def get_names_server(
    queue: _ReqQueue,
    subset: pd.DataFrame,
    log: log_setup.TempLogger,
) -> "pd.Series[str]":
    """Get labels for a subset of gene IDs using the dictionary server."""
    if subset.empty:
        return pd.Series([], dtype=str)

    gene_ids: list[str] = subset["gene_id"].tolist()
    try:
        worker_conn, server_conn = mp.Pipe()
        queue.put((gene_ids, server_conn))
        results: list[str] = worker_conn.recv()
        worker_conn.close()

        # Create a Series with the same index as the subset
        return pd.Series(results, index=subset.index)

    except Exception as e:
        log.error(f"Error communicating with dictionary server: {e}")
        return pd.Series(gene_ids, index=subset.index)


def get_names_gdt(
    gdict: "gdt.GeneDict",
    subset: pd.DataFrame,
    log: log_setup.TempLogger,
) -> "pd.Series[str]":
    """Get labels for a subset of gene IDs using the GeneDict."""
    if subset.empty:
        return pd.Series([], dtype=str)

    results: list[str] = []
    for gene_id in subset["gene_id"]:
        result = gdict.get(gene_id, None)
        if result is not None:
            results.append(result.label)
        else:
            log.warning(f"Gene ID '{gene_id}' not found in dictionary")
            results.append(gene_id)

    return pd.Series(results, index=subset.index)


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
    ext_filter: bool = False,
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
        ext_filter: Whether to use extended filtering for ORFs

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
        names_func = partial(get_names_server, req_queue)

        log.info(f"Submitiing tasks for {tsv.shape[0]} ANs with {workers} workers...")
        with cf.ProcessPoolExecutor(max_workers=workers) as executor:
            tasks = []
            for an in tsv[an_column]:
                task = executor.submit(
                    clean.clean_an,
                    log.spawn_buffer(),
                    gff_in_builder.build(an),
                    gff_out_builder.build(an),
                    names_func,
                    query_string,
                    keep_orfs,
                    ext_filter,
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
