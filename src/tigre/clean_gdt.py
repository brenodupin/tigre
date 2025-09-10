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
    import gdt  # type: ignore[import-not-found]


_ReqQueue: TypeAlias = "mp.queues.Queue[tuple[list[str], mp.connection.Connection]]"


def get_labels(
    subset: pd.DataFrame,
    queue: _ReqQueue,
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


def format_gdt(
    df: pd.DataFrame,
    queue: _ReqQueue,
    log: log_setup.TempLogger,
) -> pd.DataFrame:
    """Gather and format the DataFrame to ensure correct data types."""
    attrs = df["attributes"]

    name_up = attrs.str.extract(clean._RE_name_up, expand=False)  # type: ignore[call-overload]
    name_left = attrs.str.extract(clean._RE_name_left, expand=False)  # type: ignore[call-overload]
    name_right = attrs.str.extract(clean._RE_name_right, expand=False)  # type: ignore[call-overload]
    name_dw = attrs.str.extract(clean._RE_name_dw, expand=False)  # type: ignore[call-overload]
    name_general = attrs.str.extract(clean._RE_name, expand=False)  # type: ignore[call-overload]

    source_up = attrs.str.extract(clean._RE_source_up, expand=False)  # type: ignore[call-overload]
    source_left = attrs.str.extract(clean._RE_source_left, expand=False)  # type: ignore[call-overload]
    source_right = attrs.str.extract(clean._RE_source_right, expand=False)  # type: ignore[call-overload]
    source_dw = attrs.str.extract(clean._RE_source_dw, expand=False)  # type: ignore[call-overload]
    source_general = attrs.str.extract(clean._RE_source, expand=False)  # type: ignore[call-overload]

    df["name_left"] = name_up.where(
        name_up.notna(), name_left.where(name_left.notna(), name_general)
    )

    df["name_right"] = name_dw.where(
        name_dw.notna(), name_right.where(name_right.notna(), name_general)
    )

    df["source_left"] = source_up.where(
        source_up.notna(), source_left.where(source_left.notna(), source_general)
    )

    df["source_right"] = source_dw.where(
        source_dw.notna(), source_right.where(source_right.notna(), source_general)
    )

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

        gene_labels = get_labels(subset, queue, log)

        df.loc[og_row, "name_left"] = df.loc[og_row, "name_left"].combine_first(
            gene_labels
        )

        df.loc[og_row, "name_right"] = df.loc[og_row, "name_right"].combine_first(
            gene_labels
        )

        df.loc[og_row, "source_left"] = df.loc[og_row, "source_left"].combine_first(
            source_fallback
        )

        df.loc[og_row, "source_right"] = df.loc[og_row, "source_right"].combine_first(
            source_fallback
        )

        df.loc[og_row, "attributes"] = [
            f"name={name};source={source};"
            for name, source in zip(
                df.loc[og_row, "name_left"],
                df.loc[og_row, "source_left"],
            )
        ]

    return df


def clean_an_gdt(
    log: log_setup.TempLogger,
    gff_in: Path,
    gff_out: Path,
    req_queue: _ReqQueue,
    query_string: str,
    keep_orfs: bool,
    ext_filter: bool = False,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Clean a GFF3 file by removing unnecessary features and attributes.

    Args:
        log: Logger instance
        gff_in: Path to the input GFF3 file
        gff_out: Path to the output GFF3 file
        req_queue: Queue to communicate with the dictionary server
        clean_func: Function to clean attributes in each row
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
        header = clean.create_header(an, df.at[0, "attributes"], df.at[0, "end"])
        df["gene_id"] = df["attributes"].str.extract(gff3_utils._RE_ID, expand=False)  # type: ignore[call-overload]

        if (df["type"] == "region").sum() > 1:
            log.warning("Multiple regions found:")
            for row in df[df["type"] == "region"].itertuples():
                log.warning(f" s: {row.start} e: {row.end} | {row.attributes}")
            log.warning(" Always choosing the first region and removing the rest.")

        df_region = df.iloc[0:1].copy()
        df = df[df["type"] != "region"]
        df = format_gdt(df, req_queue, log)

        if df["end"].max() > df_region.at[0, "end"]:
            df = clean.bigger_than_region_solver(log, df, df_region.at[0, "end"])

        df = df.sort_values(by=["start", "type"], ascending=True, ignore_index=True)

        log.trace(" Overlaps solver...")
        clean._get_column_indices(df)
        df = clean.overlap_solver(log, df, df_region)

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
        log.info(f"Submitiing tasks for {tsv.shape[0]} ANs with {workers} workers...")
        with cf.ProcessPoolExecutor(max_workers=workers) as executor:
            tasks = []
            for an in tsv[an_column]:
                task = executor.submit(
                    clean_an_gdt,
                    log.spawn_buffer(),
                    gff_in_builder.build(an),
                    gff_out_builder.build(an),
                    req_queue,
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


def clean_single_gdt(
    log: log_setup.GDTLogger,
    gff_in: Path,
    gff_out: Path,
    gdict: "gdt.GeneDict",
    query_string: str = gff3_utils.QS_GENE_TRNA_RRNA_REGION,
    keep_orfs: bool = False,
    ext_filter: bool = False,
) -> None:
    """Clean a single GFF3 file using GDT for gene ID lookups.

    Args:
        log: Logger instance
        gff_in: Path to the input GFF3 file
        gff_out: Path to the output GFF3 file
        gdict: GDT GeneDict instance
        query_string: Query string to filter the GFF3 file
        keep_orfs: Whether to keep ORFs in the GFF3 file
        ext_filter: Whether to use extended filtering for ORFs

    """
    try:
        log.info("Starting dictionary server...")
        server_process, req_queue = create_dict_server(gdict, log.spawn_buffer())

        log.info(f"Processing file: {gff_in}")
        success, an, records = clean_an_gdt(
            log.spawn_buffer(),
            gff_in,
            gff_out,
            req_queue,
            query_string,
            keep_orfs,
            ext_filter,
        )
        for record in records:
            log.log(record[0], record[1])

        if not success:
            log.error(f"[{an}] Failed to process.")
            sys.exit(1)

        log.info("Processing completed successfully.")

    except Exception:
        error_msg = traceback.format_exc()
        log.error(f"Error in {gff_in}:\n{error_msg}")

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
