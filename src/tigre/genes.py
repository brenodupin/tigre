# -*- coding: utf-8 -*-
"""Module to solve overlaps in GFF3 files."""

import concurrent.futures as cf
import multiprocessing as mp
import tempfile
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Literal, TypeAlias

import polars as pl

from . import clean, clean_gdt, gff3_utils, log_setup

if TYPE_CHECKING:
    import gdt

trna_labels = {
    "TA",
    "TASX",
    "TC",
    "TD",
    "TE",
    "TER",
    "TF",
    "TFM",
    "TG",
    "TGLX",
    "TH",
    "TI",
    "TK",
    "TL1|2",
    "TM",
    "TN",
    "TP",
    "TQ",
    "TR",
    "TRMR",
    "TRNX",
    "TS1|2",
    "TSEC",
    "TSUP",
    "TT",
    "TU",
    "TV",
    "TW",
    "TX",
    "TXLE",
    "TY",
}

rrna_labels = {
    "RNR12",
    "RNR15",
    "RNR16",
    "RNR18",
    "RNR21",
    "RNR23",
    "RNR26",
    "RNR28",
    "RNR4",
    "RNR4.5",
    "RNR5",
    "RNR9",
    "RNRL",
    "RNRS",
    "RRNA",
}

_ReqQueue: TypeAlias = "mp.queues.Queue[tuple[list[str], mp.connection.Connection]]"
_NamesTypeNoPandas: TypeAlias = (
    "Callable[[list[str], log_setup.TempLogger], tuple[list[str], log_setup.TempLogger]]"
)


def get_names_server(
    queue: _ReqQueue,
    gene_ids: list[str],
    log: log_setup.TempLogger,
) -> "tuple[list[str], log_setup.TempLogger]":
    """Get labels for a subset of gene IDs using the dictionary server."""
    try:
        worker_conn, server_conn = mp.Pipe()
        queue.put((gene_ids, server_conn))
        results: list[str] = worker_conn.recv()
        worker_conn.close()

        # Create a Series with the same index as the subset
        return results, log

    except Exception as e:
        log.error(f"Error communicating with dictionary server: {e}")
        return [], log


def get_names_gdt(
    gdict: "gdt.GeneDict",
    gene_list: list[str],
    log: log_setup.TempLogger,
) -> "tuple[list[str], log_setup.TempLogger]":
    """Get labels for a subset of gene IDs using the GeneDict."""
    results: list[str] = []
    for gene_id in gene_list:
        result = gdict.get(gene_id, None)
        if result is not None:
            results.append(result.label)
        else:
            log.warning(f"Gene ID '{gene_id}' not found in dictionary")
            results.append(gene_id)

    return results, log


def pre_filter_gff(
    log: log_setup.TempLogger,
    gff_in: Path,
    gff_out: Path,
    pre_names_func: _NamesTypeNoPandas,
    keep_type: Literal["gene", "trna", "rrna"],
    clean_names_func: clean._NamesType,
    query_expr: "pl.Expr",
    keep_orfs: bool,
    ext_filter: bool = False,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Pre-filter GFF3 file to keep only relevant features for overlap cleaning."""
    log.debug(f"[{gff_in.name}] -- Start --")
    df = gff3_utils.load_gff3(
        gff_in,
        usecols=gff3_utils.GFF3_COLUMNS,
        query_expr=query_expr,
        return_polars=True,
    )
    # save first row of df (polars) as region line, then remove it from the df
    an = df.item(0, "seqid")
    header = clean.create_header(an, df.item(0, "attributes"), df.item(0, "end"))

    region_df = df.head(1)
    df = df.slice(1)

    df = df.with_columns(
        pl.col("attributes")
        .str.extract(gff3_utils._RE_ID.pattern, group_index=1)
        .alias("gene_id")
    )
    gene_ids = df["gene_id"].unique().to_list()

    log.trace(f"Extracted {len(gene_ids)} unique genes. {gene_ids}")
    results, log = pre_names_func(gene_ids, log)
    log.trace(f"Received gene labels for {len(results)} genes. {results}")
    results = [r.replace("MIT-", "").replace("PLT-", "") for r in results]
    log.trace(f"Cleaned gene labels: {results}")

    results_dict = dict(zip(gene_ids, results))
    df = df.with_columns(pl.col("gene_id").replace_strict(results_dict).alias("gene_label"))

    if keep_type == "trna":
        df = df.filter(pl.col("gene_label").is_in(trna_labels))

    elif keep_type == "rrna":
        df = df.filter(pl.col("gene_label").is_in(rrna_labels))

    elif keep_type == "gene":
        # for gene keep, remove rows that are in the trna_labels or rrna_labels sets
        set_to_remove = trna_labels.union(rrna_labels)
        log.trace(f"Filtering out features with labels in {set_to_remove}")
        df = df.filter(~pl.col("gene_label").is_in(set_to_remove))

    # write df to a tmpfile and call clean_an with the tmpfile as input
    df = df.drop(["gene_label", "gene_id"])
    df = pl.concat([region_df, df])

    with tempfile.NamedTemporaryFile(mode="w", suffix=".gff3") as tmp_gff:
        tmp_path = Path(tmp_gff.name).resolve()
        log.trace(f"Writing pre-filtered GFF3 to temporary file: {tmp_path}")
        with open(tmp_gff.name, "w") as f:
            f.write("".join(header) + "\n")
            df.write_csv(f, separator="\t", include_header=False)

        return clean.clean_an(
            log,
            tmp_path,
            gff_out,
            clean_names_func,
            query_expr,
            keep_orfs,
            ext_filter,
        )


def genes_multiple(
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
    keep_type: Literal["gene", "trna", "rrna"] = "gene",
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
        keep_type: Type of features to keep ("gene", "trna", or "rrna")

    """
    gff3_utils._ensure_spawn(log)  # Ensure spawn method is set for multiprocessing

    try:
        query_expr = gff3_utils._parse_query_to_polars(query_string)
        log.info(f"Using query expression: {query_expr}")
    except Exception as e:
        log.error(f"Error parsing query string '{query_string}': {e}")
        raise

    tsv = pl.read_csv(tsv_path, separator="\t")
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
        server_process, req_queue = clean_gdt.create_dict_server(gdict, log.spawn_buffer())
        clean_names_func = partial(clean_gdt.get_names_server, req_queue)
        pre_names_func = partial(get_names_server, req_queue)

        log.info(f"Submitiing tasks for {tsv.shape[0]} ANs with {workers} workers...")
        with cf.ProcessPoolExecutor(max_workers=workers) as executor:
            tasks = []
            for an in tsv[an_column]:
                task = executor.submit(
                    pre_filter_gff,
                    log.spawn_buffer(),
                    gff_in_builder.build(an),
                    gff_out_builder.build(an),
                    pre_names_func,
                    keep_type,
                    clean_names_func,
                    query_expr,
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
        clean_gdt.stop_dict_server(server_process, req_queue, log)
