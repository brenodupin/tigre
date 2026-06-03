# -*- coding: utf-8 -*-
"""Module to solve overlaps in GFF3 files."""

import concurrent.futures as cf
import traceback
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Literal

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


def pre_filter_gff(
    log: log_setup.TempLogger,
    gff_in: Path,
    gff_out: Path,
    keep_type: str,
    clean_names_func: clean._NamesType,
    query_expr: "pl.Expr",
    keep_orfs: bool,
    ext_filter: bool = False,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Pre-filter GFF3 file to keep only relevant features for overlap cleaning."""
    try:
        log.debug(f"[{gff_in.name}] -- Start --")
        df = gff3_utils.load_gff3(
            gff_in,
            usecols=gff3_utils.GFF3_COLUMNS,
            query_expr=query_expr,
            return_polars=True,
        )
        if not keep_orfs:
            df = gff3_utils.filter_orfs(df, extended=ext_filter)

        # We expect the first row to be the region feature.
        an, region_att, region_end = (
            df.item(0, "seqid"),
            df.item(0, "attributes"),
            df.item(0, "end"),
        )
        df = df.filter(pl.col("type") != "region")

        header = clean.create_header(an, region_att, region_end, file_name="genes.py")

        df = df.with_columns(
            pl.col("attributes")
            .str.extract(gff3_utils._RS_ID, group_index=1)
            .alias("gene_id")
        )
        gene_ids = df["gene_id"].unique().to_list()

        log.trace(f"Extracted {len(gene_ids)} unique genes. {gene_ids}")
        results, log = clean_names_func(gene_ids, log)
        log.trace(f"Received gene labels for {len(results)} genes. {results}")

        results_dict = dict(zip(gene_ids, results))
        df = df.with_columns(
            pl.col("gene_id").replace_strict(results_dict).alias("gene_label_raw")
        )
        df = df.with_columns(
            pl.col("gene_label_raw").str.replace_all(r"^MIT-|^PLT-", "").alias("gene_label")
        )
        labels = set(df["gene_label"].unique().to_list())
        log.trace(f"Processed gene labels: {labels}")

        keep_type = "gene trna rrna"

        labels_to_keep = set()

        if "trna" in keep_type:
            labels_to_keep |= trna_labels

        if "rrna" in keep_type:
            labels_to_keep |= rrna_labels

        if "gene" in keep_type:
            # for gene keep, remove rows that are in the trna_labels or rrna_labels sets
            gene_labels = set(labels) - trna_labels - rrna_labels
            labels_to_keep |= gene_labels

        df = df.filter(pl.col("gene_label").is_in(labels_to_keep))

        df = df.with_columns(
            pl.concat_str(
                ["seqid", "type", "start", "end", "strand", "gene_id"],
                separator="|",
            ).alias("_source_fallback"),
        )
        df = df.with_columns(
            pl.Series(
                [
                    f"ID={seqid}_{ftype}_{i};name_up={name};source_up={source};name_dw={name};source_dw={source};"
                    for i, (seqid, ftype, name, source) in enumerate(
                        zip(
                            df["seqid"],
                            df["type"],
                            df["gene_label_raw"],
                            df["_source_fallback"],
                        ),
                        start=1,
                    )
                ]
            ).alias("attributes")
        )
        df = df.drop(["gene_id", "gene_label", "gene_label_raw", "_source_fallback"])

        # add header to the file
        with open(gff_out, "w", encoding="utf-8") as file_handler:
            file_handler.write(header)
            df.write_csv(file_handler, separator="\t", include_header=False)

        return True, an, log.get_records()

    except Exception:
        an_error: str = an if "an" in locals() else gff_in.name
        error_msg = traceback.format_exc()
        log.error(f"Error in {an_error}:\n{error_msg}")
        return False, an_error, log.get_records()


def genes_multiple(
    log: log_setup.GDTLogger,
    tsv_path: Path,
    workers: int,
    gdict: "gdt.GeneDict",
    gff_in_ext: str = ".gff3",
    gff_in_suffix: str = "",
    gff_out_ext: str = ".gff3",
    gff_out_suffix: str = "_genes",
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
        with pl.Config(fmt_table_cell_list_len=-1):
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

        log.info(f"Submiting tasks for {tsv.shape[0]} ANs with {workers} workers...")
        with cf.ProcessPoolExecutor(max_workers=workers) as executor:
            tasks = []
            for an in tsv[an_column]:
                task = executor.submit(
                    pre_filter_gff,
                    log.spawn_buffer(),
                    gff_in_builder.build(an),
                    gff_out_builder.build(an),
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
