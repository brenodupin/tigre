# -*- coding: utf-8 -*-
"""Module to extract sequences from GFF3 files using Biopython."""

import concurrent.futures as cf
import traceback
from pathlib import Path

import polars as pl

from . import gff3_utils, log_setup

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    SeqIO = None  # type: ignore
    SeqRecord = None  # type: ignore


def getfasta_single(
    log: log_setup.TempLogger,
    gff_in: Path,
    fasta_in: Path,
    fasta_out: Path,
    bedtools_compatible: bool = False,
    skip_region: bool = False,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Extract sequences from a FASTA file using Biopython."""
    try:
        log.trace(f"biopython_getfasta: {gff_in = } | {fasta_in = } | {fasta_out = }")
        log.debug(
            f"biopython_getfasta: gff_in: {gff_in.name} | fasta_in: {fasta_in.name} | "
            f"fasta_out: {fasta_out.name}"
        )

        idx_fix = 0 if bedtools_compatible else 1
        fasta = SeqIO.read(fasta_in, "fasta").seq  # type: ignore[no-untyped-call]
        records = []

        df = gff3_utils.load_gff3(
            gff_in,  # idx 0      1       2       3
            usecols=("seqid", "type", "start", "end", "attributes"),
            return_polars=True,
        )
        df = df.with_columns(
            [
                pl.col("start").cast(pl.Int64),
                pl.col("end").cast(pl.Int64),
            ]
        )

        if skip_region:
            df = df.filter(pl.col("type") != "region")

        if df.is_empty():
            log.warning(f"Empty GFF3 file: {gff_in.name}, creating empty FASTA.")
            SeqIO.write([], fasta_out, "fasta")
            return True, gff_in.name, log.get_records()

        seqid: str = df.item(0, "seqid")
        df = df.with_columns((pl.col("start") - 1).alias("start"))

        has_ig_merged: bool = df.item(-1, "type").endswith("_merged")
        df_fast = df.slice(0, -1) if has_ig_merged else df
        seq_len = len(fasta)

        log.trace(f"\tseqid: {seqid}, seq_len: {seq_len}, has_ig_merged: {has_ig_merged}")

        start: int
        end: int
        type_val: str

        for start, end, type_val in df_fast.select(["start", "end", "type"]).iter_rows():
            records.append(
                SeqRecord(
                    seq=fasta[start:end],
                    id=f"{type_val}::{seqid}:{start + idx_fix}-{end}",
                    description="",
                )
            )

        if has_ig_merged:
            start = df.item(-1, "start")
            end = df.item(-1, "end")
            type_val = df.item(-1, "type")
            records.append(
                SeqRecord(
                    seq=fasta[start:] + fasta[: end % seq_len],
                    id=f"{type_val}::{seqid}:{start + idx_fix}-{end}",
                    description="",
                )
            )

        SeqIO.write(records, fasta_out, "fasta")

        return True, seqid, log.get_records()
    except Exception as e:
        an_error = seqid if "seqid" in locals() else gff_in.name
        error_msg = traceback.format_exc()
        log.error(f"Error in biopython_getfasta {an_error}: {e} | {error_msg}")
        return False, an_error, log.get_records()


def gefasta_multiple(
    log: log_setup.GDTLogger,
    tsv_path: Path,
    workers: int,
    gff_in_ext: str = ".gff3",
    gff_in_suffix: str = "_intergenic",
    fasta_in_ext: str = ".fasta",
    fasta_in_suffix: str = "",
    fasta_out_ext: str = ".fasta",
    fasta_out_suffix: str = "_intergenic",
    an_column: str = "AN",
    bedtools_compatible: bool = False,
    skip_region: bool = False,
    overwrite: bool = False,
) -> None:
    """Extract sequences from GFF3 files using Biopython."""
    gff3_utils._ensure_spawn(log)  # Ensure spawn method is set for multiprocessing
    tsv = pl.read_csv(tsv_path, separator="\t")

    gff_in_builder = gff3_utils.PathBuilder(gff_in_ext).use_folder_builder(
        tsv_path.parent, gff_in_suffix
    )
    fasta_in_builder = gff3_utils.PathBuilder(fasta_in_ext).use_folder_builder(
        tsv_path.parent, fasta_in_suffix
    )
    fasta_out_builder = gff3_utils.PathBuilder(fasta_out_ext).use_folder_builder(
        tsv_path.parent, fasta_out_suffix
    )

    gff3_utils.check_files(log, tsv, gff_in_builder, an_column, should_exist=True)

    gff3_utils.check_files(log, tsv, fasta_in_builder, an_column, should_exist=True)

    if not overwrite:
        gff3_utils.check_files(log, tsv, fasta_out_builder, an_column, should_exist=False)

    log.info(f"Starting processing {tsv.shape[0]} ANs with {workers} workers...")
    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        tasks = [
            executor.submit(
                getfasta_single,
                log.spawn_buffer(),
                gff_in_builder.build(an),
                fasta_in_builder.build(an),
                fasta_out_builder.build(an),
                bedtools_compatible,
                skip_region,
            )
            for an in tsv[an_column]
        ]

        for future in cf.as_completed(tasks):
            success, an, records = future.result()
            if not success:
                log.error(f"[{an}] Failed to process.")

            for record in records:
                log.log(record[0], record[1])

        log.info("All processed.")
