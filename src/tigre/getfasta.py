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


def _build_description(
    seqid: str,
    type_val: str,
    start: int,
    end: int,
    length: int,
    zero_based: bool,
) -> str:
    """Build the space-separated key=value FASTA description string."""
    description = (
        f"seqid={seqid} type={type_val} start={start + 1} end={end} length={length}"
    )
    if zero_based:
        description += f" start_0={start} end_0={end}"

    return description


def getfasta_single(
    log: log_setup.TempLogger,
    gff_in: Path,
    fasta_in: Path,
    fasta_out: Path,
    zero_based: bool = False,
    skip_region: bool = False,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Extract sequences from a FASTA file using Biopython."""
    try:
        log.trace(f"biopython_getfasta: {gff_in = } | {fasta_in = } | {fasta_out = }")
        log.debug(
            f"biopython_getfasta: gff_in: {gff_in.name} | fasta_in: {fasta_in.name} | "
            f"fasta_out: {fasta_out.name}"
        )

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
                pl.col("attributes").str.extract(gff3_utils._RS_ID, 1).alias("gff_id"),
            ]
        )

        if skip_region:
            df = df.filter(pl.col("type") != "region")

        n_rows = df.height
        df = df.filter(pl.col("end") >= pl.col("start"))
        if (n_zero := n_rows - df.height) > 0:
            log.info(f"Skipped {n_zero} zero-length region(s) in {gff_in.name}")

        if df.is_empty():
            log.warning(f"Empty GFF3 file: {gff_in.name}, creating empty FASTA.")
            SeqIO.write([], fasta_out, "fasta")
            return True, gff_in.name, log.get_records()

        if df["gff_id"].is_null().any():
            log.error(
                f"A feature in {gff_in.name} lacks the `ID` field, creating empty FASTA."
            )

            for attrs in df.filter(pl.col("gff_id").is_null())["attributes"]:
                log.error(f"  att missing ID: {attrs}")

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
        type_: str
        gff_id: str
        view = df_fast.select(["start", "end", "type", "gff_id"])

        for start, end, type_, gff_id in view.iter_rows():
            seq = fasta[start:end]

            if len(gff_id) > 50:
                log.warning(f"FASTA id exceeds 50 characters ({len(gff_id)}): {gff_id}")

            records.append(
                SeqRecord(
                    seq=seq,
                    id=gff_id,
                    description=_build_description(
                        seqid,
                        type_,
                        start,
                        end,
                        len(seq),
                        zero_based,
                    ),
                )
            )

        if has_ig_merged:
            start = df.item(-1, "start")
            end = df.item(-1, "end")
            type_ = df.item(-1, "type")
            gff_id = df.item(-1, "gff_id")
            seq = fasta[start:] + fasta[: end % seq_len]

            if len(gff_id) > 50:
                log.warning(f"FASTA id exceeds 50 characters ({len(gff_id)}): {gff_id}")

            records.append(
                SeqRecord(
                    seq=seq,
                    id=gff_id,
                    description=_build_description(
                        seqid,
                        type_,
                        start,
                        end,
                        len(seq),
                        zero_based,
                    ),
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
    zero_based: bool = False,
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
                zero_based,
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
