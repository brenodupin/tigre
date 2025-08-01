# -*- coding: utf-8 -*-
"""Utilities for working with Fasta files in the tigre package."""

import concurrent.futures as cf
from pathlib import Path
from typing import cast

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .. import gff3_utils, log_setup


def biopython_getfasta(
    log: log_setup.TempLogger,
    gff_in: Path,
    fasta_in: Path,
    fasta_out: Path,
    bedtools_compatible: bool = False,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    seqid = "Not obtained yet"
    try:
        log.trace(f"biopython_getfasta: {gff_in = } | {fasta_in = } | {fasta_out = }")
        log.debug(
            f"biopython_getfasta: gff_in: {gff_in.name} | fasta_in: {fasta_in.name} | fasta_out: {fasta_out.name}"
        )

        idx_fix = 0 if bedtools_compatible else 1
        fasta = SeqIO.read(fasta_in, "fasta").seq  # type: ignore[no-untyped-call]
        records = []

        df = gff3_utils.load_gff3(
            gff_in,  # idx 0      1       2       3
            usecols=("seqid", "type", "start", "end"),
            dtype={"start": "int64", "end": "int64"},
        )

        seqid = df.iat[0, 0]
        df["start"] = df["start"] - 1
        has_ig_merged = df.iat[-1, 1] == "intergenic_region_merged"
        df_fast = df.iloc[:-1] if has_ig_merged else df
        seq_len = len(fasta)

        log.trace(
            f"\tseqid: {seqid}, seq_len: {seq_len}, has_ig_merged: {has_ig_merged}"
        )

        for row in df_fast.itertuples():
            start = cast(int, row.start)
            end = cast(int, row.end)
            records.append(
                SeqRecord(
                    seq=fasta[start:end],
                    id=f"{row.type}::{seqid}:{start + idx_fix}-{end}",
                    description="",
                )
            )

        if has_ig_merged:
            start = cast(int, df.iat[-1, 2])
            end = cast(int, df.iat[-1, 3])
            records.append(
                SeqRecord(
                    seq=fasta[start:] + fasta[: end % seq_len],
                    id=f"intergenic_region_merged::{seqid}:{start + idx_fix}-{end}",
                    description="",
                )
            )

        SeqIO.write(records, fasta_out, "fasta")

        return True, seqid, log.get_records()
    except Exception as e:
        log.error(f"Error in biopython_getfasta: {e}")
        return False, seqid, log.get_records()


def biopython_multiple(
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
    overwrite: bool = False,
) -> None:
    """Wrapper function to execute bedtools getfasta."""
    tsv = pd.read_csv(tsv_path, sep="\t")

    gff_in_builder = gff3_utils.PathBuilder(gff_in_ext).use_folder_builder(
        tsv_path.parent, gff_in_suffix
    )
    fasta_in_builder = gff3_utils.PathBuilder(fasta_in_ext).use_folder_builder(
        tsv_path.parent, fasta_in_suffix
    )
    fasta_out_builder = gff3_utils.PathBuilder(fasta_out_ext).use_folder_builder(
        tsv_path.parent, fasta_out_suffix
    )
    gff3_utils.check_tsv(
        log, tsv, gff_in_builder, fasta_out_builder, overwrite, "AN", "GFF3", "FASTA"
    )

    gff3_utils.check_tsv(
        log,
        tsv,
        fasta_in_builder,
        fasta_out_builder,
        overwrite,
        an_column,
        "FASTA",
        "FASTA",
    )

    log.info(f"Starting processing {tsv.shape[0]} ANs with {workers} workers...")
    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        tasks = [
            executor.submit(
                biopython_getfasta,
                log.spawn_buffer(),
                gff_in_builder.build(an),
                fasta_in_builder.build(an),
                fasta_out_builder.build(an),
                bedtools_compatible,
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
