# -*- coding: utf-8 -*-
"""Utilities for working with Fasta files in the tigre package."""

from pathlib import Path
from typing import cast

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .. import gff3_utils, log_setup


def biopython_getfasta(
    log: log_setup.GDTLogger,
    an: str,
    gff_in: Path,
    fasta_in: Path,
    fasta_out: Path,
    bedtools_idx: bool = False,
) -> None:
    log.debug(f"Extracting sequences from {fasta_in} using Biopython...")

    log.debug(f"an: {an}")
    log.debug(f"gff_in: {gff_in}")
    log.debug(f"fasta_in: {fasta_in}")
    log.debug(f"fasta_out: {fasta_out}")

    idx_fix = 0 if bedtools_idx else 1
    fasta = SeqIO.read(fasta_in, "fasta").seq  # type: ignore[no-untyped-call]
    records = []

    df = gff3_utils.load_gff3(
        gff_in,  # idx 0      1       2       3
        usecols=("seqid", "start", "end", "type"),
        dtype={"start": "int64", "end": "int64"},
    )

    seqid = df.iat[0, 0]
    df["start"] = df["start"] - 1
    has_merge = df.iat[-1, 3] == "intergenic_region_merged"
    df_fast = df.iloc[:-1] if has_merge else df

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

    if has_merge:
        start = cast(int, df.iat[-1, 1])
        end = cast(int, df.iat[-1, 2])
        records.append(
            SeqRecord(
                seq=fasta[start:] + fasta[: end % len(fasta)],
                id=f"intergenic_region_merged::{seqid}:{start + idx_fix}-{end}",
                description="",
            )
        )

    SeqIO.write(records, fasta_out, "fasta")
    log.info(f"Sequences extracted to {fasta_out} using Biopython.")
