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
    an: str,
    gff_in: Path,
    fasta_in: Path,
    fasta_out: Path,
    use_bedtools_index: bool = True,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    try:
        log.debug(f"Extracting sequences from {fasta_in} using Biopython...")

        log.debug(f"an: {an}")
        log.debug(f"gff_in: {gff_in}")
        log.debug(f"fasta_in: {fasta_in}")
        log.debug(f"fasta_out: {fasta_out}")

        idx_fix = 0 if use_bedtools_index else 1
        fasta = SeqIO.read(fasta_in, "fasta").seq  # type: ignore[no-untyped-call]
        records = []

        df = gff3_utils.load_gff3(
            gff_in,  # idx 0      1       2       3
            usecols=("seqid", "start", "end", "type"),
            dtype={"start": "int64", "end": "int64"},
        )

        seqid = df.iat[0, 0]
        df["start"] = df["start"] - 1
        has_ig_merged = df.iat[-1, 3] == "intergenic_region_merged"
        df_fast = df.iloc[:-1] if has_ig_merged else df

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

        return True, an, log.get_records()
    except Exception as e:
        log.error(f"Error in biopython_getfasta: {e}")
        return False, an, log.get_records()


def biopython_multiple(
    log: log_setup.GDTLogger,
    tsv_path: Path,
    gff_in_ext: str = ".gff3",
    gff_in_suffix: str = "_intergenic",
    fasta_in_ext: str = ".fasta",
    fasta_in_suffix: str = "",
    fasta_out_ext: str = ".fasta",
    fasta_out_suffix: str = "_intergenic",
    an_column: str = "AN",
    workers: int = 0,
    use_bedtools_index: bool = True,
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
    gff3_utils.check_file_in_tsv(
        log,
        tsv,
        gff_in_builder,
        an_column,
        "GFF3",
    )

    log.info(f"Starting processing {tsv.shape[0]} ANs with {workers} workers...")
    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        tasks = [
            executor.submit(
                biopython_getfasta,
                log.spawn_buffer(),
                an,
                gff_in_builder.build(an),
                fasta_in_builder.build(an),
                fasta_out_builder.build(an),
                use_bedtools_index,
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
