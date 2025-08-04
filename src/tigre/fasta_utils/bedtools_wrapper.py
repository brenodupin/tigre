# -*- coding: utf-8 -*-
"""Bedtools wrapper functions for `bedtools getfasta``."""

import concurrent.futures as cf
import subprocess
from pathlib import Path
from typing import Optional

import pandas as pd

from .. import gff3_utils, log_setup


def get_bedtools_version(bedtools_path: Path | str = "bedtools") -> Optional[str]:
    """Check if bedtools is available in the system."""
    try:
        result = subprocess.run(
            [bedtools_path, "--version"], capture_output=True, text=True, timeout=5
        )
        return result.stdout.strip() if result.returncode == 0 else None
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None


def bedtools_getfasta(
    log: log_setup.TempLogger,
    gff_in: Path,
    fasta_in: Path,
    fasta_out: Path,
    bedtools_path: str = "bedtools",
    name_args: str = "-name+",
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Extract sequences from a FASTA file using bedtools.

    Args:
        gff_in (Path): Path to the GFF3 file containing the regions of interest.
        fasta_in (Path): Path to the input FASTA file.
        fasta_out (Path): Path to the output FASTA file where sequences will be saved.
        bedtools_path (str): Path to the bedtools executable. Defaults to "bedtools".

    """
    try:
        command = f"{bedtools_path} getfasta {name_args} -fi {fasta_in} -fo {fasta_out} -bed {gff_in}"
        log.debug(f"Running command: {command}")
        output = subprocess.run(
            command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
        )
        if output.returncode != 0:
            log.error(f"Error running bedtools getfasta: {output.stderr}")
            raise RuntimeError(f"bedtools getfasta failed with error: {output.stderr}")

        log.info(f"Sequences extracted to {fasta_out} using bedtools getfasta.")
        for line in output.stdout.splitlines():
            if line := line.strip():
                log.trace(f"b: {line}")
        return True, gff_in.stem, log.get_records()

    except Exception as e:
        log.error(f"Error in bedtools_getfasta: {e}")
        return False, gff_in.stem, log.get_records()


def bedtools_multiple(
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
    bedtools_path: str = "bedtools",
    name_args: str = "-name+",
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

    gff3_utils.check_files(log, tsv, gff_in_builder, an_column, should_exist=True)

    gff3_utils.check_files(log, tsv, fasta_in_builder, an_column, should_exist=True)

    if not overwrite:
        gff3_utils.check_files(
            log, tsv, fasta_out_builder, an_column, should_exist=False
        )

    log.trace(
        "hard stop to check some things:"
        f"builders: {gff_in_builder = } | {fasta_in_builder = } | {fasta_out_builder = }"
        f"args: {tsv_path = } | {gff_in_ext = } | {gff_in_suffix = }\n"
        f"{fasta_in_ext = } | {fasta_in_suffix = } | {fasta_out_ext = }\n"
        f"{fasta_out_suffix = } | {an_column = } | {workers = }\n"
        f"{bedtools_path = } | {name_args = } | {overwrite = }"
    )
    raise NotImplementedError(
        "bedtools still does not support circular sequences, we've submitted a PR to add this feature, "
        "but it is not yet merged. "
        "Please use the biopython process for now."
    )

    log.info(f"Starting processing {tsv.shape[0]} ANs with {workers} workers...")
    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        tasks = [
            executor.submit(
                bedtools_getfasta,
                log.spawn_buffer(),
                gff_in_builder.build(an),
                fasta_in_builder.build(an),
                fasta_out_builder.build(an),
                bedtools_path,
                name_args,
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
