# -*- coding: utf-8 -*-
"""Bedtools wrapper functions for `bedtools getfasta``."""

import subprocess
from pathlib import Path
from typing import Union

from .. import log_setup


def check_bedtools_available(bedtools_path: Union[Path, str] = "bedtools") -> bool:
    """Check if bedtools is available in the system."""
    try:
        result = subprocess.run(
            [bedtools_path, "--version"], capture_output=True, text=True, timeout=5
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


def bedtools_getfasta(
    log: log_setup.GDTLogger,
    gff_in: Path,
    fasta_in: Path,
    fasta_out: Path,
    bedtools_path: str = "bedtools",
) -> None:
    """Extract sequences from a FASTA file using bedtools.

    Args:
        gff_in (Path): Path to the GFF3 file containing the regions of interest.
        fasta_in (Path): Path to the input FASTA file.
        fasta_out (Path): Path to the output FASTA file where sequences will be saved.
        bedtools_path (str): Path to the bedtools executable. Defaults to "bedtools".

    """
    command = (
        f"{bedtools_path} getfasta -name+ -fi {fasta_in} -fo {fasta_out} -bed {gff_in}"
    )
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
