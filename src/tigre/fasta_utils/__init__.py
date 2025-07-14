# -*- coding: utf-8 -*-
"""Fasta utilities submodule with optional biopython dependency."""

from typing import TYPE_CHECKING

from .bedtools_wrapper import bedtools_getfasta, check_bedtools_available

if TYPE_CHECKING:
    from pathlib import Path

    from .. import log_setup

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    BIOPYTHON_AVAILABLE = True

    # utilities that require biopython
    from .biopython_wrapper import biopython_getfasta


except ImportError:
    BIOPYTHON_AVAILABLE = False

    # stub function if biopython is not available
    def biopython_getfasta(
        log: "log_setup.TempLogger",
        an: str,
        gff_in: "Path",
        fasta_in: "Path",
        fasta_out: "Path",
        use_bedtools_index: bool = True,
    ) -> tuple[bool, str, list[tuple[int, str]]]:
        raise ImportError(
            "Biopython is required for `biopython_getfasta`. "
            "Install it with: pip install biopython or pip install tigre[bio]"
        )


__all__ = [
    "bedtools_getfasta",
    "biopython_getfasta",
    "check_bedtools_available",
    "BIOPYTHON_AVAILABLE",
]
