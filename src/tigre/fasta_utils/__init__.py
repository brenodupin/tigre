# -*- coding: utf-8 -*-
"""Fasta utilities submodule with optional biopython dependency."""

from typing import TYPE_CHECKING

# waiting for bedtools to merge our PR adding support for circular sequences
# https://github.com/arq5x/bedtools2/pull/1127
# from .bedtools_wrapper import bedtools_getfasta, bedtools_multiple, bedtools_version


if TYPE_CHECKING:
    from pathlib import Path

    from .. import log_setup

try:
    # utilities that require biopython
    from .biopython_wrapper import biopython_getfasta, biopython_multiple

    BIOPYTHON_AVAILABLE = True


except ImportError:
    BIOPYTHON_AVAILABLE = False

    def biopython_getfasta(
        log: "log_setup.TempLogger",
        gff_in: "Path",
        fasta_in: "Path",
        fasta_out: "Path",
        bedtools_compatible: bool = False,
        skip_region: bool = False,
    ) -> tuple[bool, str, list[tuple[int, str]]]:
        """Stub for biopython_getfasta when Biopython is not available."""
        raise ImportError(
            "Biopython is required for `biopython_getfasta`. "
            "Install it with: pip install biopython or pip install tigre[bio]"
        )

    def biopython_multiple(
        log: "log_setup.GDTLogger",
        tsv_path: "Path",
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
        """Stub for biopython_multiple when Biopython is not available."""
        raise ImportError(
            "Biopython is required for `biopython_multiple`. "
            "Install it with: pip install biopython or pip install tigre[bio]"
        )


__all__ = [
    "bedtools_getfasta",
    "bedtools_multiple",
    "biopython_getfasta",
    "biopython_multiple",
    "get_bedtools_version",
    "BIOPYTHON_AVAILABLE",
]
