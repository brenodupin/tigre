# -*- coding: utf-8 -*-
"""Utilities for working with GFF3 files in the IGRE package.

The base of this module is based on the GDT package, with some modifications
to better suit the needs of IGRE.

Source
https://github.com/brenodupin/gdt/blob/master/src/gdt/gff3_utils.py
"""

import re
from pathlib import Path
from typing import Callable, Optional, Union

import pandas as pd

GFF3_COLUMNS: tuple[str, ...] = (
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
)
QS_GENE = "type == 'gene'"
QS_GENE_TRNA_RRNA = "type in ('gene', 'tRNA', 'rRNA')"

_RE_ID = re.compile(r"ID=([^;]+)")
_RE_dbxref_GeneID = re.compile(r"Dbxref=.*GeneID:")


def load_gff3(
    filename: Union[str, Path],
    sep: str = "\t",
    comment: str = "#",
    header: Optional[int] = None,
    names: tuple[str, ...] = GFF3_COLUMNS,
    usecols: tuple[str, ...] = ("type", "start", "end", "attributes"),
    query_string: Optional[str] = None,
) -> pd.DataFrame:
    """Load a GFF3 file into a pandas DataFrame, optionally filtering by a query string.

    Args:
        filename (Union[str, Path]): Path to the GFF3 file.
        sep (str): Separator used in the file.
        comment (str): Comment character in the file.
        header (int or None): Row number to use as the column names, None if no header.
        names (tuple[str, ...]): Tuple of column names to use.
        usecols (list[str]): List of columns to read from the file.
        query_string (str or None): Query string to filter the DataFrame.

    Returns:
        pd.DataFrame: DataFrame containing the filtered GFF3 data.

    """
    if query_string:
        return (
            pd.read_csv(
                filename,
                sep=sep,
                comment=comment,
                header=header,
                names=names,
                usecols=usecols,
            )
            .query(query_string)
            .sort_values(
                by=["start", "end"], ascending=[True, False], ignore_index=True
            )
        )

    return pd.read_csv(
        filename, sep=sep, comment=comment, header=header, names=names, usecols=usecols
    ).sort_values(by=["start", "end"], ascending=[True, False], ignore_index=True)


def filter_orfs(
    gff3_df: pd.DataFrame,
    orfs_strings: list[str] = ["Name=ORF", "Name=orf"],
) -> pd.DataFrame:
    """Filter out ORFs from a GFF3 DataFrame.

    Args:
        gff3_df (pd.DataFrame): DataFrame containing GFF3 data.
        orfs_strings (list): List of strings to identify ORFs.

    Returns:
        pd.DataFrame: DataFrame with ORFs removed.

    """
    return gff3_df[
        ~gff3_df["attributes"].str.contains("|".join(orfs_strings))
    ].reset_index(drop=True)


class GFFPathBuilder:
    """Flexible class to build GFF3 file paths based on accession numbers.

    This class provides a framework for building file paths with
    different strategies that can be swapped at runtime.
    """

    def __init__(self) -> None:
        """Initialize the builder with no build method set."""
        self._build_method: Optional[Callable[[str], Path]] = None
        self._str: str = "GFFPathBuilder(build=None)"

    def build(self, an: str) -> Path:
        """Build the path for a given accession number.

        This method uses the currently active build strategy.

        Args:
            an (str): Accession number.

        Returns:
            Path: Path to the GFF3 file.

        Raises:
            ValueError: If no build method has been set.

        """
        if self._build_method is None:
            raise ValueError(
                "No build method set. Use one of the 'use_*_builder' methods first."
            )
        return self._build_method(an)

    def use_standard_builder(
        self,
        base: Union[str, Path],
        ext: str = ".gff3",
        suffix: str = "",
    ) -> "GFFPathBuilder":
        """Set the standard filesystem-based build method.

        This method expects files to be stored in a flat structure
        under the specified base directory, with filenames formatted as
        "<accession_number><suffix><ext>".

        Args:
            base (Union[str, Path]): Base directory where files are stored.
            ext (str): File extension for files. Defaults to ".gff3".
            suffix (str): Suffix to append to the filename. Defaults to "".

        Returns:
            GFFPathBuilder: Returns self.

        """
        base_path = Path(base).resolve()

        def standard_builder(an: str) -> Path:
            return base_path / f"{an}{suffix}{ext}"

        self._build_method = standard_builder
        self._str = (
            f"GFFPathBuilder(build='standard', base='{base_path}', ext='{ext}', "
            f"suffix='{suffix}')"
        )
        return self

    def use_folder_builder(
        self,
        base: Union[str, Path],
        ext: str = ".gff3",
        suffix: str = "",
    ) -> "GFFPathBuilder":
        """Set the folder-based build method.

        This method expects files to be stored in subfolders named after
        the accession number, with filenames formatted as
        "<accession_number><suffix><ext>".

        Args:
            base (Union[str, Path]): Base directory where files are stored.
            ext (str): File extension for files. Defaults to ".gff3".
            suffix (str): Suffix to append to the filename. Defaults to "".

        Returns:
            GFFPathBuilder: Returns self.

        """
        base_path = Path(base).resolve()

        def folder_builder(an: str) -> Path:
            return base_path / an / f"{an}{suffix}{ext}"

        self._build_method = folder_builder
        self._str = (
            f"GFFPathBuilder(build='folder', base='{base_path}', ext='{ext}', "
            f"suffix='{suffix}')"
        )
        return self

    def use_custom_builder(
        self,
        builder_func: Callable[[str], Path],
        help_text: Optional[str] = None,
    ) -> "GFFPathBuilder":
        """Set a custom build function as a drop-in replacement.

        Args:
            builder_func (Callable[[str], Path]): A function that takes an
                accession number (str) and returns a Path object.
            help_text (Optional[str]): Optional help text to describe the
                custom builder function. Defaults to None, meaning it will print
                the function name and its signature.

        Returns:
            GFFPathBuilder: Returns self.

        Example:
            def my_custom_builder(an: str) -> Path:
                return Path(f"/custom/path/{an}_special.gff")

            gff_path = GFFPathBuilder("/base/dir")
            gff_path.use_custom_builder(my_custom_builder)
            path = gff_path.build("NC_123456.1")

        """
        self._build_method = builder_func

        if help_text is None:
            help_text = getattr(builder_func, "__name__", "anonymous_function")

        self._str = f"GFFPathBuilder(build='custom', help_text='{help_text}')"
        return self

    def __repr__(self) -> str:
        """Return a string representation of the GFFPathBuilder."""
        return self._str
