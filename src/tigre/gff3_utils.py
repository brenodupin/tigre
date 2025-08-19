# -*- coding: utf-8 -*-
"""Utilities for working with GFF3 files in the tigre package.

The base of this module is based on the GDT package, with some modifications
to better suit the needs of tigre.

Source
https://github.com/brenodupin/gdt/blob/master/src/gdt/gff3_utils.py
"""

import re
import sys
from functools import partial
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from . import log_setup

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
QS_GENE_REGION = "type in ('gene', 'region')"
QS_GENE_TRNA_RRNA_REGION = "type in ('gene', 'tRNA', 'rRNA', 'region')"

_RE_ID = re.compile(r"ID=([^;]+)")
_RE_region_taxon = re.compile(r"taxon:([^;,]+)")


def load_gff3(
    filename: str | Path,
    sep: str = "\t",
    comment: str = "#",
    header: int | None = None,
    names: tuple[str, ...] = GFF3_COLUMNS,
    usecols: tuple[str, ...] = ("type", "start", "end", "attributes"),
    query_string: str | None = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """Load a GFF3 file into a pandas DataFrame, optionally filtering by a query string.

    Args:
        filename (str | Path): Path to the GFF3 file.
        sep (str): Separator used in the file.
        comment (str): Comment character in the file.
        header (int | None): Row number to use as the column names, None if no header.
        names (tuple[str, ...]): Tuple of column names to use.
        usecols (list[str]): List of columns to read from the file.
        query_string (str | None): Query string to filter the DataFrame.
                                   Defaults to None, which means no filtering.
        **kwargs: Additional keyword arguments passed to `pd.read_csv`.

    Returns:
        pd.DataFrame: DataFrame containing the filtered GFF3 data.

    """
    if query_string:
        df: pd.DataFrame = (
            pd.read_csv(
                filename,
                sep=sep,
                comment=comment,
                header=header,
                names=names,
                usecols=usecols,
                **kwargs,
            )
            .query(query_string)
            .sort_values(
                by=["start", "end"], ascending=[True, False], ignore_index=True
            )
        )
        return df

    df = pd.read_csv(
        filename,
        sep=sep,
        comment=comment,
        header=header,
        names=names,
        usecols=usecols,
        **kwargs,
    ).sort_values(by=["start", "end"], ascending=[True, False], ignore_index=True)
    return df


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


# will be pre-compiled by the PathBuilder class
def _std_builder(
    base_path: Path,
    suffix: str,
    ext: str,
    an: str,
) -> Path:
    return base_path / f"{an}{suffix}{ext}"


# will be pre-compiled by the PathBuilder class
def _folder_builder(
    base_path: Path,
    suffix: str,
    ext: str,
    an: str,
) -> Path:
    return base_path / an / f"{an}{suffix}{ext}"


class PathBuilder:
    """Flexible class to build paths (mostly GFF/Fasta) based on accession numbers.

    This class provides a framework for building file paths with
    different strategies that can be swapped at runtime.
    """

    def __init__(self, ext: str = ".gff3") -> None:
        """Initialize the builder with no build method set."""
        self._build_method: Callable[[str], Path] | None = None
        self.ext: str = ext
        self._str: str = f"PathBuilder('{ext}', build=None)"

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

    def use_std_builder(
        self,
        base: str | Path,
        suffix: str = "",
    ) -> "PathBuilder":
        """Set the standard filesystem-based build method.

        This method expects files to be stored in a flat structure
        under the specified base directory, with filenames formatted as
        "<accession_number><suffix><self.ext>".

        Args:
            base (str | Path]): Base directory where files are stored.
            suffix (str): Suffix to append to the filename. Defaults to "".

        Returns:
            PathBuilder: Returns self.

        """
        base_path = Path(base).resolve()

        self._build_method = partial(
            _std_builder,
            base_path,
            suffix,
            self.ext,
        )
        self._str = (
            f"PathBuilder('{self.ext}', build='std', base='{base_path}', "
            f"suffix='{suffix}')"
        )
        return self

    def use_folder_builder(
        self,
        base: str | Path,
        suffix: str = "",
    ) -> "PathBuilder":
        """Set the folder-based build method.

        This method expects files to be stored in subfolders named after
        the accession number, with filenames formatted as
        "<accession_number>/<accession_number><suffix><self.ext>".

        Args:
            base (str | Path): Base directory where files are stored.
            suffix (str): Suffix to append to the filename. Defaults to "".

        Returns:
            PathBuilder: Returns self.

        """
        base_path = Path(base).resolve()

        self._build_method = partial(
            _folder_builder,
            base_path,
            suffix,
            self.ext,
        )

        self._str = (
            f"PathBuilder('{self.ext}', build='folder', base='{base_path}', "
            f"suffix='{suffix}')"
        )
        return self

    def use_custom_builder(
        self,
        builder_func: Callable[[str], Path],
        help_text: str | None = None,
    ) -> "PathBuilder":
        """Set a custom build function as a drop-in replacement.

        It should return a Path object based on the given accession number.
        The extension is a attribute of the class (defined at __init__), so it
        can be used by the custom builder function, just access it with `self.ext`.

        Args:
            builder_func (Callable[[str], Path]): A function that takes an
                accession number (str) and returns a Path object.
            help_text (str | None): Optional help text to describe the
                custom builder function. This is used by the __repr__ method
                to provide more context about the builder. If not provided,
                it defaults to the function's name.

        Returns:
            PathBuilder: Returns self.

        Example:
            def my_custom_builder(an: str) -> Path:
                return Path(f"/custom/path/{an}_special{self.ext}")

            gff_path = PathBuilder(".gff")
            gff_path.use_custom_builder(my_custom_builder)
            path = gff_path.build("NC_123456.1")
            print(path)  # Outputs: /custom/path/NC_123456.1_special.gff

        """
        self._build_method = builder_func

        if help_text is None:
            help_text = getattr(builder_func, "__name__", "anonymous_function")

        self._str = (
            f"PathBuilder('{self.ext}', build='custom', help_text='{help_text}')"
        )
        return self

    def __repr__(self) -> str:
        """Return a string representation of the PathBuilder."""
        return self._str


def check_files(
    log: log_setup.GDTLogger,
    df: pd.DataFrame,
    file_builder: PathBuilder,
    an_column: str = "AN",
    *,
    should_exist: bool = True,
) -> None:
    """Check files based on the DataFrame and the provided PathBuilder.

    Args:
        log: Logger instance
        df: DataFrame containing the data
        file_builder: PathBuilder instance for constructing file paths
        check_for_existence: If True, check that files exist (input check).
                           If False, check that files don't exist (output check).
        an_column: Column name containing the identifiers
        should_exist: If True, checks if files exist; if False,
                      checks if they do not exist.

    """
    log.trace(
        f"check_files called | builder: {file_builder} | an_column: {an_column} | "
        f"should_exist: {should_exist}"
    )

    problem_files = []
    check_type = "Missing" if should_exist else "Existing"

    for an in df[an_column]:
        file_path = file_builder.build(an)

        if file_path.is_file() != should_exist:
            problem_files.append((an, file_path))

    if problem_files:
        for an, path in problem_files:
            log.trace(f"{check_type} file for {an}, path: {path}")

        if should_exist:
            log.error(
                f"{check_type} {len(problem_files)} file(s). Please check the log for "
                "details."
            )
        else:
            log.error(
                f"{check_type} {len(problem_files)} file(s). Use `--overwrite` to "
                "overwrite them."
            )

        sys.exit(1)
