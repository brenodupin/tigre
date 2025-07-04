# -*- coding: utf-8 -*-
"""Utilities for working with GFF3 files in the IGRE package.

The base of this module is based on the GDT package, with some modifications
to better suit the needs of IGRE.

Source
https://github.com/brenodupin/gdt/blob/master/src/gdt/gff3_utils.py
"""

import re
from pathlib import Path
from typing import Optional, Union

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
    usecols: list[str] = ["type", "start", "end", "attributes"],
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
