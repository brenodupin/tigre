# -*- coding: utf-8 -*-
"""Module to combine 2 GFF3 files into one."""

import concurrent.futures as cf
from pathlib import Path

import pandas as pd

from . import gff3_utils, log_setup


def combine_pair(
    log: log_setup.TempLogger,
    gff3_file1: Path,
    gff3_file2: Path,
    gff_out: Path,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Combine two GFF3 files into one, does not ensure uniqueness of IDs (yet)."""
    try:
        header = []
        with open(gff3_file1, "r") as gff_file:
            for line in gff_file:
                if not line.startswith("#"):
                    break

                header.append(line.strip())

            df1 = gff3_utils.load_gff3(
                gff3_file1,
                usecols=gff3_utils.GFF3_COLUMNS,
            )

        df2 = gff3_utils.load_gff3(
            gff3_file2,
            usecols=gff3_utils.GFF3_COLUMNS,
        )
        # remove region from df2 (first element)
        df2 = df2.loc[1:]

        df1["gene_id"] = df1["attributes"].str.extract(gff3_utils._RE_ID, expand=False)  # type: ignore[call-overload]
        df2["gene_id"] = df2["attributes"].str.extract(gff3_utils._RE_ID, expand=False)  # type: ignore[call-overload]

        duplicated_ids = set(df1["gene_id"]).intersection(set(df2["gene_id"]))
        if duplicated_ids:
            log.warning(
                f"Warning: Found {len(duplicated_ids)} duplicated IDs between input "
                "files. The combined GFF3 file may have non-unique IDs.",
            )
            log.debug(f"Duplicated IDs: {', '.join(duplicated_ids)}")

        df1, df2 = df1.drop(columns=["gene_id"]), df2.drop(columns=["gene_id"])

        combined_df = pd.concat([df1, df2], ignore_index=True)

        with open(gff_out, "w") as f:
            f.write("\n".join(header) + "\n")
            combined_df.to_csv(f, sep="\t", index=False, header=False)

        return True, gff_out.stem, log.get_records()

    except Exception as e:
        log.error(f"Error processing files '{gff3_file1}' and '{gff3_file2}': {e}")
        return False, gff_out.stem, log.get_records()


def combine_multiple(
    log: log_setup.GDTLogger,
    tsv_path: Path,
    workers: int,
    gff1_in_ext: str = ".gff3",
    gff1_in_suffix: str = "",
    gff2_in_ext: str = ".gff3",
    gff2_in_suffix: str = "",
    gff_out_ext: str = ".gff3",
    gff_out_suffix: str = "_clean",
    an_column: str = "AN",
    overwrite: bool = False,
) -> None:
    """Orchestrates the execution of `combine_pair` across multiple GFF3 files.

    Args:
        log: Logger instance
        tsv_path: Path to the TSV file containing ANs
        workers: Number of worker processes to use
        gff1_in_ext: Extension for the first input GFF3 files
        gff1_in_suffix: Suffix for the first input GFF3 files
        gff2_in_ext: Extension for the second input GFF3 files
        gff2_in_suffix: Suffix for the second input GFF3 files
        gff_out_ext: Extension for the output GFF3 files
        gff_out_suffix: Suffix for the output GFF3 files
        an_column: Column name in the TSV file that contains ANs
        overwrite: Whether to overwrite existing output files

    """
    tsv = pd.read_csv(tsv_path, sep="\t")

    gff1_in_builder = gff3_utils.PathBuilder(gff1_in_ext).use_folder_builder(
        tsv_path.parent,
        gff1_in_suffix,
    )
    gff2_in_builder = gff3_utils.PathBuilder(gff2_in_ext).use_folder_builder(
        tsv_path.parent,
        gff2_in_suffix,
    )
    gff_out_builder = gff3_utils.PathBuilder(gff_out_ext).use_folder_builder(
        tsv_path.parent,
        gff_out_suffix,
    )

    gff3_utils.check_files(
        log,
        tsv,
        gff1_in_builder,
        an_column,
        should_exist=True,
    )
    gff3_utils.check_files(
        log,
        tsv,
        gff2_in_builder,
        an_column,
        should_exist=True,
    )

    if not overwrite:
        gff3_utils.check_files(
            log,
            tsv,
            gff_out_builder,
            an_column,
            should_exist=False,
        )

    log.info(f"Starting processing {tsv.shape[0]} ANs with {workers} workers...")
    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        tasks = [
            executor.submit(
                combine_pair,
                log.spawn_buffer(),
                gff1_in_builder.build(an),
                gff2_in_builder.build(an),
                gff_out_builder.build(an),
            )
            for an in tsv[an_column]
        ]
        log.trace("Tasks submitted, waiting for completion...")

        for future in cf.as_completed(tasks):
            success, an, records = future.result()
            if not success:
                log.error(f"[{an}] Failed to process.")

            for record in records:
                log.log(record[0], record[1])

        log.info("All processed.")
