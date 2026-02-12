# -*- coding: utf-8 -*-
"""Module to combine 2 GFF3 files into one."""

import concurrent.futures as cf
from pathlib import Path

import pandas as pd

from . import gff3_utils, log_setup


def combine_pair(
    log: log_setup.TempLogger,
    gff_file_1: Path,
    gff_file_2: Path,
    gff_out: Path,
) -> tuple[bool, str, list[log_setup._RawMsg]]:
    """Combine two GFF3 files into one.

    This combination follows the procedure:
    1. Read both GFF3 files as DataFrames.
    2. Copy region line (first row of GFF3) from the `gff-in-1` file.
    3. Remove region lines (`type == 'region'`) from both DataFrames.
    4. Check for duplicated IDs between both files and log a warning if any is found.
    5. Concatenate both DataFrames, with the region line extracted earlier as the first row.
    6. Write the combined DataFrame to a `gff-out` path, preserving the header from the
    `gff-in-1` file.

    Args:
        log: Logger instance
        gff_file_1: Path to the first GFF3 file
        gff_file_2: Path to the second GFF3 file
        gff_out: Path to the output GFF3 file

    """
    try:
        header = []
        with open(gff_file_1, "r") as gff_file:
            for line in gff_file:
                if not line.startswith("#"):
                    break

                header.append(line.strip())

            df1 = gff3_utils.load_gff3(
                gff_file_1,
                usecols=gff3_utils.GFF3_COLUMNS,
            )

        df2 = gff3_utils.load_gff3(
            gff_file_2,
            usecols=gff3_utils.GFF3_COLUMNS,
        )

        # region line should be the first line in the GFF3 file
        region = df1.iloc[0].copy()

        df1 = df1[df1["type"] != "region"]
        df2 = df2[df2["type"] != "region"]

        df1_gene_ids: "pd.Series[str]" = (
            df1["attributes"]
            .str.extract(gff3_utils._RE_ID, expand=False)  # type: ignore[call-overload]
            .astype("string")
        )
        df2_gene_ids: "pd.Series[str]" = (
            df2["attributes"]
            .str.extract(gff3_utils._RE_ID, expand=False)  # type: ignore[call-overload]
            .astype("string")
        )

        duplicated_ids = set(df1_gene_ids).intersection(df2_gene_ids)
        if duplicated_ids:
            log.warning(
                f"Warning: Found {len(duplicated_ids)} duplicated IDs between input "
                "files. The combined GFF3 file may have non-unique IDs.",
            )
            log.debug(f"Duplicated IDs: {', '.join(duplicated_ids)}")

        combined_df = pd.concat([pd.DataFrame([region]), df1, df2], ignore_index=True)

        with open(gff_out, "w") as f:
            f.write("\n".join(header) + "\n")
            combined_df.to_csv(f, sep="\t", index=False, header=False)

        return True, gff_out.stem, log.get_records()

    except Exception as e:
        log.error(
            f"Error processing files '{gff_file_1.name}' and '{gff_file_2.name}': {e}"
        )
        return False, gff_out.stem, log.get_records()


def combine_multiple(
    log: log_setup.GDTLogger,
    tsv_path: Path,
    workers: int,
    gff_in_ext_1: str = ".gff3",
    gff_in_suffix_1: str = "",
    gff_in_ext_2: str = ".gff3",
    gff_in_suffix_2: str = "",
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

        gff_in_ext_1: Extension for the first input GFF3 files
        gff_in_suffix_1: Suffix for the first input GFF3 files

        gff_in_ext_2: Extension for the second input GFF3 files
        gff_in_suffix_2: Suffix for the second input GFF3 files

        gff_out_ext: Extension for the output GFF3 files
        gff_out_suffix: Suffix for the output GFF3 files

        an_column: Column name in the TSV file that contains ANs
        overwrite: Whether to overwrite existing output files

    """
    gff3_utils._ensure_spawn(log)
    tsv = pd.read_csv(tsv_path, sep="\t")

    gff_in_builder_1 = gff3_utils.PathBuilder(gff_in_ext_1).use_folder_builder(
        tsv_path.parent,
        gff_in_suffix_1,
    )
    gff_in_builder_2 = gff3_utils.PathBuilder(gff_in_ext_2).use_folder_builder(
        tsv_path.parent,
        gff_in_suffix_2,
    )
    gff_out_builder = gff3_utils.PathBuilder(gff_out_ext).use_folder_builder(
        tsv_path.parent,
        gff_out_suffix,
    )

    gff3_utils.check_files(
        log,
        tsv,
        gff_in_builder_1,
        an_column,
        should_exist=True,
    )
    gff3_utils.check_files(
        log,
        tsv,
        gff_in_builder_2,
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
                gff_in_builder_1.build(an),
                gff_in_builder_2.build(an),
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
