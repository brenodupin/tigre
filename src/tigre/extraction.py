# -*- coding: utf-8 -*-
"""Module to extract intergenic regions from GFF3 files."""

import argparse
import subprocess
from pathlib import Path
from shutil import copy2
from typing import Optional
from typing_extensions import Final
import pandas as pd

from . import log_setup
from . import gff3_utils

import concurrent.futures

IGRegion = tuple[int, int, int, str]  # gene_id, start, end, attributes


def change_names(log: log_setup.GDTLogger, att: str, suffix: str) -> str:
    log.trace(f"change_names | att: {att} | suffix: {suffix}")
    idx = att.index("source=")
    return f"name_{suffix}={att[5:idx]}source_{suffix}={att[idx+7:]}"


def get_names(log: log_setup.GDTLogger, att: str, is_left: bool = False) -> str:
    log.trace(f"get_names | att: {att} | is_left: {is_left}")
    if is_left:
        name_start = att.index("name_left=") + 10  # len('name_left=')
        source_start = att.index("source_left=") + 12  # len('source_left=')
    else:
        name_start = att.index("name_right=") + 11  # len('name_right=')
        source_start = att.index("source_right=") + 13  # len('source_right=')

    name_end = att.find(";", name_start)
    source_end = att.find(";", source_start)
    return f"name={att[name_start:name_end]};source={att[source_start:source_end]};"


def get_names_ig(log: log_setup.GDTLogger, att: str, suffix: str) -> str:
    log.trace(f"get_names | att: {att} | suffix: {suffix}")
    name_start = att.index(f"name_{suffix}=")
    source_start = att.index(f"source_{suffix}=")
    name_end = att.find(";", name_start)
    source_end = att.find(";", source_start)
    return f"{att[name_start:name_end]};{att[source_start:source_end]};"


def exec_bedtools(
    log: log_setup.GDTLogger, fasta_in: Path, gff_in: Path, fasta_out: Path
) -> None:
    log.debug(
        f"exec_bedtools | cmd: bedtools getfasta -name+ -fi {fasta_in} -bed {gff_in} -fo {fasta_out}"
    )
    a = subprocess.run(
        [
            "bedtools",
            "getfasta",
            "-name+",
            "-fi",
            fasta_in,
            "-bed",
            gff_in,
            "-fo",
            fasta_out,
        ],
        check=True,
        capture_output=True,
    )
    if st := a.stdout.decode().strip():
        log.debug(f"bedtools stdout: {st}")
    if st := a.stderr.decode().strip():
        log.debug(f"bedtools stderr: {st}")


def save_gff(
    log: log_setup.GDTLogger,
    file_path: Path,
    ig_regions: list[IGRegion],
    seqid: str,
    source: str,
    header: list[str] = [],
) -> None:
    log.trace(f"save_gff | file_path: {file_path}")
    with open(file=file_path, mode="w") as out_file:
        [out_file.write(line) for line in header]
        for region in ig_regions:
            out_file.write(
                f"{seqid}\t{source}\tintergenic_region\t{region[1]}\t{region[2]}\t.\t+\t.\tID={region[0]};{region[3]}\n"
            )


def save_gff_with_merge(
    log: log_setup.GDTLogger,
    file_path: Path,
    ig_regions: list[IGRegion],
    seqid: str,
    source: str,
    header: list[str] = [],
) -> None:
    log.trace(f"save_gff_with_merge | file_path: {file_path}")
    with open(file=file_path, mode="w") as out_file:
        [out_file.write(line) for line in header]
        for region in ig_regions[:-1]:
            out_file.write(
                f"{seqid}\t{source}\tintergenic_region\t{region[1]}\t{region[2]}\t.\t+\t.\tID={region[0]};{region[3]}\n"
            )

        out_file.write(
            f"{seqid}\t{source}\tintergenic_region_merged\t{ig_regions[-1][1]}\t{ig_regions[-1][2]}\t.\t+\t.\tID={ig_regions[-1][0]};{ig_regions[-1][3]}\n"
        )


def save_gff_with_frag(
    log: log_setup.GDTLogger,
    file_path: Path,
    ig_regions: list[IGRegion],
    seqid: str,
    source: str,
    header: list[str] = [],
) -> None:
    log.trace(f"save_gff_with_frag | file_path: {file_path}")
    max_len = len(ig_regions) - 1
    with open(file=file_path, mode="w") as out_file:
        [out_file.write(line) for line in header]
        for i, region in enumerate(ig_regions):
            out_file.write(
                (
                    f"{seqid}\t{source}\t{'intergenic_region' if not (i == 0 or i == max_len) else 'intergenic_region_fragment'}"
                    f"\t{region[1]}\t{region[2]}\t.\t+\t.\tID={region[0]};{region[3]}\n"
                )
            )


# ('seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes')
#     0        1        2       3        4      5         6        7         8
def execution(
    log: log_setup.GDTLogger,
    gff_in: Path,
    gff_out: Path,
    fasta_in: Optional[Path],
    fasta_out: Optional[Path],
    run_bedtools: bool = False,
    new_region_start: int = 1,
) -> None:
    rs, re = False, False
    first_att = ""
    header = []
    if run_bedtools:
        if fasta_in is None or fasta_out is None:
            raise ValueError(
                "Fasta input and output files must be provided when run_bedtools is True."
            )
        fasta_in_path: Path = fasta_in
        fasta_out_path: Path = fasta_out

    with open(gff_in) as gff_file:
        for header_line in gff_file:
            if not header_line.startswith("#"):
                break
            header.append(header_line)

        # region line (first line)
        region_line: list[str] = header_line.strip().split("\t")
        is_circular = "is_circular=true" in region_line[8].lower()
        genome_size = int(region_line[4])
        seqid = region_line[0]
        source = region_line[1]
        log.debug(
            f"Genome size: {genome_size} | Is circular: {is_circular} | Seqid: {seqid} | Source: {source}"
        )
        ig_regions: list[IGRegion] = []

        # first line after region line (second line, first gene/overlap)
        gene_id = 1
        first_line = gff_file.readline().strip().split("\t")
        prev_end = 1
        if (c_start := int(first_line[3])) != new_region_start:
            c_att = first_line[8]
            log.debug(
                f"First row does not start from {new_region_start}. Start: {c_start}"
            )
            if (
                first_line[2][1] == "v"
            ):  # type = overlapping_feature_set || since the introduction of 'orf' type, we need to check for 'v' in the second character
                c_att = get_names(log, c_att, is_left=True)
                log.debug(f"Attributes get_name: {c_att}")
                prev_att = get_names(log, first_line[8], is_left=False)
            else:
                prev_att = c_att

            c_att = change_names(log, c_att, "dw")
            fatt = f"name_up=region_start;source_up=region_start;{c_att}"
            log.trace(
                f"New IG from {new_region_start} to {c_start - 1} | att: {fatt}",
            )
            ig_regions.append((gene_id, new_region_start, c_start - 1, fatt))
            gene_id += 1
            prev_end = int(first_line[4])
            rs = True
        else:
            log.debug(
                f"First row does start from {new_region_start}. No intergenic region needed."
            )
            prev_end = int(first_line[4])
            if first_line[2][1] == "v":  # type = overlapping_feature_set
                prev_att = get_names(log, first_line[8], is_left=False)
                first_att = get_names(log, first_line[8], is_left=True)
            else:
                prev_att = first_line[8]
                first_att = prev_att

        # main loop, starting at the second gene/overlap
        for loop_line in gff_file:
            line: list[str] = loop_line.strip().split("\t")
            log.trace(f"line: {line}")
            if line[0] == "":
                continue
            # check if there is a gap between the previous end and the current start
            if prev_end + 1 == (c_start := int(line[3])):
                log.trace(f"No gap between {prev_end} and {c_start}")
                prev_end = int(line[4])

                if line[2][1] == "v":  # type = overlapping_feature_set
                    prev_att = get_names(log, line[8], is_left=False)
                else:
                    prev_att = line[8]
                continue

            # if current gene is overlapping, get the attributes from the left gene
            if line[2][1] == "v":  # type = overlapping_feature_set
                att = change_names(log, prev_att, "up") + change_names(
                    log, get_names(log, line[8], is_left=True), "dw"
                )
                prev_att = get_names(log, line[8], is_left=False)
            else:
                att = change_names(log, prev_att, "up") + change_names(
                    log, line[8], "dw"
                )
                prev_att = line[8]
            log.trace(f"New IG from {prev_end + 1} to {c_start - 1} | att: {att}")
            ig_regions.append((gene_id, prev_end + 1, c_start - 1, att))
            gene_id += 1
            prev_end = int(line[4])

    # check if last gene/olverlap ends at the end of the genome
    if prev_end < genome_size:
        fatt = f"{change_names(log, prev_att, 'up')}name_dw=region_end;source_dw=region_end;"
        log.debug(f"Last row does not end at genome end {genome_size}. End: {prev_end}")
        log.trace(f"New IG from {prev_end + 1} to {genome_size} | att: {fatt}")
        ig_regions.append((gene_id, prev_end + 1, genome_size, fatt))
        re = True

    if is_circular:
        log.debug("Circular genome detected. Checking RS|RE.")
        if rs and re:  # merge the first and last regions
            log.debug("RS and RE are True, Merge required.")
            log.debug("Saving premerge gff3 file, to run bedtools on it.")
            premerge_gff = gff_out.parent / f"{gff_out.stem}_premerged.gff3"

            save_gff_with_frag(
                log, premerge_gff, ig_regions, seqid, source, header=header
            )
            if run_bedtools:
                premerge_fasta = (
                    fasta_out_path.parent / f"{fasta_out_path.stem}_premerged.fasta"
                )
                exec_bedtools(log, fasta_in_path, premerge_gff, premerge_fasta)

            # merge the first and last regions
            first_ig = ig_regions.pop(0)
            genome_size = ig_regions[-1][2]

            ig_regions[-1] = (
                ig_regions[-1][0],  # keep the gene_id of the last region
                ig_regions[-1][1],  # keep the start of the last region
                ig_regions[-1][2]
                + first_ig[2],  # end of the last region + end of the first
                get_names_ig(log, ig_regions[-1][3], suffix="up")
                + get_names_ig(log, first_ig[3], suffix="dw"),
            )  # merge attributes of the last and first regions

            log.debug("Saving merged final gff3 file.")
            save_gff_with_merge(log, gff_out, ig_regions, seqid, source, header=header)

            # update fasta file to reflect the merged region in gff
            if run_bedtools:
                log.debug("Updating fasta file to reflect the merged gff3.")
                from Bio import SeqIO
                from Bio.SeqRecord import SeqRecord

                copy2(premerge_fasta, fasta_out_path)
                log.trace(
                    f"copy2 | src: {premerge_fasta} | dst: {fasta_out_path}",
                )
                records: list[SeqRecord] = list(SeqIO.parse(fasta_out_path, "fasta"))  # type: ignore[no-untyped-call]
                drop = records.pop(0)
                new_seq = records.pop(-1)
                assert new_seq.seq is not None
                new_seq.seq += drop.seq
                new_seq.id = f"intergenic_region_merged::{seqid}:{ig_regions[-1][1] - 1}-{ig_regions[-1][2]}"
                new_seq.description = f" merged region {seqid}:{ig_regions[-1][1] - 1}-{genome_size} with {seqid}:{first_ig[1] - 1}-{first_ig[2]}"
                SeqIO.write(records + [new_seq], fasta_out_path, "fasta")

        elif re:  # update last region, no merge
            log.debug("Only RE is True, updating the last region.")

            ig_regions[-1] = (
                ig_regions[-1][0],  # keep the gene_id of the last region
                ig_regions[-1][1],  # keep the start of the last region
                ig_regions[-1][2],  # keep the end of the last region
                get_names_ig(log, ig_regions[-1][3], suffix="up")
                + change_names(log, first_att, "dw"),
            )
            save_gff(log, gff_out, ig_regions, seqid, source, header=header)

            if run_bedtools:
                exec_bedtools(log, fasta_in_path, gff_out, fasta_out_path)

        elif rs:  # update first region, no merge
            log.debug("Only RS is True. Updating the first region.")

            ig_regions[0] = (
                ig_regions[0][0],  # keep the gene_id of the first region
                ig_regions[0][1],  # keep the start of the first region
                ig_regions[0][2],  # keep the end of the first region
                change_names(log, prev_att, "up")
                + get_names_ig(log, ig_regions[0][3], suffix="dw"),
            )

            save_gff(log, gff_out, ig_regions, seqid, source, header=header)

            if run_bedtools:
                exec_bedtools(log, fasta_in_path, gff_out, fasta_out_path)

        else:  # no RS|RE
            log.debug("No RS|RE. Doing nothing about RS|RE.")
            save_gff(log, gff_out, ig_regions, seqid, source, header=header)
            if run_bedtools:
                exec_bedtools(log, fasta_in_path, gff_out, fasta_out_path)

    else:  # linear genome
        log.debug("Linear genome detected. Doing nothing about RS|RE.")
        save_gff(log, gff_out, ig_regions, seqid, source, header=header)
        if run_bedtools:
            exec_bedtools(log, fasta_in_path, gff_out, fasta_out_path)


def multiple_execution(
    log: log_setup.GDTLogger,
    tsv_path: Path,
    an_column: str = "AN",
    workers: int = 0,
    gff_ext: str = ".gff3",
    fasta_ext: str = ".fasta",
    out_suffix: str = "_intergenic",
    gff_suffix: str = "",
    fasta_suffix: str = "",
    have_fasta: bool = False,
) -> None:
    tsv = pd.read_csv(tsv_path, sep="\t")
    gff_in_builder = gff3_utils.PathBuilder(gff_ext).use_folder_builder(
        tsv_path.parent,
        gff_suffix,
    )
    gff_out_builder = gff3_utils.PathBuilder(gff_ext).use_folder_builder(
        tsv_path.parent,
        out_suffix,
    )
    gff3_utils.check_file_in_tsv(
        log,
        tsv,
        gff_in_builder,
        an_column,
    )

    if have_fasta:
        fasta_in_builder = gff3_utils.PathBuilder(fasta_ext).use_folder_builder(
            tsv_path.parent,
            fasta_suffix,
        )
        fasta_out_builder = gff3_utils.PathBuilder(fasta_ext).use_folder_builder(
            tsv_path.parent,
            out_suffix,
        )
        gff3_utils.check_file_in_tsv(
            log,
            tsv,
            fasta_in_builder,
            an_column,
        )

    log.info(f"Processing {len(tsv)} ANs with {workers} workers")
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        task = [
            executor.submit(
                execution,
                log,
                gff_in_builder.build(an),
                gff_out_builder.build(an),
                fasta_in_builder.build(an) if have_fasta else None,
                fasta_out_builder.build(an) if have_fasta else None,
                have_fasta,
            )
            for an in tsv[an_column]
        ]
    concurrent.futures.wait(task)
