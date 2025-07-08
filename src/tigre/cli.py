# -*- coding: utf-8 -*-
"""CLI Entry Point for tigre - Tool for InterGenic Region Extraction."""


import argparse
import os
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Optional, Union

from . import __version__, clean, extraction, log_setup

C_RESET = "\033[0m"

if TYPE_CHECKING:
    import pandas as pd


def deploy_gdt_support(
    log: log_setup.GDTLogger,
    gdt_path: Union[str, Path],
) -> Callable[["pd.Series", log_setup.GDTLogger], str]:
    try:
        import gdt
    except ImportError as ex:
        raise SystemExit(
            "GDT package not found. Please install gdt package to use the --gdict option."
        )
    else:
        log.debug(f"GDT package version: {gdt.__version__}")
        gdt_path = Path(gdt_path).resolve()
        if not gdt_path.is_file():
            log.error(f"GDT .gdict file not found: {gdt_path}")
            raise FileNotFoundError(f"GDT .gdict file not found: {gdt_path}")

        gdict = gdt.read_gdict(gdt_path, lazy_info=False)
        gdt.log_info(log, gdict)

        def clean_att_gdt(
            row: "pd.Series",
            log: log_setup.GDTLogger,
        ) -> str:
            """Clean attributes using GDT."""
            try:
                label = gdict[row.gene_id].label
            except KeyError:
                log.error(
                    f"Gene ID {row.gene_id} not found in GDT dictionary. Using default format."
                )
                label = row.gene_id
            return f"name={label};source={row.seqid}|{row.type}|{row.start}|{row.end}|{row.strand}|{row.gene_id};"

        return clean_att_gdt


def add_extract_single_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--gff-in",
        "-gi",
        required=True,
        type=str,
        help="GFF3 input file",
    )
    parser.add_argument(
        "--gff-out",
        "-go",
        required=True,
        type=str,
        help="GFF3 output file",
    )
    parser.add_argument(
        "--fasta-in",
        "-fi",
        required=False,
        default=None,
        type=str,
        help="Fasta input file",
    )
    parser.add_argument(
        "--fasta-out",
        "-fo",
        required=False,
        default=None,
        type=str,
        help="Fasta output file",
    )
    parser.add_argument(
        "--new-region-start",
        "-nrs",
        required=False,
        type=int,
        default=1,
        help="Check the first row for region start (in case there's a tRNA that start before region end and ends after region start)",
    )


def extract_single_command(
    args: argparse.Namespace,
    log: log_setup.GDTLogger,
) -> None:
    args.gff_in = Path(args.gff_in).resolve()
    args.gff_out = Path(args.gff_out).resolve()

    run_bedtools = bool(args.fasta_in)
    log.debug(f"run_bedtools: {run_bedtools}")
    if run_bedtools:  # check for bedtools
        try:
            a = subprocess.run(
                ["bedtools", "--version"], check=True, capture_output=True
            )
            log.debug(f"bedtools --version: {a.stdout.decode().strip()}")
        except (FileNotFoundError, subprocess.CalledProcessError) as ex:
            log.error(f"Error checking bedtools version: {ex}")
            raise SystemExit(
                "bedtools not found. Please install bedtools if you want to create fasta from extracted intergenic regions."
            )

    args.fasta_in = Path(args.fasta_in).resolve() if args.fasta_in else None
    args.fasta_out = Path(args.fasta_out).resolve() if args.fasta_out else None

    extraction.execution(
        log,
        args.gff_in,
        args.gff_out,
        args.fasta_in,
        args.fasta_out,
        run_bedtools,
        args.new_region_start,
    )


def add_extract_multiple_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--tsv",
        required=True,
        type=str,
        help="TSV file with indexed GFF3 files to standardize.",
    )
    parser.add_argument(
        "--an-column",
        required=False,
        default="AN",
        type=str,
        help="Column name for NCBI Accession Number inside the TSV. Default: AN",
    )
    parser.add_argument(
        "--gff-ext",
        required=False,
        default=".gff3",
        type=str,
        help="File Extension for GFF files. Default: '.gff3'",
    )
    parser.add_argument(
        "--gff-suffix",
        required=False,
        default="_clean",
        type=str,
        help="Suffix to be added when building GFF Paths from the TSV file. "
        "Example: '_clean' will create GFF paths like '<AN>_clean.gff3' if "
        "--gff-ext is '.gff3'. Default: '_clean'",
    )
    parser.add_argument(
        "--have-fasta",
        "-hf",
        required=False,
        type=str,
        help="If true, will look for a fasta file using GFF3 of TSV. Expected format: '<AN><fasta_ext>'.",
    )
    parser.add_argument(
        "--fasta-ext",
        required=False,
        default=".fasta",
        type=str,
        help="File Extension for Fasta files. Default: '.fasta'",
    )
    parser.add_argument(
        "--fasta-suffix",
        required=False,
        default="",
        type=str,
        help="Suffix to be added when building Fasta Paths from the TSV file. "
        "Example: '_merged' will create Fasta paths like '<AN>_merged.fasta' if "
        "`--fasta-ext '.fasta'`. Default: ''",
    )
    parser.add_argument(
        "--out-suffix",
        required=False,
        default="_intergenic",
        type=str,
        help="Suffix to be added to the output GFF3 files (and Fasta if --have-fasta is true). Default: '_intergenic'",
    )
    parser.add_argument(
        "--workers",
        required=False,
        default=0,
        type=int,
        help="Number of workers to use. "
        f"Default: 0 (use all available cores: {os.cpu_count()})",
    )


def extract_multiple_command(
    args: argparse.Namespace,
    log: log_setup.GDTLogger,
) -> None:
    """Extract intergenic regions from multiple GFF3 files, indexed via a TSV file."""
    args.tsv = Path(args.tsv).resolve()
    if not args.tsv.is_file():
        log.error(f"TSV file not found: {args.tsv}")
        raise FileNotFoundError(
            f"TSV file not found: {args.tsv}. Please provide a valid TSV file."
        )

    extraction.multiple_execution(
        log,
        args.tsv,
        args.an_column,
        workers_count(args),
        args.gff_ext,
        args.fasta_ext,
        args.out_suffix,
        args.gff_suffix,
        args.fasta_suffix,
        args.have_fasta,
    )


def add_clean_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--tsv",
        required=True,
        type=str,
        help="TSV file with the accession numbers (ANs) to be cleaned.",
    )
    parser.add_argument(
        "--gff-ext",
        required=False,
        default=".gff3",
        type=str,
        help="File Extension for GFF files. Default: '.gff3'",
    )
    parser.add_argument(
        "--gff-suffix",
        required=False,
        default="",
        type=str,
        help="Suffix to be added when building GFF Paths from the TSV file. "
        "Example: '_clean' will create GFF paths like '<AN>_clean.gff3' if "
        "--gff-ext is '.gff3'. Default: ''",
    )
    parser.add_argument(
        "--out-suffix",
        required=False,
        default="_clean",
        type=str,
        help="Suffix to be added to the output GFF3 files. `--out-suffix` replaces "
        "`--gff-suffix` on the output. Example: `--out-suffix '_clean' --gff-suffix '_merged' --gff-ext '.gff3'`, "
        "will build the input like <AN>_merged.gff3 and the output like <AN>_clean.gff3."
        " Default: '_clean'",
    )
    parser.add_argument(
        "--an-column",
        required=False,
        default="AN",
        type=str,
        help="Column name for NCBI Accession Number inside the TSV. Default: AN",
    )
    parser.add_argument(
        "--workers",
        required=False,
        default=0,
        type=int,
        help="Number of workers to use. "
        f"Default: 0 (use all available cores: {os.cpu_count()})",
    )
    parser.add_argument(
        "--gdict",
        required=False,
        default=None,
        dest="gdt",
        type=str,
        help="Path to a GDT .gdict file, used in the cleaning of GFF3 file. This requires the `gdt` package to be installed.",
    )


def clean_command(args: argparse.Namespace, log: log_setup.GDTLogger) -> None:
    """Clean command for tigre, to clean the GFF3 files."""
    tsv_path = Path(args.tsv).resolve()
    if not tsv_path.is_file():
        log.error(f"TSV file not found: {tsv_path}")
        return

    clean_func = deploy_gdt_support(log, args.gdt) if args.gdt else clean.clean_attr

    clean.orchestration(
        log,
        tsv_path,
        args.gff_ext,
        args.gff_suffix,
        args.out_suffix,
        args.an_column,
        workers_count(args, threading=True),
        clean_func,
    )


def add_main_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--version",
        action="version",
        version=f"tigre {__version__}",
        help="Show the version of the tigre package.",
    )


def add_test(
    subparser: "argparse._SubParsersAction[argparse.ArgumentParser]",
    global_args: argparse.ArgumentParser,
) -> None:
    test = subparser.add_parser(
        "test",
        help="Test command for tigre, to check if the CLI is working.",
        description="Just a simple test command to check if the CLI is working,"
        "if the logging is working, and if the version is correct.",
        parents=[global_args],
    )
    test.set_defaults(func=lambda args: print("Test command executed successfully!"))


def add_global_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--debug",
        required=False,
        default=False,
        action="store_true",
        help="Enable TRACE level in file, and DEBUG on console. "
        "Default: DEBUG level on file and INFO on console.",
    )
    parser.add_argument(
        "--log",
        required=False,
        default=None,
        type=str,
        help="Path to the log file. "
        "If not provided, a default log file will be created.",
    )
    parser.add_argument(
        "--quiet",
        required=False,
        default=False,
        action="store_true",
        help="Suppress console output. Default: console output enabled.",
    )


def workers_count(args: argparse.Namespace, threading: bool = False) -> int:
    cpu_count = os.cpu_count() or 1

    max_limit = cpu_count * 4 if threading else cpu_count

    return max_limit if args.workers <= 0 else min(args.workers, max_limit)


def cli_entrypoint() -> None:
    """Command line interface for the tigre package."""
    # Global parser to add debug, log, and quiet flags to all subcommands
    global_args = argparse.ArgumentParser(add_help=False)
    add_global_args(global_args)

    main = argparse.ArgumentParser(
        description="<PLACE HOLDER BANNER>",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[global_args],
        epilog=f"Source ~ \033[32mhttps://github.com/brenodupin/tigre{C_RESET}",
    )
    add_main_args(main)

    subs = main.add_subparsers(dest="command", required=True)

    clean = subs.add_parser(
        "clean",
        help="Clean command for tigre, to clean the GFF3 files.",
        description="This command will clean the GFF3 files in the specified folder.",
        parents=[global_args],
    )
    add_clean_args(clean)

    extract_single = subs.add_parser(
        "extract-single",
        help="Extract intergenic regions from a (clean) GFF3 file.",
        description="This command will extract intergenic regions from the GFF3 file and optionally extract sequences from a FASTA file.",
        parents=[global_args],
    )
    add_extract_single_args(extract_single)

    extract_multiple = subs.add_parser(
        "extract-multiple",
        help="Extract intergenic regions from multiple (clean) GFF3 files, indexed via a TSV file.",
        description="This command will execute `extract-single` for each GFF3 file listed in the TSV file.",
        parents=[global_args],
    )
    add_extract_multiple_args(extract_multiple)

    add_test(subs, global_args)

    args = main.parse_args()

    log = log_setup.setup_logger(args.debug, args.log, args.quiet)
    log.trace("CLI execution started")
    log.trace(f"call path: {Path().resolve()}")
    log.trace(f"cli  path: {Path(__file__).resolve()}")
    log.trace(f"args: {args}")

    if args.command == "test":
        log.debug("Executing test command")
        args.func(args)

    elif args.command == "clean":
        log.debug("Executing clean command")
        clean_command(args, log)

    elif args.command == "extract-single":
        if bool(args.fasta_in) != bool(args.fasta_out):
            log.error(f"fin: {args.fasta_in} | fout: {args.fasta_out}")
            extract_single.error(
                "--fasta-in and --fasta-out must be provided together or not at all."
            )

        log.debug("Executing extract command")
        extract_single_command(args, log)

    elif args.command == "extract-multiple":
        extract_multiple_command(args, log)

    else:
        log.error(f"Unknown command: {args.command}")
        main.print_help()
        exit(1)
