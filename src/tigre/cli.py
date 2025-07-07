# -*- coding: utf-8 -*-
"""CLI Entry Point for tigre - Tool for InterGenic Region Extraction."""


import argparse
import subprocess
from pathlib import Path

from . import __version__, clean, extraction, log_setup

C_RESET = "\033[0m"


def extract_command(args: argparse.Namespace, log: log_setup.GDTLogger) -> None:
    args.gff_in = Path(args.gff_in)
    args.gff_out = Path(args.gff_out)

    run_bedtools = bool(args.fasta_in)
    log.debug(f"run_bedtools: {run_bedtools}")
    if run_bedtools:  # check for bedtools
        try:
            args.fasta_in = Path(args.fasta_in)
            args.fasta_out = Path(args.fasta_out)

            a = subprocess.run(
                ["bedtools", "--version"], check=True, capture_output=True
            )
            log.debug(f"bedtools --version: {a.stdout.decode().strip()}")
        except (FileNotFoundError, subprocess.CalledProcessError) as ex:
            log.error(f"Error checking bedtools version: {ex}")
            raise SystemExit(
                "bedtools not found. Please install bedtools if you want to create fasta from extracted intergenic regions."
            )

    extraction.execution(args, log)


def clean_command(args: argparse.Namespace, log: log_setup.GDTLogger) -> None:
    """Clean command for tigre, to clean the GFF3 files."""
    tsv_path = Path(args.tsv)
    if not tsv_path.is_file():
        log.error(f"TSV file not found: {tsv_path}")
        return

    clean.orchestration(tsv_path, log)


def cli_entrypoint() -> None:
    """Command line interface for the tigre package."""
    # Global parser to add debug, log, and quiet flags to all subcommands
    global_flags = argparse.ArgumentParser(add_help=False)
    global_flags.add_argument(
        "--debug",
        required=False,
        default=False,
        action="store_true",
        help="Enable TRACE level in file, and DEBUG on console. "
        "Default: DEBUG level on file and INFO on console.",
    )
    global_flags.add_argument(
        "--log",
        required=False,
        default=None,
        type=str,
        help="Path to the log file. "
        "If not provided, a default log file will be created.",
    )
    global_flags.add_argument(
        "--quiet",
        required=False,
        default=False,
        action="store_true",
        help="Suppress console output. Default: console output enabled.",
    )

    main_parser = argparse.ArgumentParser(
        description="<PLACE HOLDER BANNER>",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[global_flags],
        epilog=f"Source ~ \033[32mhttps://github.com/brenodupin/tigre{C_RESET}",
    )
    main_parser.add_argument(
        "--version",
        action="version",
        version=f"tigre {__version__}",
        help="Show the version of the tigre package.",
    )

    subparsers = main_parser.add_subparsers(
        dest="command",
        required=True,
    )

    test_parser = subparsers.add_parser(
        "test",
        help="Test command for tigre, to check if the CLI is working.",
        description="Just a simple test command to check if the CLI is working,"
        "if the logging is working, and if the version is correct.",
        parents=[global_flags],
    )
    test_parser.set_defaults(
        func=lambda args: print("Test command executed successfully!")
    )

    clean_parser = subparsers.add_parser(
        "clean",
        help="Clean command for tigre, to clean the GFF3 files.",
        description="This command will clean the GFF3 files in the specified folder.",
        parents=[global_flags],
    )
    clean_parser.add_argument(
        "--tsv",
        required=True,
        type=str,
        help="TSV file with the accession numbers (ANs) to be cleaned.",
    )

    ext_parser = subparsers.add_parser(
        "extract single",
        help="Extract intergenic regions from (clean) GFF3 files.",
        description="This command will extract intergenic regions from the GFF3 files and optionally extract sequences from a FASTA file.",
        parents=[global_flags],
    )
    ext_parser.add_argument(
        "--gff-in",
        "-gi",
        required=True,
        type=str,
        help="GFF3 input file",
    )
    ext_parser.add_argument(
        "--gff-out",
        "-go",
        required=True,
        type=str,
        help="GFF3 output file",
    )
    ext_parser.add_argument(
        "--fasta-in",
        "-fi",
        required=False,
        type=str,
        help="Fasta input file",
    )
    ext_parser.add_argument(
        "--fasta-out",
        "-fo",
        required=False,
        type=str,
        help="Fasta output file",
    )
    ext_parser.add_argument(
        "--new-region-start",
        "-nrs",
        required=False,
        type=int,
        default=1,
        help="Check the first row for region start (in case there's a tRNA that start before region end and ends after region start)",
    )

    args = main_parser.parse_args()

    log = log_setup.setup_logger(args.debug, args.log, args.quiet)
    log.trace("CLI execution started")
    log.trace(f"call path: {Path().resolve()}")
    log.trace(f"cli  path: {Path(__file__)}")
    log.trace(f"args: {args}")

    if args.command == "test":
        log.debug("Executing test command")
        args.func(args)
    elif args.command == "clean":
        log.debug("Executing clean command")
        clean_command(args, log)
    elif args.command == "extract":
        if bool(args.fasta_in) != bool(args.fasta_out):
            log.error(f"fin: {args.fasta_in} | fout: {args.fasta_out}")
            ext_parser.error(
                "--fasta-in and --fasta-out must be provided together or not at all."
            )
        log.debug("Executing extract command")
        extract_command(args, log)
    else:
        log.error(f"Unknown command: {args.command}")
        main_parser.print_help()
        exit(1)
