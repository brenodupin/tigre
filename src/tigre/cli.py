# -*- coding: utf-8 -*-
"""CLI Entry Point for tigre - Tool for InterGenic Region Extraction."""


import argparse
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, TypeAlias

from . import (
    __version__,
    clean,
    clean_gdt,
    gff3_utils,
    igr,
    log_setup,
)
from .fasta_utils import BIOPYTHON_AVAILABLE, biopython_wrapper

C_RESET = "\x1b[0m"
C_BLUE = "\x1b[36m"
C_GREEN = "\x1b[92m"
C_YELLOW = "\x1b[33m"
C_BOLD = "\x1b[1m"
C_BG_GRAY = "\x1b[100m"

C_BANNER = f""" T   I   G   {C_GREEN}E   {C_YELLOW}R{C_RESET}   \\    /\\
 {C_BLUE}|   |   |   |{C_RESET} x {C_BLUE}|{C_RESET}    )  ( ')
 {C_BOLD}{C_BG_GRAY}T   I   G   {C_YELLOW}R   {C_GREEN}E{C_RESET}   (  /  )
 o   n   e   e   x    \\(__)|
 o   t   n   g   t     jgs
 l   e   i   i   r
     r   c   o   a
 f           n   c
 o               t
 r               i
                 o
                 n"""

if TYPE_CHECKING:
    pass

MAX_CPU: int = os.cpu_count() or 1

_Subparser: TypeAlias = "argparse._SubParsersAction[argparse.ArgumentParser]"
_Parser: TypeAlias = "argparse.ArgumentParser"
_Group: TypeAlias = "argparse._ArgumentGroup"


def _workers_count(workers: int, threading: bool = False) -> int:
    """Return the number of workers to use."""
    cpus = MAX_CPU * 3 if threading else MAX_CPU
    return min(workers, cpus) if workers > 0 else cpus


def handle_single_result(
    log: log_setup.GDTLogger,
    result: tuple[bool, str, list[log_setup._RawMsg]],
    error_msg: str,
) -> None:
    """Handle the result of single execution."""
    success, an, records = result

    for record in records:
        log.log(record[0], record[1])

    if not success:
        log.error(f"An error occurred while processing {an}. {error_msg}")
        sys.exit(1)

    log.info("Processing completed successfully.")


def ensure_exists(
    log: log_setup.GDTLogger,
    path: Path,
    desc: str = "file",
) -> None:
    """Ensure that the specified path exists."""
    if not path.is_file():
        log.error(f"{desc} not found, expected at {path}")
        sys.exit(1)


def ensure_overwrite(
    log: log_setup.GDTLogger,
    path: Path,
    desc: str = "file",
    overwrite: bool = False,
) -> None:
    """Ensure that the specified path can be overwritten."""
    if path.is_file() and not overwrite:
        log.error(f"{desc} already exists at {path}. Use --overwrite to replace it.")
        sys.exit(1)


def args_multiple(
    parser: _Group,
    file: str = "gff",
    io: str = "in",
    ext: str = ".gff3",
    suffix: str = "_clean",
    req: bool = False,
) -> None:
    """Add common arguments for multiple file processing."""
    up = file.upper()
    parser.add_argument(
        f"--{file}-{io}-ext",
        required=req,
        type=str,
        default=ext,
        help=f"File Extension for {io}put {up} files. Default: '{ext}'.",
    )
    parser.add_argument(
        f"--{file}-{io}-suffix",
        required=req,
        type=str,
        default=suffix,
        help=f"Suffix to append when building {io}put {up} file paths from the TSV "
        f"file. Example: '{suffix}' creates paths like '<AN>{suffix}.gff3' if the "
        f"extension is '{ext}'. Default: '{suffix}'.",
    )


def args_single(
    parser: _Group,
    file: str = "gff",
    io: str = "in",
    required: bool = True,
) -> None:
    """Add common arguments for single file processing."""
    parser.add_argument(
        f"--{file}-{io}",
        required=required,
        type=str,
        help=f"{file.upper()} {io}put file.",
    )


def args_tsv(
    parser: argparse._ArgumentGroup,
    action: str,
) -> None:
    """Add arguments for TSV file processing."""
    parser.add_argument(
        "--tsv",
        required=True,
        type=str,
        help="TSV file with a column of accession numbers, from which file paths are "
        f"derived for {action}.",
    )
    parser.add_argument(
        "--an-column",
        required=False,
        type=str,
        default="AN",
        help="Column name containing the accession number in the TSV. Default: 'AN'.",
    )
    parser.add_argument(
        "--workers",
        required=False,
        type=int,
        default=0,
        help="Number of workers to use for parallel processing. Set to 0 to use all "
        f"available cores ({MAX_CPU}). Default: 0.",
    )


def args_log(parser: _Parser) -> None:
    """Add logging arguments to the parser."""
    group = parser.add_argument_group("log options")
    group.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity level. Use multiple times for more verbose output: "
        "-v (TRACE file), -vv (DEBUG console + TRACE file), "
        "-vvv (TRACE console + TRACE file). Default: INFO console + DEBUG file.",
    )
    group.add_argument(
        "--log",
        required=False,
        type=str,
        default=None,
        help="Path to the log file. If not provided, a default log file will be "
        "created.",
    )
    group.add_argument(
        "--no-log-file",
        required=False,
        action="store_true",
        default=False,
        help="Disable file logging. Default: False.",
    )

    group.add_argument(
        "--quiet",
        required=False,
        action="store_true",
        default=False,
        help="Suppress console output. Default: False.",
    )


###############################################
########## extract command functions ##########
def extract_group(parser: _Parser) -> None:
    """Add common arguments for extract command."""
    group = parser.add_argument_group("extract options")
    group.add_argument(
        "--add-region",
        required=False,
        action="store_true",
        default=False,
        help="Add region line to the output GFF3 file. Default: False.",
    )
    group.add_argument(
        "--overwrite",
        required=False,
        action="store_true",
        default=False,
        help="Overwrite existing output files. Default: False.",
    )
    group.add_argument(
        "--feature-type",
        required=False,
        type=str,
        default="intergenic_region",
        help="Feature type name to use for intervening regions in the output GFF3 "
        "file. Features spanning the genome boundary in circular genomes will be "
        "appended with '_merged'. Default: 'intergenic_region'.",
    )


def extract_parser(
    sub: _Subparser,
    global_args: _Parser,
) -> None:
    """Create extract command parser and add it to the subparsers."""
    extract = sub.add_parser(
        "extract",
        help="`extract` command for tigre, to extract intergenic regions "
        "from GFF3 file(s).",
        parents=[global_args],
    )
    extract_group(extract)

    extract_sub = extract.add_subparsers(metavar="mode", dest="mode", required=True)

    single_parser = extract_sub.add_parser(
        "single",
        help="Process a single GFF file",
        parents=[global_args],
    )
    extract_group(single_parser)
    single = single_parser.add_argument_group("single file options")
    args_single(single, "gff", "in")
    args_single(single, "gff", "out")

    multiple_parser = extract_sub.add_parser(
        "multiple",
        help="Process multiple GFF files from TSV",
        parents=[global_args],
    )
    extract_group(multiple_parser)
    multiple = multiple_parser.add_argument_group("multiple file options")
    args_tsv(multiple, "extract")
    args_multiple(multiple, "gff", "in", ".gff3", "_clean")
    args_multiple(multiple, "gff", "out", ".gff3", "_intergenic")


def extract_command(
    args: argparse.Namespace,
    log: log_setup.GDTLogger,
) -> None:
    """Execute the extract command based on the provided arguments."""
    if args.mode == "single":
        args.gff_in = Path(args.gff_in).resolve()
        args.gff_out = Path(args.gff_out).resolve()

        ensure_exists(
            log,
            args.gff_in,
            "GFF3 input",
        )
        ensure_overwrite(
            log,
            args.gff_out,
            "GFF3 output",
            args.overwrite,
        )

        result = igr.extract_intergenic_regions(
            log.spawn_buffer(),
            args.gff_in,
            args.gff_out,
            args.add_region,
            args.feature_type,
        )
        handle_single_result(log, result, "Error extracting intergenic regions single")

    elif args.mode == "multiple":
        args.tsv = Path(args.tsv).resolve()
        ensure_exists(
            log,
            args.tsv,
            "TSV file",
        )

        igr.extract_multiple(
            log,
            args.tsv,
            _workers_count(args.workers),
            args.an_column,
            args.gff_in_ext,
            args.gff_in_suffix,
            args.gff_out_ext,
            args.gff_out_suffix,
            args.add_region,
            args.overwrite,
            args.feature_type,
        )


###############################################


################################################
########## getfasta command functions ##########
def getfasta_group(parser: _Parser) -> None:
    """Add common arguments for getfasta command."""
    group = parser.add_argument_group("getfasta options")
    group.add_argument(
        "--bedtools-compatible",
        required=False,
        action="store_true",
        default=False,
        help="Adjust FASTA headers to use 0-based indexing for bedtools compatibility. "
        "Affects headers only, not sequences. Default: False.",
    )
    group.add_argument(
        "--overwrite",
        required=False,
        action="store_true",
        default=False,
        help="Overwrite existing output files. Default: False.",
    )


def getfasta_parser(
    sub: _Subparser,
    global_args: _Parser,
) -> None:
    """Create getfasta command parser and add it to the subparsers."""
    getfasta = sub.add_parser(
        "getfasta",
        help="`getfasta` command for tigre, to get fasta from intergenic regions.",
        parents=[global_args],
    )
    getfasta_group(getfasta)

    getfasta_sub = getfasta.add_subparsers(metavar="mode", dest="mode", required=True)

    single_parser = getfasta_sub.add_parser(
        "single",
        help="Process a single GFF file",
        parents=[global_args],
    )
    getfasta_group(single_parser)
    single = single_parser.add_argument_group("single file options")
    args_single(single, "gff", "in")
    args_single(single, "fasta", "in")
    args_single(single, "fasta", "out")

    multiple_parser = getfasta_sub.add_parser(
        "multiple",
        help="Process multiple GFF files from TSV",
        parents=[global_args],
    )
    getfasta_group(multiple_parser)
    multiple = multiple_parser.add_argument_group("multiple file options")
    args_tsv(multiple, "extract")
    args_multiple(multiple, "gff", "in", ".gff3", "_intergenic")
    args_multiple(multiple, "fasta", "in", ".fasta", "")
    args_multiple(multiple, "fasta", "out", ".fasta", "_intergenic")


def getfasta_command(
    args: argparse.Namespace,
    log: log_setup.GDTLogger,
) -> None:
    """Execute the Biopython getfasta command based on the provided arguments."""
    if not BIOPYTHON_AVAILABLE:
        log.error(
            "Biopython is not available. Please install it with: "
            "pip install biopython or pip install tigre[bio]"
        )
        sys.exit(1)

    if args.mode == "single":
        args.gff_in = Path(args.gff_in).resolve()
        args.fasta_in = Path(args.fasta_in).resolve()
        args.fasta_out = Path(args.fasta_out).resolve()
        ensure_exists(log, args.gff_in, "GFF3 input")
        ensure_exists(log, args.fasta_in, "Fasta input")
        ensure_overwrite(log, args.fasta_out, "Fasta output", args.overwrite)

        result = biopython_wrapper.biopython_getfasta(
            log.spawn_buffer(),
            args.gff_in,
            args.fasta_in,
            args.fasta_out,
            args.bedtools_compatible,
        )

        handle_single_result(log, result, "Error extracting sequences single")

    elif args.mode == "multiple":
        biopython_wrapper.biopython_multiple(
            log,
            Path(args.tsv).resolve(),
            _workers_count(args.workers),
            args.gff_in_ext,
            args.gff_in_suffix,
            args.fasta_in_ext,
            args.fasta_in_suffix,
            args.fasta_out_ext,
            args.fasta_out_suffix,
            args.an_column,
            args.bedtools_compatible,
            args.overwrite,
        )


################################################


#############################################
########## Clean command functions ##########
def clean_group(parser: _Parser) -> None:
    """Add common arguments for clean command."""
    group = parser.add_argument_group("clean options")
    group.add_argument(
        "--gdict",
        required=False,
        type=str,
        dest="gdt",
        default=None,
        help="Path to a Gene Dictionary Table (GDT) file to standardize gene names. "
        "Optional.",
    )
    group.add_argument(
        "--query-string",
        required=False,
        type=str,
        default=gff3_utils.QS_GENE_TRNA_RRNA_REGION,
        help=f"pandas query string to filter features in GFF3 file. Default: "
        f"'{gff3_utils.QS_GENE_TRNA_RRNA_REGION}'.",
    )
    group.add_argument(
        "--keep-orfs",
        required=False,
        action="store_true",
        default=False,
        help="Keep ORF sequences in the output. Default: False.",
    )
    group.add_argument(
        "--overwrite",
        required=False,
        action="store_true",
        default=False,
        help="Overwrite existing output files. Default: False.",
    )


def clean_parser(
    sub: _Subparser,
    global_args: _Parser,
) -> None:
    """Create clean command parser and add it to the subparsers."""
    clean = sub.add_parser(
        "clean",
        help="`clean` command for tigre, to clean the GFF3 files.",
        parents=[global_args],
    )
    clean_group(clean)

    clean_sub = clean.add_subparsers(metavar="mode", dest="mode", required=True)

    single_parser = clean_sub.add_parser(
        "single",
        help="Process a single GFF file",
        parents=[global_args],
    )
    clean_group(single_parser)
    single = single_parser.add_argument_group("single file options")
    args_single(single, "gff", "in")
    args_single(single, "gff", "out")

    multiple_parser = clean_sub.add_parser(
        "multiple",
        help="Process multiple GFF files from TSV",
        parents=[global_args],
    )
    clean_group(multiple_parser)
    multiple = multiple_parser.add_argument_group("multiple file options")
    args_tsv(multiple, "clean")
    args_multiple(multiple, "gff", "in", ".gff3", "")
    args_multiple(multiple, "gff", "out", ".gff3", "_clean")
    multiple.add_argument(
        "--server",
        required=False,
        action="store_true",
        default=False,
        help="Force the use of the GDT server to process GFF3 files, overriding "
        "auto-detection. Default: False.",
    )
    multiple.add_argument(
        "--no-server",
        required=False,
        action="store_true",
        default=False,
        help="Disable the use of the GDT server to process GFF3 files. Default: False.",
    )


def clean_command(
    args: argparse.Namespace,
    log: log_setup.GDTLogger,
) -> None:
    """Execute the clean command based on the provided arguments."""
    if args.mode == "single":
        args.gff_in = Path(args.gff_in).resolve()
        args.gff_out = Path(args.gff_out).resolve()

        ensure_exists(
            log,
            args.gff_in,
            "GFF3 input",
        )
        ensure_overwrite(
            log,
            args.gff_out,
            "GFF3 output",
            args.overwrite,
        )

        if args.gdt:
            result = clean_gdt.clean_an_gdt(
                log.spawn_buffer(),
                args.gff_in,
                args.gff_out,
                clean_gdt.load_gdt(log, args.gdict),
                args.query_string,
                args.keep_orfs,
            )

        else:
            result = clean.clean_an(
                log.spawn_buffer(),
                args.gff_in,
                args.gff_out,
                args.query_string,
                args.keep_orfs,
            )

        handle_single_result(log, result, "Error cleaning GFF3 single")

    elif args.mode == "multiple":
        args.tsv = Path(args.tsv).resolve()
        ensure_exists(
            log,
            args.tsv,
            "TSV file",
        )

        if args.gdt:
            # solve_gdt_call will decide between ProcessPoolExecutor or
            # ThreadPoolExecutor, depending on input shape and size.
            # This can be overruled by `--server` or `--no-server`
            clean_gdt.clean_gdt_multiple(
                log,
                args,
                _workers_count(args.workers),
                _workers_count(args.workers, threading=True),
            )

        else:
            clean.clean_multiple(
                log,
                args.tsv,
                _workers_count(args.workers),
                args.gff_in_ext,
                args.gff_in_suffix,
                args.gff_out_ext,
                args.gff_out_suffix,
                args.an_column,
                args.query_string,
                args.keep_orfs,
                args.overwrite,
            )


#############################################


def cli_entrypoint() -> int:
    """Command line interface for the tigre package."""
    # Global parser to add debug, log, and quiet flags to all subcommands
    start_time = datetime.now()

    global_args = argparse.ArgumentParser(add_help=False)
    args_log(global_args)

    main = argparse.ArgumentParser(
        description=C_BANNER,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[global_args],
        epilog=f"Source ~ \033[32mhttps://github.com/brenodupin/tigre{C_RESET}",
    )
    main.add_argument(
        "--version",
        action="version",
        version=f"tigre {__version__}",
        help="Show the version of the tigre package.",
    )

    subs = main.add_subparsers(dest="cmd", required=True)

    clean_parser(subs, global_args)
    extract_parser(subs, global_args)
    getfasta_parser(subs, global_args)

    subs.add_parser(
        "test",
        help="Test command to check CLI startup and logging.",
        parents=[global_args],
    )

    args = main.parse_args()

    if args.log and args.no_log_file:
        main.error(
            "You cannot specify both --log and --no-log-file together. "
            "Please choose one or neither."
        )

    log = log_setup.setup_logger(
        args.verbose,
        args.log,
        args.quiet,
        args.no_log_file,
        args.cmd,
    )
    log.info(f"TIGRE v{__version__} - Tool for InterGenic Region Extraction")
    log.info(
        f"Command: {args.cmd} | Mode: {args.mode if hasattr(args, 'mode') else 'N/A'}"
    )
    log.trace(f"Full args: {args}")
    log.debug(f"Working directory: {Path().resolve()}")
    log.trace(f"CLI script location: {Path(__file__).resolve()}")
    log.trace(f"Python version: {sys.version.split()[0]}")

    try:
        if args.cmd == "clean":
            if args.gdt and args.mode == "multiple" and args.server and args.no_server:
                main.error("You cannot specify both --server and --no-server together.")

            clean_command(args, log)

        elif args.cmd == "extract":
            extract_command(args, log)

        elif args.cmd == "getfasta":
            getfasta_command(args, log)

    except KeyboardInterrupt:
        log.warning("Process interrupted by user (Ctrl+C)")
        return 130  # SIGINT

    except SystemExit as e:
        log.debug(f"SystemExit raised with code: {e.code}")
        raise

    except Exception as e:
        log.error(f"Unexpected error occurred: {str(e)}")
        log.debug("Full traceback:", exc_info=True)

        return 1

    total_time = datetime.now() - start_time
    log.info(f"Execution completed in {_time_formatted(total_time.total_seconds())}")
    return 0


def _time_formatted(total_seconds: float) -> str:
    """Format total seconds into a human-readable string."""
    hours = int(total_seconds // 3600)
    minutes = int((total_seconds % 3600) // 60)
    seconds = total_seconds % 60

    parts = []
    if hours > 0:
        parts.append(f"{hours}h")
    if minutes > 0:
        parts.append(f"{minutes}m")
    if seconds > 0 or not parts:  # Always show seconds if nothing else
        parts.append(f"{seconds:.2f}s")

    return " ".join(parts)
