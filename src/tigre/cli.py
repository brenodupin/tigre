# -*- coding: utf-8 -*-
"""CLI Entry Point for tigre - Tool for InterGenic Region Extraction."""


import argparse
import os
import sys
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Union

from . import __version__, clean, gff3_utils, igr, log_setup
from .fasta_utils import BIOPYTHON_AVAILABLE, bedtools_wrapper, biopython_wrapper

C_RESET = "\033[0m"

if TYPE_CHECKING:
    import gdt  # type: ignore[import-not-found]
    import pandas as pd


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
        log.error(f"{error_msg}: {an}")
        sys.exit(1)


def ensure_exists(
    log: log_setup.GDTLogger,
    path: Path,
    desc: str = "file",
) -> None:
    if not path.is_file():
        log.error(f"{desc.upper()} not found: {path}")
        sys.exit(1)


def ensure_not_exists_or_overwrite(
    log: log_setup.GDTLogger,
    path: Path,
    desc: str = "file",
    overwrite: bool = False,
) -> None:
    if path.is_file() and not overwrite:
        log.error(
            f"{desc.upper()} already exists: {path}. Use --overwrite to replace it."
        )
        sys.exit(1)


def args_multiple(
    parser: argparse._ArgumentGroup,
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
        help=f"File Extension for {io}put {up} files. Default: '{ext}'",
    )
    parser.add_argument(
        f"--{file}-{io}-suffix",
        required=req,
        type=str,
        default=suffix,
        help=f"Suffix to be added when building {io}put {up} Paths from the TSV file. "
        f"Example: '{suffix}' will create {up} paths like '<AN>{suffix}.gff3' if "
        f"--{file}-{io}-ext is '{ext}'. Default: '{suffix}'",
    )


def args_single(
    parser: argparse._ArgumentGroup,
    file: str = "gff",
    io: str = "in",
    req: bool = True,
) -> None:
    """Add common arguments for single file processing."""
    parser.add_argument(
        f"--{file}-{io}",
        required=req,
        type=str,
        help=f"{file.upper()} {io}put file",
    )


def args_tsv(
    parser: argparse._ArgumentGroup,
    action: str,
) -> None:

    parser.add_argument(
        "--tsv",
        required=True,
        type=str,
        help="TSV file with column of Accession Numbers,"
        f" which file paths are derive from to {action}",
    )
    parser.add_argument(
        "--an-column",
        required=False,
        default="AN",
        type=str,
        help="Column name for Accession Number inside the TSV. Default: 'AN'",
    )
    parser.add_argument(
        "--workers",
        required=False,
        default=0,
        type=int,
        help="Number of workers to use in parallel processing. "
        f"Default: 0 (use all available cores: {os.cpu_count()})",
    )


# pre-compile functions for GDT cleaning
def _gdt_clean(
    gdt: "gdt.GeneDictionary",
    row: "pd.Series",
    log: log_setup.TempLogger,
) -> str:
    """Clean attributes using GDT."""
    try:
        label = gdt[row.gene_id].label
    except KeyError:
        log.error(
            f"Gene ID {row.gene_id} not found in GDT dictionary. Using default format."
        )
        label = row.gene_id
    return f"name={label};source={row.seqid}|{row.type}|{row.start}|{row.end}|{row.strand}|{row.gene_id};"


def _gdt(
    log: log_setup.GDTLogger,
    gdt_path: Union[str, Path],
) -> Callable[["pd.Series", log_setup.TempLogger], str]:
    try:
        import gdt
    except ImportError:
        raise SystemExit(
            "GDT package not found. Please install gdt package to use the --gdict option."
        )

    log.debug(f"GDT package version: {gdt.__version__}")
    gdt_path = Path(gdt_path).resolve()
    if not gdt_path.is_file():
        log.error(f"GDT .gdict file not found: {gdt_path}")
        sys.exit(1)

    gdict = gdt.read_gdict(gdt_path, lazy_info=False)
    log.info(f"GDT dictionary loaded from {gdt_path}")
    gdt.log_info(log, gdict)
    gdt_clean = partial(_gdt_clean, gdict)
    log.debug("GDT support deployed successfully.")
    return gdt_clean


def extract_single_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--gff-in",
        required=True,
        type=str,
        help="GFF3 input file",
    )
    parser.add_argument(
        "--gff-out",
        required=True,
        type=str,
        help="GFF3 output file",
    )
    parser.add_argument(
        "--add-region",
        required=False,
        default=False,
        action="store_true",
        help="If true, adds the region line to the output GFF3 file. Default: False",
    )


def extract_single_command(
    args: argparse.Namespace,
    log: log_setup.GDTLogger,
) -> None:
    args.gff_in = Path(args.gff_in).resolve()
    args.gff_out = Path(args.gff_out).resolve()

    if not args.gff_in.is_file():
        log.error(f"GFF3 input file not found: {args.gff_in}")
        return

    success, an, records = igr.extract_intergenic_regions(
        log.spawn_buffer(),
        args.gff_in,
        args.gff_out,
        args.add_region,
    )

    for record in records:
        log.log(record[0], record[1])

    if not success:
        log.error(f"Error extracting intergenic regions single: {an}")
        return


def extract_multiple_args(parser: argparse.ArgumentParser) -> None:
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
        "--gff-in-ext",
        required=False,
        default=".gff3",
        type=str,
        help="File Extension for Input GFF files. Default: '.gff3'",
    )
    parser.add_argument(
        "--gff-in-suffix",
        required=False,
        default="_clean",
        type=str,
        help="Suffix to be added when building Input GFF Paths from the TSV file. "
        "Example: '_clean' will create GFF paths like '<AN>_clean.gff3' if "
        "--gff-ext is '.gff3'. Default: '_clean'",
    )
    parser.add_argument(
        "--gff-out-ext",
        required=False,
        default=".gff3",
        type=str,
        help="File Extension for Output GFF files. Default: '.gff3'",
    )
    parser.add_argument(
        "--gff-out-suffix",
        required=False,
        default="_intergenic",
        type=str,
        help="Suffix to be added when building Output GFF Paths from the TSV file. "
        "Example: '_intergenic' will create GFF paths like '<AN>_intergenic.gff3' if "
        "--gff-ext is '.gff3'. Default: '_intergenic'",
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
        "--add-region",
        required=False,
        default=False,
        action="store_true",
        help="If true, adds the region line to the output GFF3 file. Default: False",
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

    igr.extract_multiple(
        log,
        args.tsv,
        args.an_column,
        _workers_count(args.workers, threading=False),
        args.gff_in_ext,
        args.gff_in_suffix,
        args.gff_out_ext,
        args.gff_out_suffix,
        args.add_region,
    )


def clean_command(
    args: argparse.Namespace,
    log: log_setup.GDTLogger,
) -> None:
    """Clean command for tigre, to clean the GFF3 files."""

    clean_func = _gdt(log, args.gdt) if args.gdt else clean.clean_attr

    if args.mode == "single":
        args.gff_in = Path(args.gff_in).resolve()
        args.gff_out = Path(args.gff_out).resolve()

        ensure_exists(
            log,
            args.gff_in,
            "GFF3 input",
        )
        ensure_not_exists_or_overwrite(
            log,
            args.gff_out,
            "GFF3 output",
            args.overwrite,
        )

        result = clean.clean_an(
            log.spawn_buffer(),
            args.gff_in,
            args.gff_out,
            clean_func,
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

        clean.clean_multiple(
            log,
            args.tsv,
            args.gff_in_ext,
            args.gff_in_suffix,
            args.gff_out_ext,
            args.gff_out_suffix,
            args.an_column,
            _workers_count(args.workers, threading=False),
            clean_func,
            args.query_string,
            args.keep_orfs,
            args.overwrite,
        )


def getfasta_multiple_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--tsv",
        required=True,
        type=str,
        help="TSV file with indexed GFF3 files to standardize.",
    )
    extraction_group = parser.add_mutually_exclusive_group()
    extraction_group.add_argument(
        "--bedtools",
        required=False,
        default=False,
        action="store_true",
        help="If true, will use bedtools to extract sequences from the GFF3 files.",
    )
    extraction_group.add_argument(
        "--biopython",
        required=False,
        default=False,
        action="store_true",
        help="If true, will use Biopython to extract sequences from the GFF3 files.",
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
        default="_intergenic",
        type=str,
        help="Suffix to be added when building GFF Paths from the TSV file. "
        "Example: '_intergenic' will create GFF paths like '<AN>_intergenic.gff3' if "
        "--gff-ext is '.gff3'. Default: '_intergenic'",
    )
    parser.add_argument(
        "--fasta-in-ext",
        required=False,
        default=".fasta",
        type=str,
        help="File Extension for Input Fasta files. Default: '.fasta'",
    )
    parser.add_argument(
        "--fasta-in-suffix",
        required=False,
        default="",
        type=str,
        help="Suffix to be added when building Input Fasta Paths from the TSV file. "
        "Example: '_merged' will create Fasta paths like '<AN>_merged.fasta' if "
        "`--fasta-in-ext '.fasta'`. Default: ''",
    )
    parser.add_argument(
        "--fasta-out-ext",
        required=False,
        default=".fasta",
        type=str,
        help="File Extension for Output Fasta files. Default: '.fasta'",
    )
    parser.add_argument(
        "--fasta-out-suffix",
        required=False,
        default="_intergenic",
        type=str,
        help="Suffix to be added when building Output Fasta Paths from the TSV file. "
        "Example: '_merged' will create Fasta paths like '<AN>_merged.fasta' if "
        "`--fasta-in-ext '.fasta'`. Default: '_intergenic'",
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
        "--bedtools-path",
        required=False,
        default="bedtools",
        type=str,
        help="Path to the bedtools executable. Default: 'bedtools'. "
        "If bedtools is installed in your PATH, you can leave this as default. "
        "Use with --bedtools.",
    )
    parser.add_argument(
        "--name-args",
        required=False,
        default="-name+",
        type=str,
        help="Arguments to pass to bedtools getfasta. Default: '-name+' "
        "(uses feature name as sequence ID in output FASTA). "
        "Use with --bedtools.",
    )
    parser.add_argument(
        "--bedtools-compatible",
        required=False,
        default=False,
        action="store_true",
        help="If true, adjust FASTA sequence labels to use 0-based indexing (bedtools compatible). "
        "By default, labels are adjusted to match GFF3 1-based indexing. "
        "Only affects sequence labels, not the actual sequences. Use with --biopython.",
    )


def getfasta_multiple_command(
    args: argparse.Namespace, log: log_setup.GDTLogger
) -> None:
    tsv_path = Path(args.tsv).resolve()

    if not tsv_path.is_file():
        log.error(f"TSV file not found: {tsv_path}")
        return

    if bool(args.bedtools) == bool(args.biopython):
        log.error("You must specify only one of --bedtools or --biopython.")
        return

    if args.bedtools:
        bedtools_version = bedtools_wrapper.get_bedtools_version(args.bedtools_path)
        log.debug(f"bedtools version: {bedtools_version}")

        if bedtools_version == None:
            log.error("bedtools not found or not executable. Please install bedtools.")
            return

        bedtools_wrapper.bedtools_multiple(
            log,
            tsv_path,
            args.gff_ext,
            args.gff_suffix,
            args.fasta_in_ext,
            args.fasta_in_suffix,
            args.fasta_out_ext,
            args.fasta_out_suffix,
            args.an_column,
            _workers_count(args.workers, threading=False),
            args.bedtools_path,
            args.name_args,
        )

    elif args.biopython:
        if not BIOPYTHON_AVAILABLE:
            log.error(
                "Biopython is not available. Please install it with: pip install biopython or pip install tigre[bio]"
            )
            return

        biopython_wrapper.biopython_multiple(
            log,
            tsv_path,
            args.gff_ext,
            args.gff_suffix,
            args.fasta_in_ext,
            args.fasta_in_suffix,
            args.fasta_out_ext,
            args.fasta_out_suffix,
            args.an_column,
            _workers_count(args.workers, threading=False),
            args.bedtools_compatible,
        )


def getfasta_single_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--gff-in",
        required=True,
        type=str,
        help="GFF3 input file",
    )
    parser.add_argument(
        "--fasta-in",
        required=True,
        type=str,
        help="Fasta input file",
    )
    parser.add_argument(
        "--fasta-out",
        required=True,
        type=str,
        help="Fasta output file",
    )
    extraction_group = parser.add_mutually_exclusive_group()
    extraction_group.add_argument(
        "--bedtools",
        required=False,
        default=False,
        action="store_true",
        help="If true, will use bedtools to extract sequences from the GFF3 file.",
    )
    extraction_group.add_argument(
        "--biopython",
        required=False,
        default=False,
        action="store_true",
        help="If true, will use Biopython to extract sequences from the GFF3 file.",
    )
    parser.add_argument(
        "--bedtools-path",
        required=False,
        default="bedtools",
        type=str,
        help="Path to the bedtools executable. Default: 'bedtools'. "
        "If bedtools is installed in your PATH, you can leave this as default. "
        "Use with --bedtools.",
    )
    parser.add_argument(
        "--name-args",
        required=False,
        default="-name+",
        type=str,
        help="Arguments to pass to bedtools getfasta. Default: '-name+' "
        "(uses feature name as sequence ID in output FASTA). "
        "Use with --bedtools.",
    )
    parser.add_argument(
        "--bedtools-compatible",
        required=False,
        default=False,
        action="store_true",
        help="If true, adjust FASTA sequence labels to use 0-based indexing (bedtools compatible). "
        "By default, labels are adjusted to match GFF3 1-based indexing. "
        "Only affects sequence labels, not the actual sequences. Use with --biopython.",
    )


def getfasta_single_command(
    args: argparse.Namespace,
    log: log_setup.GDTLogger,
) -> None:
    args.gff_in = Path(args.gff_in).resolve()
    args.fasta_in = Path(args.fasta_in).resolve()
    args.fasta_out = Path(args.fasta_out).resolve()

    if not args.gff_in.is_file():
        log.error(f"GFF3 input file not found: {args.gff_in}")
        return

    if not args.fasta_in.is_file():
        log.error(f"Fasta input file not found: {args.fasta_in}")
        return

    if bool(args.bedtools) == bool(args.biopython):
        log.error("You must specify only one of --bedtools or --biopython.")
        return

    if args.bedtools:
        bedtools_version = bedtools_wrapper.get_bedtools_version(args.bedtools_path)
        log.debug(f"bedtools version: {bedtools_version}")

        if bedtools_version == None:
            log.error("bedtools not found or not executable. Please install bedtools.")
            return

        success, an, records = bedtools_wrapper.bedtools_getfasta(
            log.spawn_buffer(),
            args.gff_in,
            args.fasta_in,
            args.fasta_out,
            args.bedtools_path,
            args.name_args,
        )

    elif args.biopython:
        if not BIOPYTHON_AVAILABLE:
            log.error(
                "Biopython is not available. Please install it with: pip install biopython or pip install tigre[bio]"
            )
            return

        success, an, records = biopython_wrapper.biopython_getfasta(
            log.spawn_buffer(),
            args.gff_in,
            args.fasta_in,
            args.fasta_out,
            args.bedtools_compatible,
        )

    for record in records:
        log.log(record[0], record[1])

    if not success:
        log.error(f"Error extracting sequences with bedtools: {an}")
        return


def main_args(parser: argparse.ArgumentParser) -> None:
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


def args_log(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("log options")
    group.add_argument(
        "--debug",
        required=False,
        default=False,
        action="store_true",
        help="Enable TRACE level in file, and DEBUG on console. "
        "Default: DEBUG level on file and INFO on console.",
    )
    group.add_argument(
        "--log",
        required=False,
        default=None,
        type=str,
        help="Path to the log file. "
        "If not provided, a default log file will be created.",
    )
    group.add_argument(
        "--quiet",
        required=False,
        default=False,
        action="store_true",
        help="Suppress console output. Default: console output enabled.",
    )


def add_clean_group(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("clean options")
    group.add_argument(
        "--gdict",
        required=False,
        type=str,
        dest="gdt",
        help="Gene dictionary file",
    )
    group.add_argument(
        "--query-string",
        required=False,
        type=str,
        help="Query string for filtering",
    )
    group.add_argument(
        "--keep-orfs",
        required=False,
        action="store_true",
        help="Keep ORF sequences in output",
    )
    group.add_argument(
        "--overwrite",
        required=False,
        action="store_true",
        help="Overwrite existing output files. Default: False (do not overwrite).",
    )


def clean_parser(
    sub: "argparse._SubParsersAction[argparse.ArgumentParser]",
    global_args: argparse.ArgumentParser,
) -> None:
    """Create clean command parser and add it to the subparsers."""
    clean = sub.add_parser(
        "clean",
        help="Clean command for tigre, to clean the GFF3 files.",
        parents=[global_args],
    )
    add_clean_group(clean)

    clean_sub = clean.add_subparsers(metavar="mode", dest="mode", required=True)

    single_parser = clean_sub.add_parser(
        "single",
        help="Process a single GFF file",
        parents=[global_args],
    )
    add_clean_group(single_parser)
    single = single_parser.add_argument_group("single file options")
    args_single(single, "gff", "in")
    args_single(single, "gff", "out")

    multiple_parser = clean_sub.add_parser(
        "multiple",
        help="Process multiple GFF files from TSV",
        parents=[global_args],
    )
    add_clean_group(multiple_parser)
    multiple = multiple_parser.add_argument_group("multiple file options")
    args_tsv(multiple, "clean")
    args_multiple(multiple, "gff", "in", ".gff3", "")
    args_multiple(multiple, "gff", "out", ".gff3", "_clean")


def _workers_count(
    user_workers: int,
    *,
    threading: bool = False,
    process_multiplier: int = 1,
    thread_multiplier: int = 4,
) -> int:
    cpu_count = os.cpu_count() or 1

    max_limit = (
        cpu_count * thread_multiplier if threading else cpu_count * process_multiplier
    )

    return max_limit if user_workers <= 0 else min(user_workers, max_limit)


def cli_entrypoint() -> None:
    """Command line interface for the tigre package."""
    # Global parser to add debug, log, and quiet flags to all subcommands
    global_args = argparse.ArgumentParser(add_help=False)
    args_log(global_args)

    main = argparse.ArgumentParser(
        description="<PLACE HOLDER BANNER>",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[global_args],
        epilog=f"Source ~ \033[32mhttps://github.com/brenodupin/tigre{C_RESET}",
    )
    main_args(main)

    subs = main.add_subparsers(dest="command", required=True)

    clean_parser(subs, global_args)

    extract_single = subs.add_parser(
        "extract-single",
        help="Extract intergenic regions from a (clean) GFF3 file.",
        description="This command will extract intergenic regions from the GFF3 file and optionally extract sequences from a FASTA file.",
        parents=[global_args],
    )
    extract_single_args(extract_single)

    extract_multiple = subs.add_parser(
        "extract-multiple",
        help="Extract intergenic regions from multiple (clean) GFF3 files, indexed via a TSV file.",
        description="This command will execute `extract-single` for each GFF3 file listed in the TSV file.",
        parents=[global_args],
    )
    extract_multiple_args(extract_multiple)

    getfasta_multiple = subs.add_parser(
        "getfasta-multiple",
        help="Extract sequences from multiple GFF3 files, indexed via a TSV file.",
        description="This command will extract sequences from the GFF3 files listed in the TSV file, using either bedtools or Biopython.",
        parents=[global_args],
    )
    getfasta_multiple_args(getfasta_multiple)

    getfasta_single = subs.add_parser(
        "getfasta-single",
        help="Extract sequences from a GFF3 file.",
        description="This command will extract sequences from a GFF3 file, using either bedtools or Biopython.",
        parents=[global_args],
    )
    getfasta_single_args(getfasta_single)

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
        log.debug("Executing extract command")
        extract_single_command(args, log)

    elif args.command == "extract-multiple":
        log.debug("Executing extract multiple command")
        extract_multiple_command(args, log)

    elif args.command == "getfasta-multiple":
        log.debug("Executing getfasta multiple command")
        getfasta_multiple_command(args, log)

    elif args.command == "getfasta-single":
        log.debug("Executing getfasta single command")
        getfasta_single_command(args, log)

    else:
        log.error(f"Unknown command: {args.command}")
        main.print_help()
        exit(1)
