# -*- coding: utf-8 -*-
"""CLI Entry Point for IGRE - InterGenic Region Extractor."""


from . import __version__
from . import log_setup
from pathlib import Path

import argparse

C_RESET = "\033[0m"

def create_clean_gff3():
    pass


def cli_entrypoint() -> None:
    """Command line interface for the Gene Dictionary Tool (gdt)."""
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
        epilog=f"Source ~ \033[32mhttps://github.com/brenodupin/igre{C_RESET}",
    )
    main_parser.add_argument(
        "--version",
        action="version",
        version=f"igre {__version__}",
        help="Show the version of the igre package.",
    )

    subparsers = main_parser.add_subparsers(
        dest="command",
        required=True,
    )

    test_parser = subparsers.add_parser(
        "test",
        help="Test command for IGRE, to check if the CLI is working.",
        description="Just a simple test command to check if the CLI is working," \
        "if the logging is working, and if the version is correct.",
        parents=[global_flags],
    )
    test_parser.set_defaults(func=lambda args: print("Test command executed successfully!"))

    args = main_parser.parse_args()

    log = log_setup.setup_logger(args.debug, args.log, args.quiet)
    log.trace("CLI execution started")
    log.trace(f"call path: {Path().resolve()}")
    log.trace(f"cli  path: {Path(__file__)}")
    log.trace(f"args: {args}")

    if args.command == "test":
        log.debug("Executing test command")
        args.func(args)
