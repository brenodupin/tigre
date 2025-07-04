# -*- coding: utf-8 -*-
"""CLI Entry Point for IGRE - InterGenic Region Extractor."""

from . import __version__


def cli_entrypoint() -> None:
    """Entry point for the CLI."""
    print("Hello from igre CLI!")
    print(f"IGRE version: {__version__}")
