# -*- coding: utf-8 -*-
"""tigre - Tool for InterGenic Region Extraction."""

__version__ = "0.1.0"

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from . import fasta_utils

from .clean import *
from .cli import *
from .gff3_utils import *
from .igr import *
from .log_setup import *

# Cache for optional submodules
_optional_modules: dict[str, Any] = {}


def __getattr__(name: str) -> Any:
    """Lazy import optional submodules."""
    if name == "fasta_utils":
        if "fasta_utils" not in _optional_modules:
            from . import fasta_utils

            _optional_modules["fasta_utils"] = fasta_utils
        return _optional_modules["fasta_utils"]

    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
