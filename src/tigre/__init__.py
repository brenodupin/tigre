# -*- coding: utf-8 -*-
"""tigre - Tool for InterGenic Region Extraction."""

import importlib.metadata

try:
    __version__ = importlib.metadata.version(__name__)
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0-dev"  # Fallback for development mode

from .clean import *
from .clean_gdt import *
from .cli import *
from .combine import *
from .getfasta import *
from .gff3_utils import *
from .igr import *
from .log_setup import *
