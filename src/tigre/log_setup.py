# -*- coding: utf-8 -*-
"""Logger setup and management for tigre.

The base code for this module (log_setup.py from GDT package) has been
adapted for the TIGRE project.

Source
https://github.com/brenodupin/gdt/blob/master/src/gdt/log_setup.py
"""

import datetime
import glob
import logging
import os
from pathlib import Path
from typing import Any, Optional, Union, cast

TRACE = 5

_RawMsg = tuple[int, str]


class TempLogger:
    """Temporary logger to buffer log messages.

    This logger buffers messages in memory and can be used to collect logs
    during a specific operation, such as cleaning GFF files. It supports
    logging at various levels, including a custom TRACE level.

    This approach is thread-safe and allows for easy collection of log records
    without writing to a file or console immediately.
    """

    __slots__ = ("buffer", "trace_enable")

    def __init__(self, trace_enable: bool = False) -> None:
        """Initialize the buffered logger."""
        self.buffer: list[_RawMsg] = []
        self.trace_enable = trace_enable

    def _log(self, level: int, msg: str) -> None:
        """Log a message at the specified level."""
        self.buffer.append((level, msg))

    def error(self, msg: str) -> None:
        """Log an error message."""
        self._log(logging.ERROR, msg)

    def warning(self, msg: str) -> None:
        """Log a warning message."""
        self._log(logging.WARNING, msg)

    def info(self, msg: str) -> None:
        """Log an info message."""
        self._log(logging.INFO, msg)

    def debug(self, msg: str) -> None:
        """Log a debug message."""
        self._log(logging.DEBUG, msg)

    def trace(self, msg: str) -> None:
        """Log a trace message."""
        if self.trace_enable:
            self._log(TRACE, msg)

    def get_records(self) -> list[_RawMsg]:
        """Return the list of log records."""
        return self.buffer.copy()

    def clear(self) -> None:
        """Clear the buffer."""
        self.buffer.clear()


class GDTLogger(logging.Logger):
    """Extended logger class for GDT with TRACE (5) level support."""

    def __init__(self, name: str, level: int = logging.NOTSET) -> None:
        """Initialize the GDT logger with custom TRACE level."""
        super().__init__(name, level)
        self._trace_enabled = self.level <= TRACE

    def setLevel(self, level: Union[int, str]) -> None:  # noqa: N802
        """Override setLevel to update cached trace status."""
        super().setLevel(level)
        self._trace_enabled = self.level <= TRACE

    def trace(self, message: Any, *args: Any, **kwargs: Any) -> None:
        """Log 'msg % args' with severity 'TRACE'.

        Trace is a custom level below DEBUG (10), but above NOTSET (0), valued at 5.

        To pass exception information, use the keyword argument exc_info with
        a true value, e.g.

        logger.trace("Houston, we have a %s", "thorny problem", exc_info=1)
        """
        if self.isEnabledFor(TRACE):
            self._log(TRACE, message, args, **kwargs)

    def spawn_buffer(self) -> TempLogger:
        """Create a buffered logger to collect log records."""
        return TempLogger(self._trace_enabled)


logging.addLevelName(TRACE, "TRACE")
logging.setLoggerClass(GDTLogger)

_logging_levels: dict[str, int] = {
    "DISABLE": logging.CRITICAL + 1,  # above CRITICAL, used to disable logging
    "CRITICAL": logging.CRITICAL,
    "ERROR": logging.ERROR,
    "WARNING": logging.WARNING,
    "INFO": logging.INFO,
    "DEBUG": logging.DEBUG,
    "TRACE": TRACE,
    "NOTSET": logging.NOTSET,
}


def _cleanup_logs(log_dir: Path, max_files: int = 10) -> None:
    """Remove old log files, keeping only the most recent ones.

    Args:
        log_dir (Path): Directory where log files are stored.
        max_files (int): Maximum number of log files to keep. Defaults to 10.

    """
    log_files = sorted(glob.glob(str(log_dir / "tigre_*.log")))
    for old_file in log_files[: -(max_files - 1)]:
        try:
            os.remove(old_file)
        except Exception as e:
            print(f"Error removing old log file {old_file}: {e}")
            raise


def create_dev_logger(
    console_level: str = "INFO",
    file_level: str = "DEBUG",
    log_file: Optional[Path] = None,
) -> GDTLogger:
    """Set up the logger for the GDT package.

    Args:
        console_level (str): Logging level for console output.
        file_level (str): Logging level for file output.
        log_file (Optional[Path]): Path to the log file. If None, a default
                                   log file will be created at the project root,

    Returns:
        GDTLogger: Configured logger instance.

    """
    console_level_int = _logging_levels.get(console_level, logging.INFO)
    file_level_int = _logging_levels.get(file_level, logging.DEBUG)

    if log_file:
        log_file_path = Path(log_file).resolve()
        log_file_path.touch(exist_ok=True)

    else:
        project_root = Path(__file__).parent.parent.parent
        log_dir = project_root / "logs"

        log_dir.mkdir(exist_ok=True)

        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H_%M_%S")
        log_file_path = log_dir / f"tigre_{timestamp}.log"
        _cleanup_logs(log_dir)

    # Create and configure logger
    log = cast(GDTLogger, logging.getLogger("gdt"))
    log.setLevel(min(console_level_int, file_level_int))

    # Remove any existing handlers (in case logger was already configured)
    for handler in log.handlers[:]:
        log.removeHandler(handler)

    # Create console handler
    # (StreamHandler defaults to sys.stderr, can be changed to sys.stdout)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(console_level_int)
    console_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console_handler.setFormatter(console_formatter)
    log.addHandler(console_handler)

    # Create file handler
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(file_level_int)
    file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(file_formatter)
    log.addHandler(file_handler)

    log.propagate = False

    log.debug("Dev log setup complete")
    log.debug(f"Console logging level {console_level}")
    log.debug(f"File logging level {file_level} at {log_file_path}")
    return log


def create_simple_logger(
    print_to_console: bool = True,
    console_level: str = "INFO",
    save_to_file: bool = True,
    file_level: str = "DEBUG",
    log_file: Union[Path, str, None] = None,
) -> GDTLogger:
    """Create a simple logger with optional console and file output.

    Args:
        print_to_console (bool): Whether to print logs to console. Defaults to True.
        console_level (Optional[str]): Log level for console output.
        save_to_file (bool): Whether to save logs to a file.
        file_level (Optional[str]): Log level for file output.
        log_file (Optional[Path]): Path to the log file.

    Returns:
        GDTLogger: Configured logger instance.

    """
    log = cast(GDTLogger, logging.getLogger("gdt"))
    levels = []
    # Remove any existing handlers
    for handler in log.handlers[:]:
        log.removeHandler(handler)

    if print_to_console:
        console_level_int = _logging_levels.get(console_level, logging.INFO)
        console_handler = logging.StreamHandler()
        console_handler.setLevel(console_level_int)
        console_formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s"
        )
        console_handler.setFormatter(console_formatter)
        log.addHandler(console_handler)
        levels.append(console_level_int)

    if save_to_file:
        if log_file is None:
            print(
                "No log file specified, even though save_to_file is True. "
                "Log will be save in the current directory in 'gdt_default.log'."
            )
            log_file = "gdt_default.log"

        log_file_path = Path(log_file).resolve()
        log_file_path.touch(exist_ok=True)

        file_level_int = _logging_levels.get(file_level, logging.DEBUG)
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(file_level_int)
        file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(file_formatter)
        log.addHandler(file_handler)
        levels.append(file_level_int)

    log.setLevel(min(levels) if levels else _logging_levels["DISABLE"])
    log.propagate = False
    log.debug("Simple log setup complete.")
    log.debug(f"Console logging level {console_level if print_to_console else 'None'}")
    log.debug(
        f"File logging level {file_level} at {log_file_path if log_file else 'None'}"
    )

    return log


def setup_logger(
    debug: bool,
    log_file: Union[Path, str, None],
    quiet: bool,
) -> GDTLogger:
    """Set up logger based on command line arguments."""
    console_level = "DISABLE" if quiet else ("DEBUG" if debug else "INFO")
    file_level = "TRACE" if debug else "DEBUG"

    # Create logger based on log file preference
    if log_file:
        log_file = Path(log_file).resolve()
        log = create_simple_logger(
            print_to_console=not quiet,
            console_level=console_level,
            file_level=file_level,
            log_file=log_file,
        )
    else:
        log = create_dev_logger(console_level=console_level, file_level=file_level)

    log.trace("Logger setup complete.")
    return log
