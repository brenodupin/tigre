# -*- coding: utf-8 -*-
"""Logger setup and management for tigre.

The base code for this module (log_setup.py from GDT package) has been
adapted for the TIGRE project.

Source
https://github.com/brenodupin/gdt/blob/master/src/gdt/log_setup.py
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any, cast

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

    def setLevel(self, level: int | str) -> None:  # noqa: N802
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


def create_simple_logger(
    print_to_console: bool = True,
    console_level: str = "INFO",
    save_to_file: bool = True,
    file_level: str = "DEBUG",
    log_file: Path | str | None = None,
) -> GDTLogger:
    """Create a simple logger with optional console and file output.

    Args:
        print_to_console (bool): Whether to print logs to console. Defaults to True.
        console_level (str): Log level for console output.
        save_to_file (bool): Whether to save logs to a file.
        file_level (str): Log level for file output.
        log_file (Path | str | None): Path to the log file.

    Returns:
        GDTLogger: Configured logger instance.

    """
    log = cast(GDTLogger, logging.getLogger("tigre"))
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
                "Log will be saved in the current directory in 'tigre_default.log'."
            )
            log_file = "tigre_default.log"

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
    log_file: Path | str | None,
    quiet: bool,
    no_log_file: bool,
    command: str,
) -> GDTLogger:
    """Set up logger based on command line arguments."""
    console_level = "DISABLE" if quiet else ("DEBUG" if debug else "INFO")
    file_level = "DISABLE" if no_log_file else ("TRACE" if debug else "DEBUG")

    if log_file is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = Path(f"tigre_{command}_{timestamp}.log")

    log_file = Path(log_file).resolve()
    log = create_simple_logger(
        print_to_console=not quiet,
        console_level=console_level,
        save_to_file=not no_log_file,
        file_level=file_level,
        log_file=log_file,
    )

    log.trace("Logger setup complete.")
    return log
