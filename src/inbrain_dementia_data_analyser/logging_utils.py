"""Logging helpers.

The package never configures the root logger on import. Call
:func:`enable_logging` (or configure :mod:`logging` yourself) to see progress
messages. The ported pipelines report through :func:`make_print`, which turns
``print``-style calls into ``logger.info`` records.
"""
from __future__ import annotations

import logging
from typing import Callable

PACKAGE_LOGGER = "inbrain_dementia_data_analyser"


def enable_logging(level: int | str = "INFO", fmt: str = "%(message)s") -> logging.Logger:
    """Send the package's log records to stderr.

    Parameters
    ----------
    level:
        Logging level name or number (``"INFO"`` shows progress, ``"WARNING"``
        only problems).
    fmt:
        Format string for the handler.

    Returns
    -------
    logging.Logger
        The package logger.
    """
    logger = logging.getLogger(PACKAGE_LOGGER)
    if isinstance(level, str):
        level = getattr(logging, level.upper())
    logger.setLevel(level)
    if not any(getattr(h, "_inbrain", False) for h in logger.handlers):
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter(fmt))
        handler._inbrain = True  # type: ignore[attr-defined]
        logger.addHandler(handler)
    logger.propagate = False
    return logger


def make_print(logger: logging.Logger) -> Callable[..., None]:
    """Return a ``print``-compatible function that writes to ``logger``.

    Lines that start with ``[WARN]``, ``[ERROR]`` or ``[SKIP]`` (after leading
    spaces) are logged at WARNING/ERROR level so that they are visible even
    when INFO messages are hidden.
    """

    def _print(*args: object, sep: str = " ", end: str = "\n", **_: object) -> None:
        msg = sep.join(str(a) for a in args)
        head = msg.lstrip()
        if head.startswith("[ERROR]"):
            logger.error(msg)
        elif head.startswith(("[WARN]", "[SKIP]", "[ACTION]")):
            logger.warning(msg)
        else:
            logger.info(msg)

    return _print


def quiet_pandas_warnings(fn: Callable) -> Callable:
    """Silence routine pandas warnings (downcasting, chained assignment, date
    inference) while ``fn`` runs; the notebooks silenced them globally."""
    import functools
    import warnings

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        import pandas as pd

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=FutureWarning)
            warnings.simplefilter("ignore", category=UserWarning)
            warnings.simplefilter("ignore", category=pd.errors.SettingWithCopyWarning)
            warnings.simplefilter("ignore", category=pd.errors.PerformanceWarning)
            return fn(*args, **kwargs)

    return wrapper
