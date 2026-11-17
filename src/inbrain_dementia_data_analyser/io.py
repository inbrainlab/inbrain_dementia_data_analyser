"""Reading tables safely and with bounded memory."""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

EMPTY_TOKENS = ("", "nan", "none", "nat", "na")

_CSV_STRATEGIES = (
    dict(engine="c", on_bad_lines="error", encoding="utf-8"),
    dict(engine="c", on_bad_lines="skip", encoding="utf-8"),
    dict(engine="python", on_bad_lines="skip", encoding="utf-8"),
    dict(engine="python", on_bad_lines="skip", encoding="latin-1"),
)


class TableReadError(IOError):
    """Raised when a table exists but cannot be parsed."""


def read_header(path: str | os.PathLike, sheet: str | int | None = None) -> list[str]:
    """Column names of a CSV or Excel file, without reading its rows."""
    path = Path(path)
    if path.suffix.lower() in (".xlsx", ".xls"):
        return [str(c) for c in pd.read_excel(path, sheet_name=sheet or 0, nrows=0).columns]
    last = None
    for enc in ("utf-8", "latin-1"):
        try:
            return [str(c) for c in pd.read_csv(path, nrows=0, encoding=enc).columns]
        except (UnicodeDecodeError, pd.errors.ParserError) as exc:
            last = exc
    raise TableReadError(f"Cannot read the header of {path}: {last}")


def read_table(path: str | os.PathLike, *, sheet: str | int | None = None,
               usecols: Sequence[str] | None = None, dtype=None) -> pd.DataFrame:
    """Read a CSV (with fallbacks for malformed lines and encodings) or an Excel sheet.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    TableReadError
        If every parsing strategy fails.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Table not found: {path}")
    if path.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(path, sheet_name=sheet or 0, usecols=usecols, dtype=dtype)
    last = None
    for i, kw in enumerate(_CSV_STRATEGIES):
        try:
            extra = {} if kw["engine"] == "python" else {"low_memory": False}
            df = pd.read_csv(path, usecols=usecols, dtype=dtype, **kw, **extra)
            if i > 0:
                log.warning("  [WARN] %s: read with fallback strategy %d (%s)", path.name, i + 1, kw)
            return df
        except (UnicodeDecodeError, pd.errors.ParserError, ValueError) as exc:
            last = exc
    raise TableReadError(f"Cannot parse {path}: {last}")


def read_csv_lean(path: str | os.PathLike, usecols: Iterable[str], *,
                  text_cols: Iterable[str] | None = None, chunksize: int = 4000) -> pd.DataFrame:
    """Read selected columns of a very wide CSV in row chunks.

    ``pd.read_csv(..., low_memory=False)`` tokenises the whole file before it
    applies ``usecols``; on the expanded master that costs gigabytes even for a
    few columns. Reading in chunks keeps memory proportional to the columns
    used. Every column is read as text; columns whose non-empty values are all
    numeric are converted back to numbers, except ``text_cols``.
    """
    usecols = list(dict.fromkeys(usecols))
    keep_text = set(text_cols or ())
    parts = [chunk for chunk in pd.read_csv(path, usecols=usecols, dtype=str, chunksize=chunksize)]
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=usecols)
    for c in df.columns:
        if c in keep_text:
            continue
        s = df[c]
        nn = s.dropna()
        if nn.empty:
            continue
        num = pd.to_numeric(s, errors="coerce")
        if num.notna().sum() == len(nn):
            df[c] = num
    return df


def is_empty_token(x: object) -> bool:
    """True for NaN/None and for the text placeholders used inside lists."""
    if x is None:
        return True
    if isinstance(x, float) and np.isnan(x):
        return True
    return str(x).strip().lower() in EMPTY_TOKENS


def file_fingerprint(path: str | os.PathLike) -> dict:
    """Name, size and modification time of a file (for provenance records)."""
    p = Path(path)
    st = p.stat()
    return {"path": str(p), "name": p.name, "bytes": st.st_size,
            "modified": pd.Timestamp(st.st_mtime, unit="s").isoformat()}
