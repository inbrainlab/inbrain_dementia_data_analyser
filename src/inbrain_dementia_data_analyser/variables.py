"""Applying ``variables.yaml`` to raw tables.

Dataset modules ask for the variables configured for one of their tables
(:func:`variable_sources`) and call :func:`materialise` to obtain the column
that holds the cleaned values. The rest of the series construction (subject ID
normalisation, dates, sort order) stays with the dataset module, because each
source has its own visit and date conventions.

Tables declared with ``generic: true`` in a dataset file are not known to the
dataset module; :func:`build_generic_variables` reads them and builds every
variable that points to them, so a new instrument can be added without code.
"""
from __future__ import annotations

import logging
import re
from typing import Iterator

import numpy as np
import pandas as pd

from .config import ProjectConfig, SourceSpec, load_config
from .transforms import TRANSFORMS

log = logging.getLogger(__name__)

HELPER_PREFIX = "__var__"


def variable_sources(dataset: str, table: str | None = None,
                     config: ProjectConfig | None = None) -> list[tuple[str, SourceSpec]]:
    """``(variable_name, source)`` pairs configured for ``dataset`` (and ``table``)."""
    cfg = config or load_config()
    return cfg.sources(dataset, table)


def resolve_column(df: pd.DataFrame, src: SourceSpec) -> str | None:
    """First candidate column of ``src`` present in ``df`` (case-insensitive if asked)."""
    lower = {c.lower(): c for c in df.columns} if src.case_insensitive else None
    for c in src.column:
        if c in df.columns:
            return c
        if lower is not None and c.lower() in lower:
            return lower[c.lower()]
    return None


def clean_values(values: pd.Series, src: SourceSpec) -> pd.Series:
    """Apply ``valid_codes``, ``missing_values`` and ``valid_range`` of a source."""
    out = values
    if src.valid_codes is not None:
        num = pd.to_numeric(out, errors="coerce")
        out = num.where(num.isin(list(src.valid_codes)))
    if src.missing_values:
        num = pd.to_numeric(out, errors="coerce")
        miss = num.isin([v for v in src.missing_values if not isinstance(v, str)])
        txt = out.astype(str).str.strip().isin([str(v) for v in src.missing_values])
        out = out.where(~(miss | txt))
    if src.valid_range is not None:
        lo, hi = src.valid_range
        num = pd.to_numeric(out, errors="coerce")
        ok = pd.Series(True, index=out.index)
        if lo is not None:
            ok &= (num >= lo) if src.range_inclusive else (num > lo)
        if hi is not None:
            ok &= (num <= hi) if src.range_inclusive else (num < hi)
        out = num.where(ok)
    return out


def materialise(df: pd.DataFrame, name: str, src: SourceSpec) -> tuple[pd.DataFrame | None, str | None]:
    """Return ``(frame, column)`` holding the values of variable ``name``.

    * Plain sources return the table itself and the raw column (no copy).
    * Sources with a ``transform`` or cleaning rules write the result to a
      helper column ``__var__<name>`` of ``df`` (in place) and return it.
    * Sources with a ``row_filter`` return a filtered view of the table.

    ``(None, None)`` means the variable cannot be built from this table; a
    warning names the missing column(s).
    """
    if df is None or df.empty:
        return None, None
    frame = df
    if src.row_filter:
        col = src.row_filter["column"]
        if col not in frame.columns:
            log.warning("    [WARN] %s: row_filter column %r not found -> skipped", name, col)
            return None, None
        frame = frame[frame[col].astype(str).str.contains(src.row_filter["pattern"], case=False,
                                                          na=False, regex=True)]
        if frame.empty:
            return None, None
    if src.transform:
        values = TRANSFORMS[src.transform](frame, src)
        if values is None:
            log.warning("    [WARN] %s: inputs of transform %r not found (%s) -> skipped",
                        name, src.transform, ", ".join(src.raw_columns()) or "-")
            return None, None
    else:
        col = resolve_column(frame, src)
        if col is None:
            log.warning("    [WARN] '%s' not found -> %s skipped", "' / '".join(src.column), name)
            return None, None
        if not src.needs_preparation:
            return frame, col
        values = frame[col]
    helper = f"{HELPER_PREFIX}{name}"
    values = clean_values(values, src)
    if frame is df:
        df[helper] = values
        return df, helper
    frame = frame.copy()
    frame[helper] = values
    return frame, helper


# ---------------------------------------------------------------------------
# generic tables
# ---------------------------------------------------------------------------
def generic_table_keys(dataset: str, config: ProjectConfig | None = None) -> list[str]:
    """Tables of ``dataset`` that are read by the generic engine."""
    cfg = config or load_config()
    ds = cfg.dataset(dataset)
    return [k for k, t in ds.tables.items() if t.generic]


def iter_generic_series(dataset: str, tables: dict[str, pd.DataFrame],
                        config: ProjectConfig | None = None) -> Iterator[pd.DataFrame]:
    """Scored series for every variable that points to a generic table."""
    from .series import scored_series_from_long

    cfg = config or load_config()
    ds = cfg.dataset(dataset)
    for key in generic_table_keys(ds.name, cfg):
        spec = ds.tables[key]
        df = tables.get(key)
        if df is None or df.empty:
            continue
        if spec.id_column not in df.columns or spec.date_column not in df.columns:
            log.warning("    [WARN] generic table %s: id/date column %r/%r missing -> skipped",
                        key, spec.id_column, spec.date_column)
            continue
        for name, src in cfg.sources(ds.name, key):
            frame, col = materialise(df, name, src)
            if col is None:
                continue
            date_col = src.date_column or spec.date_column
            if ds.time_axis in ("days", "months", "visit"):
                parse = False
            else:
                parse = True
            part = scored_series_from_long(frame, spec.id_column, date_col, col, name,
                                           parse_dates=parse)
            if not part.empty:
                yield part


def attach_generic_variables(master: pd.DataFrame, dataset: str, tables: dict[str, pd.DataFrame],
                             config: ProjectConfig | None = None, id_col: str = "ID") -> pd.DataFrame:
    """Left-join the series of generic tables onto a per-subject table.

    Subject IDs are compared as text and the master's ID values are kept. The
    number of subjects present only in the generic tables is logged (they are
    not added, so the cohort's subject list stays the one defined by its
    clinical tables).
    """
    parts = list(iter_generic_series(dataset, tables, config))
    if not parts:
        return master
    extra = parts[0]
    for p in parts[1:]:
        extra = extra.join(p, how="outer")
    extra.index = extra.index.map(lambda x: str(x).strip())
    key = master[id_col].map(lambda x: str(x).strip())
    known = set(key)
    dropped = sum(1 for i in extra.index if i not in known)
    if dropped:
        log.warning("    [WARN] %s generic variables: %d subjects not in the cohort table were ignored",
                    dataset, dropped)
    clash = [c for c in extra.columns if c in master.columns]
    if clash:
        log.warning("    [WARN] %s generic variables overwrite existing columns: %s", dataset, clash[:8])
        master = master.drop(columns=clash)
    joined = extra.reindex(key.values)
    joined.index = master.index
    return pd.concat([master, joined], axis=1)


def load_generic_tables(dataset: str, ctx) -> dict[str, pd.DataFrame]:
    """Read the generic tables of a dataset through its run context."""
    from .io import read_table

    out = {}
    for key in generic_table_keys(dataset, ctx.config):
        path = ctx.path(key)
        if path is None:
            continue
        spec = ctx.dataset.tables[key]
        out[key] = read_table(path, sheet=spec.sheet)
        log.info("  [OK]   %s (generic): %d rows", key, len(out[key]))
    return out
