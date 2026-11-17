"""Per-subject longitudinal series in the stage-1 output format.

A scored measure ``X`` is written as::

    all_X_date | all_X_score      '|'-joined, chronological, aligned
    1st_X_score, 1st_X_date       first non-missing value
    Lst_X_score, Lst_X_date       last non-missing value
    nb_of_X_visits                number of non-missing values
"""
from __future__ import annotations

from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


def safe_date(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """Parse ``df[col]`` to datetime in place (invalid values -> NaT); no-op if absent."""
    if col in df.columns:
        df[col] = pd.to_datetime(df[col], errors="coerce")
    return df


def list_str(values: Iterable[Any]) -> str | float:
    """Pipe-separated string of the non-missing values; NaN when nothing is left."""
    clean = [str(v) for v in values
             if pd.notna(v) and str(v).strip() not in ("", "nan", "NaT")]
    return "|".join(clean) if clean else np.nan


def days_from(d1: Any, d2: Any) -> float:
    """``(d2 - d1).days`` as a float; NaN when either date is missing."""
    if pd.isna(d1) or pd.isna(d2):
        return np.nan
    try:
        return float((d2 - d1).days)
    except Exception:  # noqa: BLE001 - mixed types: keep the legacy NaN result
        return np.nan


def d_fmt(d: Any) -> str | float:
    """Format a date as ``YYYY-MM-DD``; numbers (day offsets, visits) as text; NaN stays NaN."""
    if pd.isna(d):
        return np.nan
    return str(d.date()) if hasattr(d, "date") else str(d)


def _d_fmt_text(d: Any) -> str:
    """ADNI/AIBL variant: missing dates become the text 'NaT'/'nan'."""
    return str(d.date()) if hasattr(d, "date") and pd.notna(d) else str(d)


def build_scored_columns(scores_by_subj: Mapping[Any, Sequence[tuple[Any, Any]]],
                         test_name: str, *, missing_date_as_text: bool = False) -> pd.DataFrame:
    """Turn ``{subject: [(date, value), ...]}`` into the scored-family columns.

    Parameters
    ----------
    scores_by_subj:
        Chronologically ordered ``(date, value)`` pairs per subject. Values that
        are missing are dropped; their dates are dropped with them, so the date
        and score lists stay aligned.
    test_name:
        Output stem ``X``.
    missing_date_as_text:
        Keep the ADNI/AIBL behaviour where a missing first/last date is written
        as the text ``NaT`` instead of an empty cell.

    Returns
    -------
    pandas.DataFrame
        Indexed by subject ID (empty frame when there is nothing to write).
    """
    fmt = _d_fmt_text if missing_date_as_text else d_fmt
    rows = []
    for subj, visits in scores_by_subj.items():
        clean = [(d, v) for d, v in visits if pd.notna(v)]
        if not clean and not visits:
            continue
        rows.append({
            "ID": subj,
            f"all_{test_name}_date": ("|".join("NA" if pd.isna(d) else str(fmt(d)) for d, _ in clean)
                                       or np.nan),
            f"all_{test_name}_score": ("|".join(str(v) for _, v in clean) or np.nan),
            f"1st_{test_name}_score": clean[0][1] if clean else np.nan,
            f"1st_{test_name}_date": fmt(clean[0][0]) if clean else np.nan,
            f"Lst_{test_name}_score": clean[-1][1] if clean else np.nan,
            f"Lst_{test_name}_date": fmt(clean[-1][0]) if clean else np.nan,
            f"nb_of_{test_name}_visits": len(clean),
        })
    return pd.DataFrame(rows).set_index("ID") if rows else pd.DataFrame()


def scored_series_from_long(df: pd.DataFrame, id_col: str, date_col: str, value_col: str,
                            name: str, *, parse_dates: bool = True,
                            missing_date_as_text: bool = False) -> pd.DataFrame:
    """Build the scored family from a long table (one row per visit)."""
    if value_col not in df.columns:
        return pd.DataFrame()
    sub = df[[id_col, date_col, value_col]].copy()
    if parse_dates:
        safe_date(sub, date_col)
    sub = sub.sort_values([id_col, date_col])
    scores = {subj: list(zip(g[date_col], g[value_col])) for subj, g in sub.groupby(id_col)}
    return build_scored_columns(scores, name, missing_date_as_text=missing_date_as_text)
