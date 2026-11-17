"""Named transforms that variables in ``variables.yaml`` can refer to.

A transform receives the raw table and the :class:`~.config.SourceSpec` and
returns one value per row (a :class:`pandas.Series` aligned with the table), or
``None`` when none of the columns it needs is present.

Register a new transform with the :func:`register` decorator::

    from inbrain_dementia_data_analyser.transforms import register

    @register("my_score")
    def my_score(df, src):
        return df[src.items[0]] * 2

and refer to it in ``variables.yaml`` with ``transform: my_score``.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import pandas as pd

TransformFn = Callable[[pd.DataFrame, "object"], "pd.Series | None"]
TRANSFORMS: dict[str, TransformFn] = {}


def register(name: str) -> Callable[[TransformFn], TransformFn]:
    """Decorator that makes a transform available to ``variables.yaml``."""

    def deco(fn: TransformFn) -> TransformFn:
        TRANSFORMS[name] = fn
        return fn

    return deco


def _present(df: pd.DataFrame, cols) -> list[str]:
    return [c for c in cols if c in df.columns]


def _numeric(df: pd.DataFrame, cols, src=None) -> pd.DataFrame:
    """Items as numbers; codes in ``params.item_missing_values`` (e.g. 9 =
    "don't know") become missing before the items are combined."""
    num = df[cols].apply(pd.to_numeric, errors="coerce")
    miss = (src.params.get("item_missing_values") if src is not None else None) or ()
    if miss:
        num = num.mask(num.isin(list(miss)))
    return num


@register("sum")
def t_sum(df: pd.DataFrame, src) -> pd.Series | None:
    """Row sum of ``items``; ``params.min_count`` non-missing items needed (default 1);
    ``params.item_missing_values`` lists item codes treated as missing."""
    cols = _present(df, src.items)
    if not cols:
        return None
    return _numeric(df, cols, src).sum(axis=1, min_count=int(src.params.get("min_count", 1)))


@register("mean")
def t_mean(df: pd.DataFrame, src) -> pd.Series | None:
    """Row mean of ``items``; ``params.min_count`` non-missing items needed (default 1);
    ``params.item_missing_values`` lists item codes treated as missing."""
    cols = _present(df, src.items)
    if not cols:
        return None
    num = _numeric(df, cols, src)
    out = num.mean(axis=1)
    return out.where(num.notna().sum(axis=1) >= int(src.params.get("min_count", 1)))


@register("max")
def t_max(df: pd.DataFrame, src) -> pd.Series | None:
    """Row maximum of ``items`` (e.g. 'any medication' from 0/1 flags)."""
    cols = _present(df, src.items)
    if not cols:
        return None
    return _numeric(df, cols, src).max(axis=1, skipna=True)


@register("coalesce")
def t_coalesce(df: pd.DataFrame, src) -> pd.Series | None:
    """First non-missing value across ``items`` in the given order."""
    cols = _present(df, src.items)
    if not cols:
        return None
    out = df[cols[0]]
    for c in cols[1:]:
        out = out.combine_first(df[c])
    return out


@register("difference")
def t_difference(df: pd.DataFrame, src) -> pd.Series | None:
    """``items[0] - items[1]`` (e.g. RAVLT forgetting = trial 5 - delayed recall)."""
    if len(src.items) != 2 or not set(src.items) <= set(df.columns):
        return None
    num = _numeric(df, list(src.items))
    return num[src.items[0]] - num[src.items[1]]


@register("ratio")
def t_ratio(df: pd.DataFrame, src) -> pd.Series | None:
    """``items[0] / items[1]``; division by zero gives NaN."""
    if len(src.items) != 2 or not set(src.items) <= set(df.columns):
        return None
    num = _numeric(df, list(src.items))
    den = num[src.items[1]].replace(0, np.nan)
    return num[src.items[0]] / den


@register("bmi")
def t_bmi(df: pd.DataFrame, src) -> pd.Series | None:
    """Body-mass index from weight and height.

    ``params``: ``weight`` and ``height`` (column names), and either fixed
    units (``weight_unit: kg|lb``, ``height_unit: cm|in|m``) or unit-code
    columns (``weight_unit_column`` / ``height_unit_column``; codes 1 = lb/in,
    2 = kg/cm). Uses :func:`~.clinical.compute_bmi`.
    """
    from .clinical import compute_bmi, find_col

    p = src.params
    w, h = p.get("weight"), p.get("height")
    if src.case_insensitive:
        w, h = find_col(df, w) if w else None, find_col(df, h) if h else None
    if not w or not h or w not in df.columns or h not in df.columns:
        return None
    wu = df[p["weight_unit_column"]] if p.get("weight_unit_column") in df.columns else p.get("weight_unit", "kg")
    hu = df[p["height_unit_column"]] if p.get("height_unit_column") in df.columns else p.get("height_unit", "cm")
    return compute_bmi(df[w], df[h], weight_unit=wu, height_unit=hu)


@register("bmi_reported_or_computed")
def t_bmi_reported(df: pd.DataFrame, src) -> pd.Series | None:
    """Reported BMI (``params.bmi_column``) when plausible (10-80), else computed (``bmi``)."""
    calc = t_bmi(df, src)
    col = src.params.get("bmi_column")
    if calc is None:
        calc = pd.Series(np.nan, index=df.index)
    if col and col in df.columns:
        rep = pd.to_numeric(df[col], errors="coerce")
        rep = rep.where(rep.between(10, 80))
        return rep.combine_first(calc)
    return calc


@register("sum_timepoint_items")
def t_sum_timepoint_items(df: pd.DataFrame, src) -> pd.Series | None:
    """Wide timepoint tables only (4RTNI): sum of ``<instrument>_<tp>_<item>`` per timepoint.

    The 4RTNI pipeline evaluates this transform itself, one timepoint at a
    time; on an ordinary long table it is not meaningful and returns ``None``.
    """
    return None
