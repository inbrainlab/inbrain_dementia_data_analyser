import numpy as np
import pandas as pd
import pytest

from inbrain_dementia_data_analyser.config import SourceSpec
from inbrain_dementia_data_analyser.series import build_scored_columns
from inbrain_dementia_data_analyser.variables import materialise


def src(**kw):
    base = dict(variable="V", dataset="ADNI", table="t")
    base.update(kw)
    for k in ("column", "items", "missing_values"):
        if k in base and not isinstance(base[k], tuple):
            base[k] = tuple(base[k])
    if "valid_codes" in base:
        base["valid_codes"] = tuple(base["valid_codes"])
    return SourceSpec(**base)


def test_plain_column_is_used_as_is():
    df = pd.DataFrame({"A": [1, 2]})
    frame, col = materialise(df, "V", src(column=("A",)))
    assert frame is df and col == "A"


def test_first_candidate_column_present_wins():
    df = pd.DataFrame({"VSBPDIAS": [80]})
    _, col = materialise(df, "BP_Dias", src(column=("VSBPDIA", "VSBPDIAS")))
    assert col == "VSBPDIAS"


def test_missing_column_returns_none():
    assert materialise(pd.DataFrame({"A": [1]}), "V", src(column=("B",))) == (None, None)


def test_sum_transform_needs_one_item():
    df = pd.DataFrame({"I1": [1, np.nan], "I2": [2, np.nan]})
    frame, col = materialise(df, "TOTAL", src(transform="sum", items=("I1", "I2")))
    assert frame[col].tolist()[0] == 3 and np.isnan(frame[col].tolist()[1])


def test_coalesce_uds_v3_then_v2():
    df = pd.DataFrame({"TRAILA": [30, np.nan], "OTRAILA": [99, 45]})
    frame, col = materialise(df, "TMT_A", src(transform="coalesce", items=("TRAILA", "OTRAILA")))
    assert frame[col].tolist() == [30, 45]


def test_cleaning_rules():
    df = pd.DataFrame({"S": [0, 1, 2, 9, -4], "W": [150, 888, 0, 9999, 120]})
    f, c = materialise(df, "SEIZ", src(column=("S",), valid_codes=(0, 1, 2)))
    assert f[c].tolist()[:3] == [0, 1, 2] and f[c].isna().sum() == 2
    f, c = materialise(df, "W", src(column=("W",), missing_values=(888, 9999), valid_range=(0, None),
                                    range_inclusive=False))
    assert f[c].tolist()[0] == 150 and f[c].isna().sum() == 3


def test_row_filter_and_case_insensitive():
    df = pd.DataFrame({"PACCNAME": ["PACC composite", "Digit"], "paccrn": [0.5, 3]})
    f, c = materialise(df, "PACC", src(column=("PACCRN",), case_insensitive=True,
                                       row_filter={"column": "PACCNAME", "pattern": "composite"}))
    assert len(f) == 1 and f[c].iloc[0] == 0.5


def test_bmi_transform_units():
    df = pd.DataFrame({"W": [70, 154.3], "H": [175, 68.9], "WU": [2, 1], "HU": [2, 1]})
    f, c = materialise(df, "BMI", src(transform="bmi", params={
        "weight": "W", "height": "H", "weight_unit_column": "WU", "height_unit_column": "HU"}))
    assert f[c].round(0).tolist() == [23.0, 23.0]


def test_scored_columns_stay_aligned():
    s = {"S1": [(pd.Timestamp("2020-01-01"), 28), (pd.Timestamp("2021-01-01"), np.nan),
                (pd.Timestamp("2022-01-01"), 26)]}
    out = build_scored_columns(s, "MMSE").loc["S1"]
    assert out["all_MMSE_date"] == "2020-01-01|2022-01-01"
    assert out["all_MMSE_score"] == "28|26"
    assert out["nb_of_MMSE_visits"] == 2 and out["Lst_MMSE_score"] == 26


def test_item_codes_excluded_before_summing():
    """ECog: 9 = "don't know" must not enter the total."""
    df = pd.DataFrame({"I1": [1, 9], "I2": [4, 2]})
    f, c = materialise(df, "ECOG", src(transform="sum", items=("I1", "I2"),
                                       params={"item_missing_values": [9]}))
    assert f[c].tolist() == [5, 2]
