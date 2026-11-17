"""Regression tests for the data bugs fixed after v0.1.0 (see docs/REFACTOR_REPORT.md)."""
import numpy as np
import pandas as pd
import yaml

from inbrain_dementia_data_analyser.config import SourceSpec, load_config
from inbrain_dementia_data_analyser.datasets import a4
from inbrain_dementia_data_analyser.variables import materialise


def src(**kw):
    base = dict(variable="V", dataset="NIFD", table="t")
    base.update(kw)
    for k in ("column", "items", "missing_values", "raw_missing_values"):
        if k in base and not isinstance(base[k], tuple):
            base[k] = tuple(base[k])
    return SourceSpec(**base)


# --- missing codes ---------------------------------------------------------

def test_dataset_raw_codes_never_reach_a_score():
    df = pd.DataFrame({"MMSE_TOT": [28, -5, 30, -9]})
    f, c = materialise(df, "MMSE", src(column=("MMSE_TOT",), raw_missing_values=(-5, -9)))
    assert f[c].tolist()[0] == 28 and f[c].isna().tolist() == [False, True, False, True]


def test_raw_codes_are_masked_in_transform_items_not_in_the_result():
    """A difference may be negative; only its inputs are cleaned."""
    df = pd.DataFrame({"A": [5, -1, 3], "B": [7, 2, 3]})
    f, c = materialise(df, "D", src(transform="difference", items=("A", "B"),
                                    raw_missing_values=(-1, -4)))
    assert f[c].tolist()[0] == -2 and np.isnan(f[c].tolist()[1]) and f[c].tolist()[2] == 0


def test_nifd_npi_logical_skip_counts_as_zero():
    """NIFD NPI severity -6 = symptom absent -> 0; -9 -> missing; 12 items needed."""
    items = [f"I{i}" for i in range(12)]
    df = pd.DataFrame([[-6] * 12, [-6] * 10 + [2, 3], [-9] + [-6] * 11], columns=items)
    s = src(transform="sum", items=tuple(items), raw_missing_values=tuple(range(-9, 0)),
            params={"item_recode": {-6: 0}, "min_count": 12})
    f, c = materialise(df, "NPI", s)
    assert f[c].tolist()[:2] == [0, 5] and np.isnan(f[c].tolist()[2])


def test_item_valid_range_drops_nacc_item_codes():
    df = pd.DataFrame({f"REY{i}REC": [5, 5] for i in range(1, 6)})
    df.loc[1, "REY3REC"] = 95
    s = src(transform="sum", items=tuple(df.columns),
            params={"min_count": 5, "item_valid_range": [0, 15]})
    f, c = materialise(df, "RAVLT_I", s)
    assert f[c].tolist()[0] == 25 and np.isnan(f[c].tolist()[1])


def test_variable_valid_range_is_inherited_and_can_be_switched_off(config_copy):
    p = config_copy / "variables.yaml"
    data = yaml.safe_load(p.read_text())
    data["variables"]["MMSE"]["sources"]["A4"] = {"table": "mmse", "column": "MMSCORE",
                                                  "valid_range": None}
    p.write_text(yaml.safe_dump(data, sort_keys=False))
    cfg = load_config(config_copy)
    by_ds = {ds: s[0] for ds, s in cfg.variables["MMSE"].sources.items()}
    assert by_ds["NACC"].valid_range == (0, 30)
    assert by_ds["A4"].valid_range is None


def test_dataset_codes_reach_the_sources(cfg):
    nifd = cfg.variables["MMSE"].sources["NIFD"][0]
    adni = cfg.variables["MMSE"].sources["ADNI"][0]
    assert -5 in nifd.raw_missing_values and set(adni.raw_missing_values) == {-1, -4}
    assert cfg.variables["MMSE"].sources["A4"][0].raw_missing_values == ()


# --- variable mapping ------------------------------------------------------

def test_nifd_moca_is_the_moca_total(cfg):
    s = cfg.variables["MoCA"].sources["NIFD"][0]
    assert s.column == ("MOCA_TotWithEduc",) and s.date_column == "MOCA_DCDATE"
    assert cfg.variables["Benson_copy"].sources["NIFD"][0].column == ("MODREY",)


def test_ravlt_measures_are_not_pooled(cfg):
    v = cfg.variables
    assert v["RAVLT_T1"].sources["ADNI"][0].column == ("AVTOT1",)
    assert v["RAVLT_I"].sources["ADNI"][0].items == tuple(f"AVTOT{i}" for i in range(1, 6))
    assert v["RAVLT_I"].sources["NACC"][0].items == tuple(f"REY{i}REC" for i in range(1, 6))
    assert v["RAVLT_I"].sources["PREVENTAD"][0].column == ("List_A_Sum_Trials_1_5",)
    assert v["RAVLT_D"].sources["ADNI"][0].column == ("AVDEL30MIN",)
    assert v["RAVLT_D"].sources["NACC"][0].column == ("REYDREC",)
    assert v["RAVLT_Recog"].sources["ADNI"][0].column == ("AVDELTOT",)


# --- A4 days ---------------------------------------------------------------

def test_a4_unknown_examday_stays_missing():
    reg = pd.DataFrame({"BID": ["B1", "B1"], "VISCODE": [1, 3], "EXAMDAY": [0, 50]})
    df = pd.DataFrame({"BID": ["B1", "B1", "B2"], "VISCODE": [1, 4, 4], "X": [1, 2, 3]})
    out = a4.attach_examday(df, reg)
    assert out["_examday"].tolist()[0] == 0
    assert out["_examday"].isna().tolist() == [False, True, True]     # no 2970 / 0 placeholder


def test_a4_day_lists_keep_unknown_days_aligned():
    assert a4._na_list([0.0, np.nan, 50.0]) == "0.0|NA|50.0"
    assert a4._na_list([1.02, np.nan]) == "1.02|NA"


def test_a4_not_done_registry_visits_are_not_counted():
    reg = pd.DataFrame({"BID": ["B1"] * 3 + ["B2"], "VISCODE": [1, 2, 3, 1],
                        "VISTYPE": [1, 0, 1, 0], "EXAMDAY": [0, np.nan, 40, np.nan]})
    a4._CTX = None
    out = a4.build_diagnosis(reg)
    assert out.loc["B1", "all_diag_date"] == "0.0|40.0" and out.loc["B1", "nb_of_visits"] == 2
    assert out.loc["B2", "1st_diag"] == "CU" and out.loc["B2", "nb_of_visits"] == 0
