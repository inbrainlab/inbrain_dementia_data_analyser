"""v0.4.0: diagnosis as two axes (cognitive stage + primary etiology)."""
import numpy as np
import pandas as pd

from inbrain_dementia_data_analyser import diagnosis as dx
from inbrain_dementia_data_analyser.datasets import adni, aibl, mcsa, nifd, preventad


# ----------------------------------------------------------------- labels
def test_harmonised_label():
    assert dx.harmonised_label("CN", "none") == "CN"
    assert dx.harmonised_label("MCI", "AD") == "MCI"            # etiology kept apart
    assert dx.harmonised_label("IMPNOMCI", "unknown") == "IMPNOMCI"
    assert dx.harmonised_label("DEM", "AD") == "AD"
    assert dx.harmonised_label("DEM", "unknown") == "Dementia-NOS"
    assert dx.harmonised_label("DEM", "Other") == "Dementia-Other"
    assert dx.harmonised_label("DEM", "AD+VaD") == "AD+VaD"
    assert dx.harmonised_label("Unknown", "PSP") == "PSP"
    assert dx.harmonised_label("Unknown", "unknown") == "Unknown"


def test_join_etiologies():
    assert dx.join_etiologies(["VaD", "AD", "AD"]) == "AD+VaD"
    assert dx.join_etiologies([]) == "unknown"


# ----------------------------------------------------------------- NACC
def test_nacc_uses_naccudsd_and_naccetpr():
    row = pd.Series({"NACCUDSD": 3, "NACCETPR": 1})
    assert dx.nacc_stage_etiology(row)[:2] == ("MCI", "AD")       # MCI due to AD stays MCI
    row = pd.Series({"NACCUDSD": 4, "NACCETPR": 8})
    assert dx.nacc_stage_etiology(row)[:2] == ("DEM", "VaD")
    row = pd.Series({"NACCUDSD": 1, "NACCETPR": 88})
    assert dx.nacc_stage_etiology(row)[:2] == ("CN", "none")
    row = pd.Series({"NACCUDSD": 2, "NACCETPR": 99})
    assert dx.nacc_stage_etiology(row)[:2] == ("IMPNOMCI", "unknown")


# ----------------------------------------------------------------- raw UDS D1 (OASIS-3, ABVIB)
def test_uds_v3_lowercase_names_are_read():
    row = pd.Series({"DEMENTED": 1, "alzdis": 1, "alzdisif": 1, "epilep": 1, "epilepif": 3})
    stage, etio, epy = dx.uds_d1_stage_etiology(row)
    assert (stage, etio, epy) == ("DEM", "AD", True)             # epilepsy -> covariate only


def test_uds_contributing_etiology_not_primary():
    row = pd.Series({"DEMENTED": 1, "PROBAD": 1, "PROBADIF": 2, "CVD": 1, "CVDIF": 1})
    assert dx.uds_d1_stage_etiology(row)[:2] == ("DEM", "VaD")


def test_uds_impnomci_and_mci():
    assert dx.uds_d1_stage_etiology(pd.Series({"IMPNOMCI": 1}))[:2] == ("IMPNOMCI", "unknown")
    row = pd.Series({"MCIAMEM": 1, "PROBAD": 1, "PROBADIF": 1})
    assert dx.uds_d1_stage_etiology(row)[:2] == ("MCI", "AD")
    assert dx.uds_d1_stage_etiology(pd.Series({"NORMCOG": 1}))[:2] == ("CN", "none")


# ----------------------------------------------------------------- other cohorts
def test_adni_mci_due_to_ad_is_mci():
    assert adni._stage_adni(pd.Series({"DIAGNOSIS": 2, "DXMDUE": 1}))[:2] == ("MCI", "AD")
    assert adni._stage_adni(pd.Series({"DIAGNOSIS": 3, "DXAD": 1}))[:2] == ("DEM", "AD")
    assert adni._stage_adni(pd.Series({"DIAGNOSIS": 1}))[:2] == ("CN", "none")


def test_aibl_mcsa_nifd():
    assert aibl._stage_aibl(pd.Series({"DXCURREN": 3}))[:2] == ("DEM", "AD")
    assert aibl._stage_aibl(pd.Series({"DXCURREN": 2}))[:2] == ("MCI", "unknown")
    assert mcsa._stage_mcsa(pd.Series({"FCOGIMPR": 1}))[0] == "MCI"
    assert mcsa._stage_mcsa(pd.Series({"FCOGIMPR": 4}))[0] == "Unknown"   # undocumented code
    st, et, _ = nifd._stage_nifd(pd.Series({"DX": "BV"}))
    assert st == "Unknown" and et == "bvFTD"


def test_preventad_scd_is_cn():
    assert preventad._classify_clinical_dx("Subjective cognitive decline") == "CN"
    assert preventad._classify_clinical_dx("MCI") == "MCI"


# ----------------------------------------------------------------- trajectory
def _dx_frame(stages, etios, legacy, dates):
    d = pd.DataFrame({"DX_STAGE": stages, "DX_ETIOLOGY": etios, "DX_LEGACY": legacy,
                      "DATE": pd.to_datetime(dates)})
    d["DX_LABEL"] = [dx.harmonised_label(s, e) for s, e in zip(stages, etios)]
    return d


def test_trajectory_counts_stage_changes_only():
    g = _dx_frame(["CN", "CN", "MCI", "MCI", "DEM"], ["none", "none", "AD", "unknown", "AD"],
                  ["CN", "CN", "AD+MCI", "MCI", "AD"],
                  ["2010-01-01", "2011-01-01", "2012-01-01", "2013-01-01", "2014-01-01"])
    rec = dx.trajectory_record("P1", g, g["DATE"].tolist(), fmt=lambda d: str(d.date()))
    assert rec["1st_diag"] == "CN" and rec["Lst_diag"] == "AD"
    assert rec["nb_of_change"] == 2                      # CN->MCI, MCI->DEM (etiology change ignored)
    assert rec["1st_stage"] == "CN" and rec["Lst_stage"] == "DEM" and rec["max_stage"] == "DEM"
    assert rec["all_stage"].split("|") == ["CN", "CN", "MCI", "MCI", "DEM"]
    assert rec["1st_diag_legacy"] == "CN"


def test_trajectory_skips_unknown_stage():
    g = _dx_frame(["CN", "Unknown", "CN"], ["none", "unknown", "none"], ["CN", "x", "CN"],
                  ["2010-01-01", "2011-01-01", "2012-01-01"])
    rec = dx.trajectory_record("P2", g, g["DATE"].tolist(), fmt=lambda d: str(d.date()))
    assert rec["nb_of_change"] == 0 and rec["Lst_stage"] == "CN"


# ----------------------------------------------------------------- expand + anchor
def test_expand_stage_windows_and_conversion(tmp_path):
    from inbrain_dementia_data_analyser.integration.expand import expand_master
    from inbrain_dementia_data_analyser.cohorts.anchor import create_anchor_cohort
    master = pd.DataFrame({
        "uid": ["ADNI_1", "ADNI_2"], "Databank": ["ADNI", "ADNI"], "ID": ["1", "2"],
        "1st_diag": ["CN", "MCI"], "Lst_diag": ["AD", "CN"],
        "1st_diag_date": ["2010-01-01", "2010-01-01"],
        "all_diag": ["CN|MCI|AD", "MCI|CN"],
        "all_stage": ["CN|MCI|DEM", "MCI|CN"],
        "all_etiology": ["none|AD|AD", "AD|none"],
        "all_diag_date": ["2010-01-01|2011-01-01|2012-06-01", "2010-01-01|2011-01-01"],
        "Lst_stage": ["DEM", "CN"], "Lst_etiology": ["AD", "none"],
        "all_MRI_T1w_date": ["2011-01-20", "2010-01-10"],
        "all_MRI_T1w_ID": ["I1", "I2"],
    })
    mp, ep = tmp_path / "master.csv", tmp_path / "expanded.csv"
    master.to_csv(mp, index=False)
    expand_master(mp, ep)
    e = pd.read_csv(ep).set_index("uid")
    assert e.loc["ADNI_1", "stage_lt6"] == "CN" and e.loc["ADNI_1", "stage_lt12"] == "MCI"
    assert e.loc["ADNI_1", "conversion_to_stage"] == "MCI"
    assert e.loc["ADNI_1", "conversion_direction"] == "progression"
    assert e.loc["ADNI_2", "conversion_direction"] == "reversion"
    res = create_anchor_cohort(ep, anchor="MRI_T1w", required=(), out_dir=tmp_path)
    t = res.table.set_index("uid")
    assert t.loc["ADNI_1", "diag_at_anchor"] == "MCI"
    assert t.loc["ADNI_1", "stage_at_anchor"] == "MCI"
    assert t.loc["ADNI_1", "etiology_at_anchor"] == "AD"
    assert t.loc["ADNI_2", "stage_at_anchor"] == "MCI"
