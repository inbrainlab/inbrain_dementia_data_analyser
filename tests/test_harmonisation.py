import numpy as np
import pandas as pd

from inbrain_dementia_data_analyser.integration import merge
from inbrain_dementia_data_analyser.integration.common import configure, vrs_measure_stems, vrs_time


def test_race_code_5_depends_on_the_cohort(cfg):
    """Trap 9: codebooks are per cohort."""
    merge.configure(cfg)
    adni = merge.harmonise_demographics(pd.DataFrame({"race": [5, "1|5", 7]}), "ADNI")
    nacc = merge.harmonise_demographics(pd.DataFrame({"race": [5, 1, 99]}), "NACC")
    assert adni["race"].tolist()[:2] == ["White", "Multiracial"] and pd.isna(adni["race"].iloc[2])
    assert adni["race_detail"].iloc[1] == "AIAN+White"
    assert nacc["race"].tolist()[:2] == ["Asian", "White"] and pd.isna(nacc["race"].iloc[2])
    assert nacc["race_raw"].tolist() == [5, 1, 99]


def test_nacc_apoe_codes(cfg):
    merge.configure(cfg)
    df = pd.DataFrame({"APOE_genotype": [np.nan] * 4, "APOE_code": [1, 2, 6, 9]})
    out = merge.harmonise_demographics(df, "NACC")
    assert out["APOE_genotype"].tolist()[:3] == ["e3/e3", "e3/e4", "e2/e2"]
    assert pd.isna(out["APOE_genotype"].iloc[3])
    assert out["APOE_e4_count"].tolist()[:3] == [0, 1, 0]


def test_sex_and_education_sentinels(cfg):
    merge.configure(cfg)
    mcsa = merge.harmonise_demographics(pd.DataFrame({"gender": [1, 0], "education": [16, 99]}), "MCSA")
    assert mcsa["gender"].tolist() == ["M", "F"]
    assert mcsa["education"].iloc[0] == 16 and pd.isna(mcsa["education"].iloc[1])


def test_family_rename_never_double_prefixes(cfg):
    """Trap 8: 'T2w' -> 'MRI_T2w' must leave 'all_MRI_T2w_FLAIR_date' alone."""
    merge.configure(cfg)
    df = pd.DataFrame(columns=["all_T2w_date", "all_MRI_T2w_FLAIR_date", "nb_total_T1w"])
    out = merge.apply_family_rename(df, "ADNI")
    assert set(out.columns) == {"all_MRI_T2w_date", "all_MRI_T2w_FLAIR_date", "nb_total_MRI_T1w"}


def test_diagnosis_lists_are_mapped_token_by_token(cfg):
    merge.configure(cfg)
    assert merge.normalize_diag("CON|CON|MCI") == "CN|CN|MCI"
    assert merge.normalize_diag("AD+VAD") == "AD+VaD"


def test_exact_stems(cfg):
    """Trap 8: MRI_T2w is not MRI_T2w_FLAIR, CDR is not CDR_SB."""
    configure(cfg)
    header = [f"{s}_lt6{x}" for s in ("MRI_T2w", "MRI_T2w_FLAIR", "CDR", "CDR_SB") for x in ("", "_acq")]
    assert vrs_measure_stems("MRI_T2w", header) == ["MRI_T2w"]
    assert vrs_measure_stems("CDR", header) == ["CDR"]
    pet = ["ABETA_PIB_lt6", "ABETA_PIB_lt6_acq", "ABETA_AV45_pre6", "ABETA_AV45_pre6_acq"]
    assert set(vrs_measure_stems("PET_ABETA", pet)) == {"ABETA_PIB", "ABETA_AV45"}


def test_time_axis_kinds(cfg):
    configure(cfg)
    t, kind = vrs_time(pd.Series(["2010-01-01", "129", "3"]), pd.Series(["ADNI", "OASIS", "MCSA"]))
    assert kind.tolist()[:2] == ["date", "visit"] and t.iloc[1] == 129
    assert pd.isna(t.iloc[2])                       # MCSA visit numbers are not a time
