"""v0.4.1: long plasma tables with a test map (ADNI FNIH Biomarkers Consortium,
A4/LEARN External Data). Synthetic values only."""
import re

import pandas as pd
import pytest
import yaml

from inbrain_dementia_data_analyser.datasets import a4
from inbrain_dementia_data_analyser.fluid import plasma_from_spec, plasma_join


def _long(rows):
    return pd.DataFrame(rows, columns=["PTID", "TEST", "VALUE", "_examday"])


def test_test_map_keeps_only_listed_tests_with_their_names():
    df = _long([("S1", "AlzPath_ptau217", 0.4, 0), ("S1", "C2N_ptau217", 1.1, 0),
                ("S2", "AlzPath_ptau217", 0.6, 30), ("S2", "Unlisted", 9.0, 30)])
    spec = dict(long=True, test_col="TEST", value_col="VALUE",
                test_map={"AlzPath_ptau217": "Plasma_pTau217_AlzPath",
                          "C2N_ptau217": "Plasma_pTau217_C2N"})
    out = plasma_join(plasma_from_spec(df, spec, "PTID", "_examday", None))
    names = {m.group(1) for c in out.columns if (m := re.match(r"all_(.+)_score$", c))}
    assert names == {"Plasma_pTau217_AlzPath", "Plasma_pTau217_C2N"}
    # two platforms of one analyte no longer collapse into one name (v0.4.0 crash)
    assert out.columns.is_unique


def test_adni_fnihbc_spec_has_unique_output_names():
    from pathlib import Path
    y = yaml.safe_load((Path(__file__).resolve().parents[1] / "configs/datasets/adni.yaml").read_text())
    specs = y["settings"]["plasma_specs"]
    fn = specs["FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES"]
    assert fn["long"] and fn["test_map"]
    assert len(set(fn["test_map"].values())) == len(fn["test_map"])


@pytest.mark.parametrize("raw,norm", [("003", "3"), ("3", "3"), ("000", "0"),
                                      (" BL ", "bl"), ("m24", "m24")])
def test_viscode_normalisation(raw, norm):
    assert a4._norm_viscode(pd.Series([raw])).iloc[0] == norm


def test_a4_plasma_specs_days_and_registry(tmp_path, cfg):
    ext = tmp_path / "Clinical" / "External Data"
    ext.mkdir(parents=True)
    # Lilly p-tau217: own collection day; '<LLOQ' has no number and is not counted
    pd.DataFrame({"BID": ["B1", "B1", "B2"], "VISCODE": [6, 66, 6],
                  "TESTCD": "PTAU217", "ORRES": ["0.20", "0.30", "<LLOQ"],
                  "METHOD": "MSD", "COLLECTION_DATE_DAYS_CONSENT": [70, 1800, 65]}
                 ).to_csv(ext / "biomarker_pTau217.csv", index=False)
    # Roche: test codes with trailing spaces, value in LABRESN
    pd.DataFrame({"BID": ["B1", "B1"], "VISCODE": ["001", "001"],
                  "LBTESTCD": ["GFAP ", "APOE4"], "LABRESN": ["0.09", "1"],
                  "LBMETHOD": "ECLIA", "LABD_DAYS_CONSENT": [0, 0]}
                 ).to_csv(ext / "biomarker_Plasma_Roche_Results.csv", index=False)
    # Araclon: zero-padded VISCODE, day from the registry
    pd.DataFrame({"BID": ["B2"], "VISCODE": ["003"], "LBTESTCD": ["TP42/TP40"],
                  "LBORRES": ["0.09"], "LBMETHOD": "ELISA"}
                 ).to_csv(ext / "biomarker_AB_Test.csv", index=False)
    registry = pd.DataFrame({"BID": ["B2", "B2"], "VISCODE": ["1", "3"], "EXAMDAY": [0, 50]})

    a4.PATHS = {}
    a4.A4_PLASMA_DIR = str(tmp_path)
    a4.A4_PLASMA_FILE_RX = re.compile(cfg.datasets["A4"].settings["plasma_file_regex"], re.I)
    a4.A4_PLASMA_SPECS = dict(cfg.datasets["A4"].settings["plasma_specs"])
    out = a4.build_plasma_a4(registry)

    days = lambda v: [float(t) for t in str(v).split("|")]
    assert days(out.loc["B1", "all_Plasma_pTau217_Lilly_date"]) == [70, 1800]
    assert pd.isna(out.loc["B2", "all_Plasma_pTau217_Lilly_score"])
    assert days(out.loc["B1", "all_Plasma_GFAP_Roche_score"]) == [0.09]
    assert not any("APOE" in c for c in out.columns)          # not in the test map
    assert days(out.loc["B2", "all_Plasma_AB42_40_Araclon_date"]) == [50]
