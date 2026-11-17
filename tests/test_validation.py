import pandas as pd

from inbrain_dementia_data_analyser.validation import validate_dataset


def _write(p, cols, rows=1):
    p.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{c: 1 for c in cols}] * rows).to_csv(p, index=False)


def _nifd(root, cols=("LONI_ID", "CLINICAL_LINKDATE", "MMSE_TOT")):
    folder = root / "Tables_NIFD"
    folder.mkdir(parents=True)
    pd.DataFrame([{c: 1 for c in cols}]).to_excel(folder / "NIFD_Clinical_Data_20200724_updated.xlsx",
                                                  index=False)
    return folder


def test_missing_folder_is_an_error(tmp_path, cfg):
    rep = validate_dataset(tmp_path, "AIBL", config=cfg)
    assert not rep.ok and "folder not found" in rep.errors[0]


def test_missing_required_table_is_an_error(tmp_path, cfg):
    folder = tmp_path / "Tables_ABVIB"
    for name in ("UDS_D1_CLINICIAN_DIAGNOSIS_13Feb2026.csv", "UDS_B4_CDR_19Feb2026.csv"):
        _write(folder / name, ["PTID", "VISITYR"])
    rep = validate_dataset(tmp_path, "ABVIB", config=cfg)
    assert not rep.ok
    assert any("'gds'" in e for e in rep.errors)
    assert {t.key: t.status for t in rep.tables}["imaging"] == "optional-missing"


def test_missing_id_column_and_variable_columns(tmp_path, cfg):
    _nifd(tmp_path, cols=("SUBJECT", "CLINICAL_LINKDATE"))
    rep = validate_dataset(tmp_path, "NIFD", config=cfg)
    t = rep.tables[0]
    assert t.status == "missing-columns" and "LONI_ID" in t.missing_columns
    assert "MMSE_TOT" in t.missing_variable_columns        # warning, variable skipped
    assert any("IDA" in e for e in rep.errors)              # NIFD needs its IDA export


def test_complete_dataset_passes(tmp_path, cfg):
    folder = _nifd(tmp_path)
    pd.DataFrame({"Subject ID": ["1_S_1"], "Study Date": ["1/1/2010"], "Image ID": [1]}).to_csv(
        folder / "idaSearch_11_16_2026_NIFD.csv", index=False)
    rep = validate_dataset(tmp_path, "NIFD", config=cfg)
    assert rep.ok, rep.summary()
    assert rep.files()["clinical"].name.startswith("NIFD_Clinical_Data")
