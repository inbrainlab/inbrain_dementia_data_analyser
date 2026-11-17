import os
import time

import pandas as pd

from inbrain_dementia_data_analyser.io import read_csv_lean
from inbrain_dementia_data_analyser.paths import pattern_to_regex, resolve_file


def test_pattern_dates_and_copy_prefix():
    rx = pattern_to_regex("All_Subjects_DXSUM_{date}.csv")
    assert rx.match("All_Subjects_DXSUM_08Sep2025.csv")
    assert rx.match("Copy of All_Subjects_DXSUM_08Sep2025.csv")
    assert rx.match("all_subjects_dxsum_2025-03-21.csv")
    assert not rx.match("All_Subjects_DXSUM.csv")
    assert pattern_to_regex("aibl_mmse_{date}.csv").match("aibl_mmse_01-Jun-2018.csv")


def test_pattern_does_not_swallow_longer_table_names():
    """TAUSUVR must not match TAUSUVR_PETSURFER (two different A4 tables)."""
    rx = pattern_to_regex("TAUSUVR_{date}.csv")
    assert rx.match("TAUSUVR_09Jun2026.csv")
    assert not rx.match("TAUSUVR_PETSURFER_09Jun2026.csv")


def test_resolve_picks_the_newest_download(tmp_path):
    for name in ("A4_VMRI_PRV2_26Jan2026.csv", "A4_VMRI_PRV2_09Jun2026.csv"):
        (tmp_path / name).write_text("BID\n1\n")
    m = resolve_file(tmp_path, "A4_VMRI_PRV2_{date}.csv")
    assert m.path.name == "A4_VMRI_PRV2_09Jun2026.csv"
    assert m.ambiguous


def test_resolve_prefers_original_over_drive_copy(tmp_path):
    (tmp_path / "Copy of X_22Mar2026.csv").write_text("a\n")
    time.sleep(0.01)
    (tmp_path / "X_22Mar2026.csv").write_text("a\n")
    os.utime(tmp_path / "Copy of X_22Mar2026.csv", None)
    assert resolve_file(tmp_path, "X_{date}.csv").path.name == "X_22Mar2026.csv"


def test_read_csv_lean_matches_full_read(tmp_path):
    """Trap 5: chunked reading of selected columns gives the same values."""
    df = pd.DataFrame({"uid": [f"A_{i}" for i in range(50)], "x": range(50),
                       "img_id": ["1321166", None] * 25, "y": ["a"] * 50})
    p = tmp_path / "wide.csv"
    df.to_csv(p, index=False)
    lean = read_csv_lean(p, ["uid", "x", "img_id"], text_cols=["img_id"], chunksize=7)
    assert list(lean["x"]) == list(range(50))
    assert lean["img_id"].dropna().unique().tolist() == ["1321166"]  # trap 4: never '1321166.0'
