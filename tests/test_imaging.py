import numpy as np
import pandas as pd

from inbrain_dementia_data_analyser.datasets.a4 import normalise_a4_visit
from inbrain_dementia_data_analyser.datasets.rtni4 import classify_4rtni_description
from inbrain_dementia_data_analyser.imaging.ida import (ida_combine_blocks, ida_find, ida_load,
                                                        ida_match_index)
from inbrain_dementia_data_analyser.imaging.sessions import session_align_ids


def test_one_id_per_session_date():
    """Trap 1: dates and IDs are paired per session, never 'one date, six IDs'."""
    dates = ["2010-01-05", "2010-01-05", "2011-02-01", "2010-01-05"]
    ids = ["1321166", "1321170", "1400001", "1321168.0"]
    sess, rep, allids = session_align_ids(dates, ids)
    assert sess == ["2010-01-05", "2011-02-01"]
    assert rep == ["1321166", "1400001"]                 # lowest numeric ID of the date
    assert allids[0] == "1321166/1321170/1321168"         # trap 4: no '.0'


def test_ida_ids_cast_to_clinical_dtype():
    """Trap 2: text IDA IDs joined to integer clinical IDs duplicated every subject."""
    ida = pd.DataFrame({"nb_total_MRI_T1w": [1, 2]}, index=pd.Index(["311", "12"], name="ID"))
    out = ida_match_index(ida, pd.Series([12, 311, 400]))
    assert out.index.dtype.kind in "iu"
    clinical = pd.DataFrame({"CDR": [0.5, 0]}, index=pd.Index([12, 311], name="ID"))
    assert len(clinical.join(out, how="outer")) == 2


def test_blocks_are_taken_whole_from_one_source():
    """Trap 3: never mix the date list of one source with the ID list of another."""
    T = "MRI_T1w"
    p = pd.DataFrame({f"nb_total_{T}": [2, np.nan], f"all_{T}_date": ["2010-01-01|2011-01-01", np.nan],
                      f"{T}_IDs": [np.nan, np.nan]}, index=["A", "B"])
    s = pd.DataFrame({f"nb_total_{T}": [1, 1], f"all_{T}_date": ["2009-05-05", "2012-01-01"],
                      f"{T}_IDs": ["999", "777"]}, index=["A", "B"])
    out = ida_combine_blocks(p, s)
    assert out.loc["A", f"all_{T}_date"] == "2010-01-01|2011-01-01"
    assert pd.isna(out.loc["A", f"{T}_IDs"])             # not borrowed from the other source
    assert out.loc["B", f"{T}_IDs"] == "777"


def test_ida_export_found_by_content_and_sentinels_removed(tmp_path):
    """Trap 6 and 11: any idaSearch_* file with Study Date + Image ID; age 0 / sex X are missing."""
    pd.DataFrame({"Subject ID": ["002_S_1", "002_S_1"], "Sex": ["X", "F"], "Visit": ["Baseline"] * 2,
                  "Study Date": ["1/05/2010", "1/05/2010"], "Age": [0, 71], "Modality": ["MRI", "MRI"],
                  "Description": ["MPRAGE", "MPRAGE"], "Type": ["Original"] * 2,
                  "Imaging Protocol": ["Manufacturer=SIEMENS;Field Strength=3.0"] * 2,
                  "Image ID": [11, 12]}).to_csv(tmp_path / "idaSearch_1_2_2030_ANY.csv", index=False)
    pd.DataFrame({"Subject ID": [1], "Age": [70]}).to_csv(tmp_path / "idaSearch_old.csv", index=False)
    found = ida_find(str(tmp_path))
    assert [p.split("/")[-1] for p in found["new"]] == ["idaSearch_1_2_2030_ANY.csv"]
    ida = ida_load(str(tmp_path), legacy_fill=False)
    assert ida["age"].isna().sum() == 1 and ida["sex"].isna().sum() == 1


def test_a4_visit_matches_registry_viscode():
    """Trap 10: 'SCV4 (Screening MRI)' is registry VISCODE 4."""
    assert normalise_a4_visit("SCV4 (Screening MRI)") == "4"


def test_ge_sequence_names():
    """Trap 12: 4RTNI GE scans."""
    assert classify_4rtni_description("IR-SPGR") == "MRI_T1w"
    assert classify_4rtni_description("T2 FSE") == "MRI_T2w"
