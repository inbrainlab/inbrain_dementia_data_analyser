import pytest
import yaml

from inbrain_dementia_data_analyser.config import ConfigError, canonical_dataset_name, load_config
from inbrain_dementia_data_analyser.transforms import TRANSFORMS


def test_default_config_loads(cfg):
    assert set(cfg.datasets) == {"ADNI", "A4", "ABVIB", "AIBL", "MCSA", "NACC", "NIFD", "OASIS",
                                 "PREVENTAD", "4RTNI"}
    assert len(cfg.variables) > 150
    for name, var in cfg.variables.items():
        for ds, srcs in var.sources.items():
            for s in srcs:
                table = s.table or cfg.datasets[ds].default_table
                assert table in cfg.datasets[ds].tables, (name, ds, table)
                assert s.transform is None or s.transform in TRANSFORMS


def test_expand_features_match_the_notebook_list(cfg):
    feats = cfg.expand_features()
    for stem in ("MMSE", "CDR", "CDR_SB", "BMI", "SEIZURES", "ABETA42", "Plasma_pTau217_simoa"):
        assert feats[stem] == "scored"
    for stem in ("MRI_T1w", "PET_FDG", "ABETA_PIB", "TAU_MK6240"):
        assert feats[stem] == "imaging"
    assert "MRI_T1w_1.5T" not in feats
    assert len(feats) == 64


def test_dataset_aliases(cfg):
    assert canonical_dataset_name("PREVENT-AD", cfg.datasets) == "PREVENTAD"
    assert canonical_dataset_name("4trni", cfg.datasets) == "4RTNI"
    with pytest.raises(ConfigError):
        canonical_dataset_name("UKBB", cfg.datasets)


def _edit_variables(config_dir, fn):
    p = config_dir / "variables.yaml"
    data = yaml.safe_load(p.read_text())
    fn(data["variables"])
    p.write_text(yaml.safe_dump(data, sort_keys=False))


@pytest.mark.parametrize("entry,message", [
    ({"sources": {"ADNI": {"table": "nope", "column": "X"}}}, "table 'nope' is not declared"),
    ({"sources": {"ADNI": {"table": "mmse", "transform": "nope", "items": ["A"]}}}, "unknown transform"),
    ({"sources": {"ADNI": {"table": "mmse"}}}, "give 'column'"),
    ({"sources": {"ADNI": {"table": "mmse", "colum": "MMSCORE"}}}, "unknown key"),
    ({"sources": {"UKBB": {"column": "X"}}}, "Unknown dataset"),
    ({"label": "no sources"}, "no sources given"),
])
def test_invalid_variable_entries_name_the_entry(config_copy, entry, message):
    _edit_variables(config_copy, lambda v: v.__setitem__("BROKEN_VAR", entry))
    with pytest.raises(ConfigError) as exc:
        load_config(config_copy, refresh=True)
    assert "BROKEN_VAR" in str(exc.value) or "UKBB" in str(exc.value)
    assert message in str(exc.value)


def test_new_variable_in_two_datasets(config_copy):
    """The README example: one variable, different column names per dataset."""
    def add(v):
        v["RAVLT_forgetting"] = {
            "label": "RAVLT forgetting", "domain": "memory", "expand": True,
            "sources": {"ADNI": {"table": "neurobat", "column": "Ravlt foregetting"},
                        "PREVENTAD": {"table": "ravlt", "column": "Ravlt esq"}}}
    _edit_variables(config_copy, add)
    cfg = load_config(config_copy, refresh=True)
    assert [n for n, _ in cfg.sources("ADNI", "neurobat")][-1] == "RAVLT_forgetting"
    assert cfg.sources("PREVENTAD", "ravlt")[-1][1].column == ("Ravlt esq",)
    assert cfg.expand_features()["RAVLT_forgetting"] == "scored"
