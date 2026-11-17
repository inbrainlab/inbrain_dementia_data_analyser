import pandas as pd
import yaml

from inbrain_dementia_data_analyser.cli import main
from inbrain_dementia_data_analyser.config import load_config
from inbrain_dementia_data_analyser.variables import attach_generic_variables, iter_generic_series


def test_generic_table_needs_no_code(config_copy):
    """A new instrument: declare the table and the variable, no Python."""
    p = config_copy / "datasets" / "adni.yaml"
    d = yaml.safe_load(p.read_text())
    d["tables"]["ravlt_extra"] = {"pattern": "RAVLT_EXTRA_{date}.csv", "required": False, "generic": True,
                                  "id_column": "PTID", "date_column": "EXAMDATE"}
    p.write_text(yaml.safe_dump(d, sort_keys=False))
    v = config_copy / "variables.yaml"
    data = yaml.safe_load(v.read_text())
    data["variables"]["RAVLT_forgetting"] = {"sources": {"ADNI": {
        "table": "ravlt_extra", "transform": "difference", "items": ["T5", "DELAYED"]}}}
    v.write_text(yaml.safe_dump(data, sort_keys=False))
    cfg = load_config(config_copy, refresh=True)
    raw = pd.DataFrame({"PTID": ["S1", "S1", "S2"], "EXAMDATE": ["2020-01-01", "2021-01-01", "2020-06-01"],
                        "T5": [12, 10, 9], "DELAYED": [10, 5, 9]})
    part = next(iter_generic_series("ADNI", {"ravlt_extra": raw}, cfg))
    assert part.loc["S1", "all_RAVLT_forgetting_score"] == "2|5"
    master = pd.DataFrame({"ID": ["S1", "S3"], "x": [1, 2]})
    out = attach_generic_variables(master, "ADNI", {"ravlt_extra": raw}, cfg)
    assert out.loc[0, "nb_of_RAVLT_forgetting_visits"] == 2 and pd.isna(out.loc[1, "nb_of_RAVLT_forgetting_visits"])


def test_cli_check_config(capsys):
    assert main(["check-config"]) == 0
    assert "Configuration OK" in capsys.readouterr().out
