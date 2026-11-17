"""Cohort figures on a small synthetic cohort + expanded master."""
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("matplotlib")

from inbrain_dementia_data_analyser.figures import DxGrouper, load_figure_settings, make_cohort_figures  # noqa: E402
from inbrain_dementia_data_analyser.figures.data import (LOST, NO_VISIT, build_figure_data,  # noqa: E402
                                                         trajectories, visit_table)


def _synthetic(n=120, seed=0):
    rng = np.random.default_rng(seed)
    banks = rng.choice(["ADNI", "NACC", "OASIS", "NIFD"], n)
    first = rng.choice(["CN", "MCI", "AD", "bvFTD", "LBD", "AD+VaD", "Unknown", "PSP", "CON"], n,
                       p=[.35, .25, .15, .06, .05, .05, .04, .03, .02])
    rows, exp = [], []
    for i in range(n):
        uid = f"{banks[i]}_{i}"
        nv = int(rng.integers(1, 6))
        labs = [first[i]] * nv
        if first[i] == "CN" and nv > 2 and rng.random() < .3:
            labs[2:] = ["MCI"] * (nv - 2)
        if first[i] == "MCI" and nv > 1 and rng.random() < .4:
            labs[1:] = ["AD"] * (nv - 1)
        if banks[i] == "OASIS":                                     # day offsets
            dates = [str(365 * k) for k in range(nv)]
            base = "0"
        else:
            dates = [str((pd.Timestamp("2010-01-01") + pd.DateOffset(months=12 * k)).date()) for k in range(nv)]
            base = dates[0]
        has = {m: rng.random() < p for m, p in
               (("CDR", .9), ("MMSE", .7), ("MoCA", .4), ("PTAU", .3), ("ABETA42", .3))}
        rows.append({"uid": uid, "Databank": banks[i], "ID": str(i), "1st_diag": labs[0], "Lst_diag": labs[-1],
                     "diag_at_anchor": labs[0], "1st_diag_date": base, "Lst_diag_date": dates[-1],
                     "MRI_T1w_loniID": f"I{i}", "MRI_T1w_date": base,
                     "CDR": 0.5 if has["CDR"] else np.nan, "MMSE": 28 if has["MMSE"] else np.nan,
                     "MoCA": 25 if has["MoCA"] else np.nan})
        e = {"uid": uid, "Databank": banks[i], "ID": str(i), "1st_diag": labs[0], "all_diag": "|".join(labs),
             "all_diag_date": "|".join(dates), "1st_diag_date": base, "Lst_diag_date": dates[-1]}
        for k in range(nv):
            w = f"lt{12 * k + 6}"
            e[f"MRI_T1w_{w}"] = f"I{i}_{k}"; e[f"MRI_T1w_{w}_acq"] = dates[k]
            for m in ("CDR", "MMSE", "MoCA", "PTAU", "ABETA42"):
                if has[m] and (k == 0 or rng.random() < .6):
                    e[f"{m}_{w}"] = 1.0; e[f"{m}_{w}_acq"] = dates[k]
        exp.append(e)
    return pd.DataFrame(rows), pd.DataFrame(exp)


def test_grouping_rules(cfg):
    g = DxGrouper(load_figure_settings(config=cfg))
    assert g("CON") == "CN" and g("bvFTD") == "FTD / tauopathy" and g("LBD+PD") == "LBD"
    assert g("AD+VaD") == "Mixed" and g("Unknown") is None and g(np.nan) is None
    assert g.transition("CN", "MCI") == "Progression" and g.transition("AD", "MCI") == "Reversion"
    assert g.transition("AD", "LBD") == "Change of aetiology"


def test_trajectories_use_groups_and_time_axes(cfg):
    s = load_figure_settings(config=cfg)
    h = pd.DataFrame({"uid": ["A_1", "O_2"], "Databank": ["ADNI", "OASIS"],
                      "all_diag": ["CN|CN|MCI", "LBD|LBD+PD|AD"],
                      "all_diag_date": ["2010-01-01|2011-01-01|2012-01-01", "0|200|730"]})
    tr = trajectories(visit_table(h, DxGrouper(s)))
    assert tr.loc["A_1", "conversion_to_group"] == "MCI"
    assert abs(tr.loc["A_1", "months_to_conversion"] - 24.0) < 0.1
    assert tr.loc["O_2", "conversion_to_group"] == "AD"          # LBD -> LBD+PD is not a conversion
    assert abs(tr.loc["O_2", "months_to_conversion"] - 730 / 30.4375) < 0.1


def test_modalities_found_in_cohort_and_expanded(cfg):
    coh, exp = _synthetic()
    fd = build_figure_data(coh, settings=load_figure_settings(config=cfg), expanded=exp, config=cfg)
    assert fd.reference == "T1w" and fd.N == len(fd.sets["T1w"])
    assert "cohort column" in fd.modality_source["CDR"]
    assert "expanded master" in fd.modality_source["p-tau"]          # no PTAU column in the cohort
    assert fd.n_removed == int((coh["diag_at_anchor"] == "Unknown").sum())
    st = fd.states_by_bin()
    assert set(st[0].dropna()) <= set(fd.order)
    assert st.iloc[:, 1:].isin([LOST, NO_VISIT] + fd.order + ["Other / unknown"]).all().all()


def test_all_figures_are_written(cfg, tmp_path):
    coh, exp = _synthetic()
    rep = make_cohort_figures(coh, tmp_path, expanded=exp, config=cfg,
                              spec={"longitudinal": {"upset_min_n": 1, "sankey_min_flow": 1}})
    assert not rep.skipped, rep.skipped
    names = set(rep.files)
    assert {"modality_connectome_graph", "composite_cohort_lattice", "mmse_moca_venn",
            "dx_transition_sankey_6m", "followup_availability_area_by_dx",
            "converters_transition_upset"} <= names
    assert all(p[0].stat().st_size > 10_000 for p in rep.files.values())
    assert rep.tables_path.exists() and "mmse_moca_venn" in rep.summary_path.read_text()


def test_cohort_only_without_expanded(cfg, tmp_path):
    coh, _ = _synthetic()
    coh["converted"] = 0
    rep = make_cohort_figures(coh, tmp_path, config=cfg, figures=["dx_group_counts_log", "modality_venn",
                                                                  "non_converters_stability_upset",
                                                                  "dx_transition_sankey"])
    assert "dx_group_counts_log" in rep.files and "mmse_moca_venn" in rep.files
    assert "non_converters_stability_upset" in rep.files
    assert "dx_transition_sankey_6m" in rep.skipped               # needs the expanded master


def test_spec_overrides(cfg, tmp_path):
    coh, exp = _synthetic()
    spec = {"modalities": {"T1w": {"column": "MRI_T1w_loniID"}, "CDR": "CDR", "MMSE": "MMSE"},
            "bundles": {}, "venn": ["CDR", "MMSE"], "longitudinal": {"availability": ["CDR"]}}
    rep = make_cohort_figures(coh, tmp_path, expanded=exp, config=cfg, spec=spec,
                              figures=["composite_cohort_dx_stack", "modality_venn"])
    assert "cdr_mmse_venn" in rep.files
    tab = pd.read_excel(rep.tables_path, sheet_name="composite_cohort_dx_stack", index_col=0)
    assert list(tab.index) == ["T1w", "T1w + CDR", "T1w + MMSE", "T1w + CDR + MMSE"]


def test_greedy_lattice_with_many_steps(cfg, tmp_path):
    coh, exp = _synthetic()
    rep = make_cohort_figures(coh, tmp_path, expanded=exp, config=cfg, spec={"bundles": {}},
                              figures=["composite_cohort_lattice"])
    tab = pd.read_excel(rep.tables_path, sheet_name="composite_cohort_lattice")
    assert list(tab["n_modalities_added"].value_counts().sort_index()) == [1, 5, 4, 3, 2, 1]
