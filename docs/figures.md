# Cohort figures

`inbrain_dementia_data_analyser.figures` draws 14 figures that describe a
cohort built with the cohort tools. They work with every tool:

* cross-sectional cohorts: anchor (`*_summary.csv`, `*_pairs.csv`) and multimodal;
* longitudinal cohorts (`*_longitudinal.csv`) and longitudinal multimodal.

The cohort decides **who** is in the figures. The expanded master
(`MASTER_expanded_240m.csv`) gives **their history**: every diagnosis with its
date, and the 6-month windows of each measure.

```python
from inbrain_dementia_data_analyser.figures import make_cohort_figures
rep = make_cohort_figures("out/T1w_CDR_BMI_summary.csv", "out/figures",
                          expanded=".../Integrated_ALL/MASTER_expanded_240m.csv",
                          spec=None,            # YAML file or dict overriding configs/figures.yaml
                          figures="all")        # or "cross_sectional", "longitudinal", or a list of names
print(rep.files, rep.skipped)
```

```bash
inbrain-dementia figures out/T1w_CDR_BMI_summary.csv $ROOT --spec my_figures.yaml
inbrain-dementia figures out/T1w_CDR_BMI_summary.csv --expanded MASTER_expanded_240m.csv \
                 --figures dx_transition_sankey followup_availability_area_by_dx
```

Each run writes one PNG per figure (also PDF / SVG via `style.formats`),
`figure_source_tables.xlsx` with the numbers behind each figure, and
`FIGURES.md`. That file lists the figures, where each modality was found and
the diagnosis labels that were left out. matplotlib is required
(`pip install ".[figures]"`).

## The figures

| Name | Type | What it shows | Purpose |
|---|---|---|---|
| `modality_connectome_graph` | Network graph with a diagnosis ring (donut) around each node | Modalities as nodes (size = subjects with the modality), edges = subjects sharing two modalities (width / colour), ring = % of each diagnostic group | Which modalities coexist in the same subjects, and how the diagnostic profile changes between modalities |
| `modality_overlap_matrix` | Annotated lower-triangle heatmap | n and % (of the reference) of shared subjects per pair; diagonal = n per modality | Exact overlap values behind the graph |
| `composite_cohort_lattice` | Hasse diagram, diagnosis strip in each box | Composite cohorts that add one modality at a time, with n, % and diagnostic mix | Sample cost of requiring each additional modality |
| `composite_cohort_dx_stack` | 100 % stacked horizontal bars | Diagnosis in each composite cohort | Selection bias of each combination |
| `dx_group_counts_log` | Vertical bars, log scale | Subjects and % per diagnostic group | Group imbalance, rare groups included |
| `dx_subgroup_labels` | Horizontal bars coloured by group | Original labels inside groups that pool several labels | Transparency of the diagnostic grouping |
| `dx_by_databank_heatmap` | Contingency heatmap (log shading, totals row) | Diagnostic group x source databank | Which databanks contribute each diagnosis |
| `<a>_<b>_venn` (default `mmse_moca_venn`) | Two-set Venn, area-proportional | Subjects with only A, only B, or both (n, % of the union) | Overlap of two measures |
| `venn_region_by_databank` | 100 % stacked bars, grey scale | Source databank in each Venn region | Whether the use of A or B depends on each databank's protocol |
| `venn_region_by_dx` | 100 % stacked bars | Diagnosis in each Venn region | Whether the A and B groups differ clinically |
| `converters_transition_upset` | UpSet plot (bars + matrix with arrows + set bars) | Transitions of converters (from -> to) by type (progression, reversion, change of aetiology), median (IQR) months to conversion | Conversion trajectories and time to event |
| `non_converters_stability_upset` | UpSet plot (hatched bars + dot matrix) | Non-converters per diagnosis, median (IQR) follow-up, fraction with a single assessment | Diagnostic stability and real follow-up |
| `dx_transition_sankey_<w>m` | Longitudinal Sankey / alluvial | Diagnosis in each window since baseline (plus "No visit" and "Lost to follow-up"); flow width = subjects moving between consecutive windows | Stability, conversions, reversions and attrition over time |
| `followup_availability_area_by_dx` | Stacked area chart, one panel per measure | Subjects with each measure in each window since baseline, by baseline diagnosis | Attrition per measure and the feasible horizon for longitudinal analyses |

## Inputs and rules

**Subjects.** One row per subject of the cohort table (`uid`, or
`Databank` + `ID`). Longitudinal multimodal cohorts have one row per visit;
the first row of each subject is used.

**Diagnosis.** The cross-sectional figures use `diag_at_anchor` (anchor
cohorts), or else `1st_diag` (`dx_column=` to choose). Labels are grouped as
set in `configs/figures.yaml` (`diagnosis_groups`). By default: CN (incl.
NIFD `CON`), MCI, AD, LBD (incl. LBD+PD), FTD / tauopathy (bvFTD, FTD, svPPA,
nfvPPA, PSP, CBS), VaD, and Mixed (any other compound label). Labels without
a group (Unknown, Dementia-NOS, OtherDx, PD, ...) are left out. Their number
is printed in the figures and in `FIGURES.md`.

**Modalities.** Each modality is found, in this order:

1. `file:` - another cohort CSV (subject listed there);
2. `column:` - a column of the cohort table (presence = non-missing);
3. automatically: `<measure>`, `<measure>_value`, `<measure>_date`,
   `<measure>_loniID`, `match_<measure>_date`, `<measure>_V<k>_..._date`;
4. the expanded master, if the cohort table has none of these columns
   (`source: ever` = any acquisition, `source: baseline` = first window).
   `source:` also forces this lookup when a column exists.

`reference` is the modality every subject has (the anchor). Percentages are
"% of <reference>", and subjects without it are left out. `bundles` adds
several modalities as one step of the composite cohorts (default
`p-tau + Aβ42`). `composites: auto` draws every combination when there
are at most 4 steps. With more, each column adds every remaining step to the
largest composite of the previous column. A list of step combinations draws
exactly those (see `configs/figures/thesis_modality_files.yaml`).

**Time.** Diagnoses come from `all_diag` / `all_diag_date` of the expanded
master. Time is counted in months from the first diagnosis, with the same
rules as the cohort tools: calendar dates, or day offsets for OASIS and A4.
Times are compared only within a subject and the same time axis. A
**conversion** is the first later diagnosis in a different *group*. A change
inside one group (e.g. LBD -> LBD+PD) is not a conversion. The direction is
taken from `stage_rank`: CN < MCI < dementia; dementia -> dementia = change of
aetiology.

* Sankey: the state in each window is the last diagnosis of that window.
  "No visit" means no diagnosis in the window but one later. "Lost to
  follow-up" means no later diagnosis. "Other / unknown" is a label without
  a group.
* Follow-up availability uses the windows of the expanded master. They count
  from the first diagnosis, or from the first acquisition when a subject has
  no diagnosis date. The dashed line is the number of subjects still
  followed.
* Without the expanded master, the conversion figures use the cohort columns
  `1st_diag`, `conversion_to`, `converted`, `months_to_conversion`, and the
  Sankey and follow-up figures are skipped. Subjects of the cohort that are
  missing from the expanded master (e.g. an older master) also fall back to
  the cohort columns, and the figure note says how many.

## Spec examples

`configs/figures/` has three specs, each matched to a cohort spec in `configs/cohorts/`:

* `multimodal_t1_fdg.yaml` - PET / MRI modalities from the `match_*` columns,
  CDR from the expanded master, Venn of amyloid vs tau PET;
* `longitudinal_t1w.yaml` - clinical measures looked up in the first window,
  12-month windows up to 72 months;
* `thesis_modality_files.yaml` - one cohort file per modality and the ten
  composite cohorts of the thesis figure. On the thesis files the
  cross-sectional figures give the same numbers as the original `make_figs.py`.

A minimal spec:

```yaml
modalities:
  Amyloid PET: PET_ABETA          # short form = measure
  p-tau217: PLASMA_PTAU217
  CDR: CDR
reference: Amyloid PET
bundles: {}
venn: [p-tau217, CDR]
longitudinal: {window_months: 12, max_months: 96, availability: [Amyloid PET, CDR]}
style: {formats: [png, pdf]}
```
