# Output columns

## Per-dataset tables (`Tables_<DS>/Integrated_<DS>/<DS>_longitudinal.csv`)

One row per subject. Lists are text joined by `|` in chronological order;
`NA` stands for a missing token inside a list so that positions stay aligned.

**Scored measure `X`** (test, scale, vital sign, lab value, fluid biomarker):

| Column | Content |
|---|---|
| `all_X_date` | dates (or day offsets / months, depending on the dataset) |
| `all_X_score` | values, aligned with the dates |
| `1st_X_score`, `1st_X_date` | first non-missing value |
| `Lst_X_score`, `Lst_X_date` | last non-missing value |
| `nb_of_X_visits` | number of values |
| `nb_of_X_censored` | fluid markers: values reported as <LoQ / >ULoQ (kept at the limit) |
| `all_X_batch`, `all_X_rundate` | fluid markers: assay batch / run date, aligned |
| `X_lab`, `X_kit` | fluid markers: laboratory and assay kit |

Fluid names: `Plasma_<Analyte>_<Platform>`, `CSF_<Analyte>_<Platform>`,
`Lab_<Analyte>`; the main CSF panels keep `ABETA42`, `ABETA40`, `TAU`, `PTAU`.

**Imaging type `T`** (`MRI_T1w`, `MRI_T2w`, `MRI_T2w_FLAIR`, `MRI_T2star`,
`MRI_DTI`, `MRI_fMRI`, `MRI_ASL`, `PET_FDG`, `ABETA_AV45`, `ABETA_PIB`,
`ABETA_FBB`, `ABETA_Flutemetamol`, `ABETA_NAV4694`, `TAU_Flortaucipir`,
`TAU_MK6240`, `TAU_PI2620`; AIBL also `MRI_T1w_1.5T`). One token per
**session** (subject x type x date):

| Column | Content |
|---|---|
| `all_T_date` / `all_T_visit` | session date (or visit token) |
| `T_IDs` | one representative image ID per session (lowest ID; `Original` before `Pre-processed`) |
| `all_T_image_ids` | every image ID of each session, joined by `/` |
| `all_T_manufacturer`, `all_T_model`, `all_T_field_strength` | scanner per session |
| `all_T_suvr`, `all_T_centiloid` | PET values per session, when the source has them |
| `nb_total_T`, `nb_images_T` | sessions, image series |
| `T_manufacturer` | vendors seen for the subject |
| `1st_T_ID`, `Lst_T_ID`, `1st_T_visit`, `Lst_T_visit` | first / last session |

Invariant: for a subject and type every `all_T_*` list and `T_IDs` have
exactly `nb_total_T` tokens.

Summaries: `MRI_type_disponible`, `PET_tracer_disponible`, `IDA_legacy_types`,
`IDA_sex`, `Plasma_panels_disponible`, `Plasma_analytes_disponible`,
`CSF_panels_disponible`.

**Diagnosis trajectory:** `1st_diag`, `Lst_diag`, `1st_diag_date`,
`Lst_diag_date`, `all_diag`, `all_diag_date`, `nb_of_visits`, `nb_of_change`,
`1st_time_of_change`, `Lst_time_of_change`, `all_time_of_change`, `Diag_total`.
Since v0.4.0 (see `docs/diagnosis.md`): `1st_stage`, `Lst_stage`, `max_stage`,
`all_stage`, `1st_etiology`, `Lst_etiology`, `all_etiology`, `epilepsy_any`,
`nb_of_epilepsy_visits`, `nb_of_label_change` and the pre-0.4.0 labels
`1st_diag_legacy`, `Lst_diag_legacy`, `all_diag_legacy`. `all_diag` holds the
harmonised label; `all_stage` and `all_etiology` are aligned with it.

**Demographics:** `ID`, `gender`, `race`, `hispanic`, `education`,
`birth_date`, `age_at_1st_diag`, an age-at-entry column (name varies),
`APOE_genotype`, `APOE_e4_count`, plus dataset extras. These are the RAW codes
of each dataset; they are harmonised at the merge.

## Merged master (`Integrated_ALL/MASTER_all_cohorts.csv`)

All per-dataset columns after renaming, stacked; plus `Databank`,
`uid = <Databank>_<ID>`, `age_at_baseline`, and the harmonised demographics
(`gender` M/F, `race`, `race_detail`, `hispanic` 1/0, `education` in years,
`handedness` R/L/A, `APOE_genotype` e.g. `e3/e4`, `APOE_e4_count`,
`APOE_e4_any`), with the originals in `<column>_raw`. Diagnosis labels are
canonical (`CN`, `MCI`, `AD`, `Dementia-NOS`, `LBD`, `bvFTD`, `svPPA`,
`nfvPPA`, `PSP`, `CBS`, `VaD`, ...; compound labels joined by `+`).
`MASTER_column_coverage.csv`: non-missing values per column and dataset.

## Expanded master (`Integrated_ALL/MASTER_expanded_240m.csv`)

Windows relative to the first-diagnosis date (or the measure's own first
acquisition when the subject has no diagnosis date): `pre6 ... pre240`,
`pregt240` before; `lt6 ... lt240`, `gt240` after.

| Column | Content |
|---|---|
| `<stem>_<w>` | value (imaging: 1 = present); a Python-list string when the window holds several acquisitions |
| `<stem>_<w>_acq` | date, day offset or visit token of the acquisition |
| `<stem>_<w>_id` | representative image ID (imaging) |
| `<stem>_<w>_image_ids`, `_manufacturer`, `_model`, `_suvr`, `_centiloid` | imaging channels |
| `<stem>_<w>_method`, `_batch`, `_rundate` | fluid channels |
| `diag_<w>`, `diag_<w>_acq` | diagnosis (harmonised label) in the window |
| `stage_<w>`, `etiology_<w>` | stage and primary etiology in the window (same visits as `diag_<w>`) |
| `months_to_conversion`, `conversion_to`, `converted`, `conversion_to_stage`, `conversion_direction` | first change of STAGE (progression / reversion); `conversion_to` = label at that visit |
| `baseline_anchor` | `diagnosis` or `first_acq` |

## Cohort tools

**Anchor pairing** (`<prefix>_pairs.csv`, one row per anchor event):
`uid`, kept subject columns, `anchor_event_n`, `anchor_kind`,
`<anchor>_window/_value/_acq/_id/...`, and per matched measure `M`:
`M_value`, `M_acq`, `M_id`, `M_window`, `M_delta_days` (+ `M_stem` for groups,
channels, `M_kit` / `M_lab` for fluids); `n_optional_matched`,
`diag_at_anchor`, `stage_at_anchor`, `etiology_at_anchor`, `diag_source`, `months_from_1st_diag`
(+ `_source` exact/window), `age_at_anchor` (+ `_source`).
`<prefix>_summary.csv`: the same, compact (loniID or value + date per measure).
`<prefix>_loniIDs.csv`: `uid`, `Databank`, `ID`, `anchor_event_n`,
`anchor_acq`, `feature`, `stem`, `loniID`, `acq`, `delta_days`, `required`.

**Multimodal** (`<prefix>_multimodal.csv`, one row per subject): per modality
`L`: `match_L_date`, `match_L_id`, `match_L_image_ids`, `match_L_value`,
`match_L_tracer`, `match_L_manufacturer`, `match_L_model`, `match_L_suvr`,
`match_L_centiloid`, `match_L_offset_days`; `match_span_days`,
`match_n_modalities`, `match_modalities`, `match_time_axis`.

**Longitudinal** (`<prefix>_longitudinal.csv`, one row per subject): per visit
`<mod>_V<n>_<t>m_` + `date`, `id`, `all_ids`, `n_series`, `offset_days`,
`manufacturer`, `model`, `matched`; `<mod>_baseline_date`, `_n_matched`,
`_n_sessions`, `_n_series_total`, `_n_sessions_unused`, `_time_axis`,
`_scanner_change`.

**Longitudinal multimodal** (`<prefix>_longitudinal_multimodal.csv`, one row
per subject x matched visit): `visit` (`V1_0m` ...), `target_months`,
`offset_days`, then every column of the anchor event of that visit (as in the
anchor pairing). `<prefix>_longitudinal_multimodal_subjects.csv`: one row per
subject (`n_visits_matched`, `n_complete_sessions`, `followup_days`, `visits`).
