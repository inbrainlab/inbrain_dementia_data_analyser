# Configuration reference

All files are YAML and live in `configs/` (or in the folder given by
`config_dir=` / `--config-dir` / the environment variable `INBRAIN_CONFIG_DIR`).
They are validated when loaded: a problem stops the run with the file, the
entry and what is wrong. `inbrain-dementia check-config` runs only that check.

## `datasets/<name>.yaml`

| Key | Required | Meaning |
|---|---|---|
| `name` | yes | dataset key used everywhere (`ADNI`, `PREVENTAD`, `4RTNI`, ...); aliases such as `PREVENT-AD` are accepted in function calls |
| `label`, `portal`, `description` | | documentation |
| `folder` | yes | raw folder under the root path (`Tables_ADNI`) |
| `output_folder`, `output_file` | yes | where the per-dataset table is written, inside `folder` |
| `subject_id` | yes | subject identifier of the dataset |
| `time_axis` | yes | `date`, `days`, `months` or `visit` (documentation and generic tables) |
| `module` | yes | pipeline module in `src/.../datasets/` |
| `default_table` | | table used by variables that give no `table` |
| `allow_copy_prefix` | | accept the `Copy of ` prefix (default true) |
| `ida` | | `{enabled: bool, required: bool}` - LONI IDA search export |
| `tables` | yes | raw tables, see below |
| `settings` | | dataset-specific settings read by the pipeline |

Table entries:

| Key | Meaning |
|---|---|
| `pattern` | file name with `{date}` (any portal date format) and `*` |
| `required` | missing file = validation error (default true) |
| `recursive` | also search sub-folders |
| `id_columns` | at least one must be in the header |
| `date_columns` | at least one must be in the header |
| `columns` | all must be in the header |
| `sheet` | Excel sheet |
| `generic`, `id_column`, `date_column` | table read by the generic engine (no pipeline code); every variable pointing to it is built as a long table (one row per visit) |

Settings in use:

* ADNI: `plasma_specs`, `lab_specs`, `csf_extra_specs` (tables found by name
  stem under `Tables_ADNI/`, markers `[raw column, output name, (lab), (kit)]`),
  `medhist_tables`, `csf_main_lab`, `csf_main_kit`;
* ABVIB: `table_date_columns` (tables with their own date column);
* A4: `plasma_file_regex`;
* NACC: `chunk_size`.

## `variables.yaml`

See the header of the file and section 5 of the README. Entry keys: `label`,
`domain`, `kind` (`score`, `categorical`, `flag`), `unit`, `expand`,
`builtin`, `description`, `sources`. Source keys: `table`, `column`,
`date_column`, `transform`, `items`, `params`, `missing_values`,
`valid_range`, `range_inclusive`, `valid_codes`, `row_filter`,
`case_insensitive`, `instrument`, `suffix`, `date_suffix`, `note`.

Order of the cleaning rules: `valid_codes`, then `missing_values`, then
`valid_range`. Cleaning and transforms never modify the raw column; the
result goes to a helper column `__var__<name>`.

How each pipeline dates a series (not configurable per variable unless noted):

| Dataset | Date of a series |
|---|---|
| ADNI | `EXAMDATE`, else `VISDATE`; `date_column` overrides |
| A4 | registry `EXAMDAY` of the visit (days since consent) |
| ABVIB | UDS forms `VISITYR/VISITMO/VISITDAY`; `table_date_columns` for COGSCOR, LABORATORY, SUBCAT; `date_column` overrides |
| AIBL | `EXAMDATE` |
| MCSA | `VisitDate`; `date_column` overrides (imaging metrics: `ImagingDate`) |
| NACC | `VISITYR/VISITMO/VISITDAY` |
| NIFD | `CLINICAL_LINKDATE`; `date_column` overrides |
| OASIS | `days_to_visit` |
| PREVENT-AD | `Date_taken`, ordered by candidate age; `date_column` overrides |
| 4RTNI | `<instrument>_<tp>_<date_suffix>`, else the timepoint in months |

Transforms (`transforms.py`; add your own with `@register("name")`):

| Name | Result |
|---|---|
| `sum` | row sum of `items`; `params.min_count` answered items (default 1) |
| `mean` | row mean of `items`; `params.min_count` |
| `max` | row maximum (e.g. any of several 0/1 medication flags) |

`sum`, `mean` and `max` accept `params.item_missing_values`: item codes set to
missing before the items are combined (ECog: `[9]` = "don't know").
| `coalesce` | first non-missing of `items` in order |
| `difference` | `items[0] - items[1]` |
| `ratio` | `items[0] / items[1]` |
| `bmi` | from `params.weight`, `params.height` and fixed units (`weight_unit`, `height_unit`) or unit-code columns (`weight_unit_column`, `height_unit_column`: 1 = lb/in, 2 = kg/cm); implausible BMI (<10, >80) -> missing |
| `bmi_reported_or_computed` | `params.bmi_column` when 10-80, else `bmi` |
| `sum_timepoint_items` | 4RTNI only: sum of `<instrument>_<tp>_<item>` per timepoint |

## `harmonization.yaml`

`global_rename` (all datasets), `column_rename` (per dataset, exact names),
`family_rename` (per dataset, a whole measure family: every column that
contains the stem as a `_`-delimited token; a stem that already carries the
canonical prefix is left alone and a rename that would create an existing
name is skipped and reported), `age_priority` (first available column becomes
`age_at_baseline`).

## `demographics.yaml`

`sex_standard` and `sex_codes`, `race_codes` (per dataset) and `race_text`,
`hispanic_codes` and `hispanic_default`, `handedness_codes` and
`handedness_text`, `education_sentinels` and `education_max_years`,
`nacc_apoe_codes` and `nacc_old_apoe_fix` (correction for NACC tables written
with the old, wrong NACCAPOE map). Keys are raw codes as lower-case text;
`null` means "unknown" (missing). Codes not listed become missing and are
reported at the end of the merge.

## `diagnosis.yaml`

`diag_map` (upper-case raw label -> canonical label) and `label_columns`.

## `imaging.yaml`

`types` (name, modality, description, `expand`), `measure_groups`
(`PET_ABETA`, `PET_TAU`) and `fluid_group_patterns` (regular expressions on
stems: `PLASMA_PTAU217`, `CSF_AB42`, ...). Values of different platforms of a
fluid group are not on the same scale; the cohort tools report which stem was
matched.

## `integration.yaml`

`merge_order`, `output_folder`, `expand.*` (window size and horizon, channels,
carried subject-level columns, datasets with numeric time) and `time_axis.*`
(unit of numeric acquisition times per dataset for the cohort tools; datasets
whose image IDs carry the day offset).

## `cohorts/*.yaml`

Arguments of one cohort tool plus `tool:` (`anchor`, `cross_sectional`,
`multimodal`, `longitudinal`, `longitudinal_multimodal`); run with
`inbrain-dementia cohort <file> <root>`.
