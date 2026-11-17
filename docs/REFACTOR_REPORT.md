# Refactor report (notebooks -> package)

Date: 2026-11-16. Source: `multicohort_explorator_plasma_v3_manufacturer.ipynb`
(32 cells) and `cohort_merge_integration_unified_v3_manufacturer.ipynb`
(13 cells), final versions.

## 1. Regression comparison: ran, outputs identical

The package was run on the lab's data snapshot of 2026-11-16 and every output
was compared value by value with the output of the notebooks on the same input
(same rows by key, same column set, same values; numbers compared to 1e-9).

| Output | Notebook | Package | Result |
|---|---|---|---|
| ADNI_longitudinal.csv | 5,029 x 439 | 5,029 x 439 | identical (before the ECog fix, section 3.7) |
| A4_longitudinal.csv | 6,945 x 391 | 6,945 x 391 | identical |
| ABVIB_longitudinal.csv | 311 x 366 | 311 x 366 | identical |
| AIBL_longitudinal.csv | 1,463 x 230 | 1,463 x 230 | identical |
| MCSA_longitudinal.csv | 6,226 x 427 | 6,226 x 427 | identical |
| NIFD_longitudinal.csv | 355 x 364 | 355 x 364 | identical |
| OASIS_longitudinal.csv | 1,378 x 443 | 1,378 x 443 | identical |
| PREVENTAD_longitudinal.csv | 348 x 533 | 348 x 533 | identical |
| 4RTNI_longitudinal.csv | 135 x 508 | 135 x 508 | identical |
| NACC_longitudinal.csv, first 30,000 rows of the main file | 17,064 x 467 | 17,064 x 467 | identical |
| MASTER_all_cohorts.csv (9 datasets) | 22,190 x 2,015 | 22,190 x 2,015 | identical |
| MASTER_column_coverage.csv | 2,016 x 10 | 2,016 x 10 | identical |
| MASTER_expanded_240m.csv | 22,190 x 7,185 | 22,190 x 7,185 | identical |
| Anchor pairing, T1w + CDR/CDR-SB/BMI, most_optional (pairs, summary, loniIDs csv/txt) | 9,339 subjects | 9,339 | identical |
| Anchor pairing, T1w + CDR (+T2w optional), all events | 23,071 events / 10,609 subjects | same | identical |
| Multimodal T1w + FDG, 12 months (table, loniIDs csv/txt) | 1,902 | 1,902 | identical |
| Longitudinal T1w 0/12/24 m, >= 2 visits (table, loniIDs csv/txt) | 3,784 | 3,784 | identical |

Run times on a 2.9 GB Linux VM (Python 3.10, pandas 2.3): ADNI 95 s, MCSA
45 s, A4 28 s, AIBL 14 s, OASIS 11 s, others < 7 s; merge 30 s; expansion
92 s (peak 1.6 GB); cohort tools 14-48 s each.

`pytest`: 47 tests on synthetic data, all passing (one per trap of the brief,
configuration loading and validation, transforms and cleaning rules, the
generic-table path, multimodal search against exhaustive search on 300 random
cases, visit matching).

To repeat the comparison on another machine: copy the notebook outputs to a
folder, run the package, then `python tools/compare_outputs.py <root> <folder>`
(see the docstring of the script).

## 2. Not verified

* **NACC end to end.** The full main file (720 MB, ~210,000 visits) does not
  fit the time limit of the sandbox; a 30,000-row sample was run through both
  the notebook cell and the package (identical). The full run should be done
  once in Colab; the per-run peak was 0.6 GB for the sample.
* **Google Colab.** Not executed there; installation was tested with a regular
  (non-editable) wheel install and with the source tree.
* **ADNI plasma specs.** The plasma tables were not in the snapshot; the code
  path is the notebook's, unchanged.
* **NACC subset utility** (notebook cell 31: MMSE/CDR/MoCA and T1w
  availability for a list of NACCIDs) was not ported; the same selection can
  be made from the merged master.

## 3. Intentional differences from the notebooks

1. **Configuration instead of edited cells.** Raw file names are patterns
   (any download date, `Copy of ` prefix), scored-variable lists are in
   `variables.yaml`, renames / codebooks / diagnosis map / imaging groups /
   expansion settings are YAML. The values are those of the notebooks.
2. **Validation before processing.** A missing required table, required
   column or (NIFD) IDA export raises `ValidationError` unless
   `allow_missing=True`. The notebooks printed `[SKIP]` and continued.
3. **Column order** inside a per-dataset table can differ from the notebook
   for scored measures (they now follow the order of `variables.yaml`). The
   column set and values are the same; the merged master was already sorted.
4. **Messages** go through `logging` and are in English. The anchor-pairing
   Excel sheets are named `summary`, `pairs`, `loniIDs` (were `resumo`,
   `pares`, `loniIDs`).
5. **Warnings** (pandas date inference, downcasting, chained assignment) are
   silenced only while a pipeline or tool runs; the notebooks silenced all
   warnings globally.
6. **Provenance records** (`<output>.provenance.json`) next to every output.
7. **ADNI ECog (fixed after the first comparison).** `ecog_pt` now reads
   `All_Subjects_ECOG12PT_*` (items `ECOGPT1..12`, date `COLDATE`); the
   notebook pointed to `All_Subjects_ECOGPT_*`, whose items are named
   `MEMORY1`, `LANG1`, ..., so `ECOG_PT` was never built. It now exists for
   1,021 subjects (1,119 visits). In `ECOG_PT` and `ECOG_SP` the item code 9
   ("don't know") is no longer added to the total: the notebook summed it,
   giving totals up to 74 on a 12-48 scale. `ECOG_SP` changed for 203 subjects
   (one subject whose items were all "don't know" lost its only visit); every
   other ADNI column is unchanged. The usual ECog score is the mean of the
   answered items (`transform: mean`); the total is kept as a sum to stay
   close to the notebook.
8. **New:** generic tables (`generic: true`), named transforms and cleaning
   rules for variables, `create_cross_sectional_cohort`,
   `create_longitudinal_multimodal_cohort`, a CLI.

## 4. Suspected bugs kept (behaviour of the notebooks reproduced)

| Where | What | One-line fix if the lab agrees |
|---|---|---|
| NIFD `MoCA` | Source column `MODREY`, which the NIFD/4RTNI batteries use for the modified Rey figure copy (4RTNI maps it to `Benson_copy`). | point `MoCA -> NIFD` to the NIFD MoCA total column, or remove it |
| A4 times | Visits whose VISCODE has no registry day get 99 x 30 = 2970 days for every measure, except ADLPQ and PACC, which get a missing time (the notebook attached EXAMDAY twice to those two tables). | treat unknown VISCODEs as missing in `attach_examday` (`datasets/a4.py`) and drop `_LEGACY_DOUBLE_ATTACH` |
| `RAVLT_I` | ADNI `AVTOT1` and NACC `REY1REC` are trial 1; PREVENT-AD `List_A_Sum_Trials_1_5` is the sum of trials 1-5. Pooled under one name. | split into `RAVLT_T1` and `RAVLT_T1_5` in `variables.yaml` |
| NACC vitals | Negative-code cleaning applies only to the columns listed in `NPS_NUMERIC_COLS` (`datasets/nacc.py`); a new NACC variable needs `valid_range` / `missing_values` in its entry. | - (documented) |

## 5. Decisions for the lab

* **Definition of the longitudinal multimodal cohort** (new): visits are the
  acquisitions of a schedule modality that have every required measure within
  the margin; complete visits closer than `session_window_days` are merged
  (most optional matches wins); targets are matched as in the longitudinal
  tool. Alternative definitions: require each measure only at some visits, or
  build the schedule from the clinical visits instead of an imaging modality.
* NIFD race coding (file uses 50 / 99, decoded as UDS; the dictionary says
  5 = multiracial, 13 subjects affected).
* AIBL: 601 imaging-only subjects kept without diagnosis; 1.5 T T1w kept apart
  from `MRI_T1w`.
* OASIS centiloid column `Centiloid_fSUVR_rsf_TOT_CORTMEAN`.
* Representative image of a session: lowest image ID, `Original` before
  `Pre-processed`.
* Marital status / employment / living situation are not harmonised across
  datasets yet (OASIS writes them as categorical series); ADNI ethnicity is
  not extracted; the diagnosis token `EPY` is not in `diagnosis.yaml`.
* Licence of the code and authorship (`LICENSE`, `CITATION.cff`).
* Outputs are written inside the raw-data folders (`Tables_<DS>/Integrated_<DS>/`),
  as the notebooks did. A separate output root is an easy addition if wanted.

## 6. Design notes

* **Port, not rewrite.** Each dataset module is its notebook cell with the
  configuration lifted out; the dataset-specific logic (diagnosis classifiers,
  description rules for imaging types, IDA visit mappings, fluid builders)
  stays in code. The shared helpers of the notebooks (IDA loader and
  classifier, session alignment, scanner channels, fluid series, BMI and UDS
  helpers, scored-series builder) exist once in the package.
* **Module state.** The dataset modules keep the notebooks' module-level
  names (`PATHS`, `OUT_LONG`, ...), filled by `configure(ctx)` at the start of
  `run(ctx)`. One dataset runs at a time; this is not thread-safe. Turning
  each module into a class is the natural next step once the lab is happy with
  the outputs.
* **Variables.** The dataset code asks `variables.yaml` which variables to
  build from each table and calls `materialise()`, which resolves the column,
  applies transforms and cleaning into a helper column, and leaves dating and
  ordering to the dataset (each source has its own visit conventions).
* **Memory.** Wide files are read in row chunks; the NACC main file only with
  the columns used by the configured variables; the expanded master is written
  in blocks.
