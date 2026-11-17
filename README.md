<p align="center">
  <img src="docs/assets/inbrain_logo.jpg" alt="InBrain" width="320">
</p>

# inbrain_dementia_data_analyser

Harmonises ten public dementia and ageing cohorts into one table and builds
study cohorts from it: sets of subjects who have particular MRI / PET images,
cognitive scores and fluid biomarkers close together in time. A cohort comes
out as a table of subjects (or visits) with their clinical values and the
LONI image IDs to download.

Datasets: **ADNI, A4/LEARN, ABVIB, AIBL, MCSA, NACC, NIFD, OASIS-3,
PREVENT-AD, 4RTNI**.

```
raw tables (Tables_<DATASET>/)                                  configs/*.yaml
        │  validate_dataset()  - folders, files, columns               │
        ▼                                                              │
 preprocess_dataset()  ──►  Tables_<DS>/Integrated_<DS>/<DS>_longitudinal.csv
        │                   one row per subject, every measure as a dated list
        ▼
 integrate()  ──►  Integrated_ALL/MASTER_all_cohorts.csv        (merge + harmonisation)
              ──►  Integrated_ALL/MASTER_expanded_240m.csv      (6-month windows)
        ▼
 create_cross_sectional_cohort()      one row per subject
 create_multimodal_cohort()           2, 3, 4 ... modalities within a window
 create_longitudinal_cohort()         one modality at 0 / 12 / 24 ... months
 create_longitudinal_multimodal_cohort()  a schedule where every visit has a full set
```

The package is a refactor of the lab notebooks
`multicohort_explorator_plasma_v3_manufacturer.ipynb` and
`cohort_merge_integration_unified_v3_manufacturer.ipynb`. On the same input it
writes the same files (checked on the 2026-11-16 data snapshot; see
[docs/REFACTOR_REPORT.md](docs/REFACTOR_REPORT.md)).

> **No data in this repository.** All ten sources are distributed under data
> use agreements. Download them yourself from the portals listed in
> [docs/datasets.md](docs/datasets.md); never commit subject-level files (the
> `.gitignore` blocks CSV/XLSX files and the `Tables_*` / `Integrated_*` folders).

---

## 1. Installation

Python 3.10 or newer.

```bash
git clone https://github.com/vini-rodrigues-2003/inbrain_dementia_data_analyser.git
cd inbrain_dementia_data_analyser
pip install -e .            # editable: uses the configs/ folder of this clone
pip install -e ".[test]"    # optional, to run the tests: pytest
```

**Google Colab**

```python
from google.colab import drive
drive.mount("/content/drive")

!git clone -q https://github.com/vini-rodrigues-2003/inbrain_dementia_data_analyser.git /content/inbrain_dementia_data_analyser
!pip install -q /content/inbrain_dementia_data_analyser    # dependencies + package + default configs

import inbrain_dementia_data_analyser as ida
ida.enable_logging()
ROOT = "/content/drive/MyDrive/VRS_inbrain_documents/Tables_projects_dementia"
```

To edit the configuration in Colab, keep a copy on Drive
(`!inbrain-dementia init-config /content/drive/MyDrive/inbrain_configs`) and
load it with `ida.load_config("/content/drive/MyDrive/inbrain_configs")` (or
set `os.environ["INBRAIN_CONFIG_DIR"]`).

Memory: everything wide is read in row chunks; the whole pipeline (NACC
included) is meant to fit a standard Colab session. Pre-process the datasets
one at a time if RAM is tight.

## 2. Where the raw data go

One root folder, one sub-folder per dataset. The names are fixed (including
the historical spellings `Tables_4TRNI` and `Tables_PREVENT-AD`):

```
<root_path>/
├── Tables_ADNI/        All_Subjects_DXSUM_<date>.csv, PTDEMOG_<date>.csv, ...
├── Tables_A4/          A4_REGISTRY_PRV2_<date>.csv, ..., idaSearch_<m>_<d>_<yyyy>_A4.csv
├── Tables_ABVIB/       UDS_*_<date>.csv, COGSCOR_<date>.csv, ..., idaSearch_*.csv
├── Tables_AIBL/        aibl_<table>_01-Jun-2018.csv, idaSearch_*.csv
├── Tables_MCSA/        MCSA_Data_<date>.csv, idaSearch_*.csv
├── Tables_NACC/        investigator_ftldlbd_nacc<v>.csv, investigator_scan_*.csv, investigator_fcsf_*.csv
├── Tables_NIFD/        NIFD_Clinical_Data_<date>_updated.xlsx, idaSearch_*.csv
├── Tables_OASIS/       OASIS3_UDS*.csv, OASIS3_demographics.csv, OASIS3_MR_json.csv, ...
├── Tables_PREVENT-AD/  <Instrument>_Registered_PREVENTAD_<date>.csv
├── Tables_4TRNI/       4RTNI_DATA_<date>.csv, idaSearch_*.csv
└── Integrated_ALL/     written by integrate()
```

* The exact tables and file patterns of each dataset are listed in
  [docs/datasets.md](docs/datasets.md) (generated from `configs/datasets/`).
* Download dates in the names do not matter: `{date}` in a pattern accepts
  `08Sep2025`, `01-Jun-2018`, `2025-03-21`, `20250321`; the `Copy of ` prefix
  added by Google Drive is accepted; with several matches the newest wins.
* **LONI IDA search export** (A4, ABVIB, AIBL, MCSA, NIFD, 4RTNI): put the
  `idaSearch_*.csv` downloaded from the IDA image search in the dataset folder
  (any name starting with `idaSearch_` that has the columns `Study Date` and
  `Image ID`). Older per-visit exports can stay in `old_ida_search/`.
* Each pipeline writes its table to `Tables_<DS>/Integrated_<DS>/` (as the
  notebooks did), next to a `.provenance.json` record.

## 3. Quick start

```python
import inbrain_dementia_data_analyser as ida

ida.enable_logging()                      # progress messages
ROOT = "/path/to/Tables_projects_dementia"

# 1. validation: folders, files and columns
reports = ida.validate_all(ROOT)
for r in reports.values():
    print(r.summary())

# 2. pre-processing: one table per dataset
ida.preprocess_dataset("ADNI", ROOT)
ida.preprocess_all(ROOT, ["AIBL", "OASIS", "PREVENT-AD"])   # several, one after the other

# 3. integration: merged master + expanded master
paths = ida.integrate(ROOT)

# 4. cohorts
res = ida.create_cross_sectional_cohort(paths["expanded"], anchor="PET_ABETA",
                                        required=["PLASMA_PTAU217", "CDR"], margin_months=6)
print(res.n_subjects, res.files)
```

The same from the command line:

```bash
inbrain-dementia validate   $ROOT
inbrain-dementia preprocess $ROOT --datasets ADNI AIBL
inbrain-dementia integrate  $ROOT
inbrain-dementia cohort configs/cohorts/longitudinal_t1w.yaml $ROOT
```

[examples/quickstart.ipynb](examples/quickstart.ipynb) walks through the
whole flow.

## 4. Step by step

### 4.1 Validation

```python
rep = ida.validate_dataset(ROOT, "NIFD")
rep.ok            # False when a required table, a required column or the IDA export is missing
rep.errors        # what makes the output incomplete
rep.warnings      # optional tables, variable columns not found, several matching files
rep.to_frame()    # one row per table: status, matched file, missing columns
```

`preprocess_dataset` runs the validation first and stops on errors
(`ValidationError`) unless `allow_missing=True`. A notebook that could not
find a table printed `[SKIP]` and still wrote an output that looked complete;
the package refuses to do that silently.

### 4.2 Pre-processing

```python
df = ida.preprocess_dataset("A4", ROOT)          # returns the table and writes it
```

Output: one row per subject. Every scored measure `X` (test, scale, lab value,
fluid biomarker) becomes `all_X_date | all_X_score | 1st_X_* | Lst_X_* |
nb_of_X_visits`; every imaging type `T` becomes one token per session in
`all_T_date`, `T_IDs`, `all_T_image_ids`, `all_T_manufacturer`, ... See
[docs/output_schema.md](docs/output_schema.md).

### 4.3 Integration

```python
master = ida.merge_cohorts(ROOT)                  # MASTER_all_cohorts.csv (+ coverage report)
ida.expand_master(ROOT + "/Integrated_ALL/MASTER_all_cohorts.csv",
                  ROOT + "/Integrated_ALL/MASTER_expanded_240m.csv")
# or both: ida.integrate(ROOT)
```

The merge renames columns (`configs/harmonization.yaml`), recodes sex, race,
ethnicity, education, handedness and APOE per dataset
(`configs/demographics.yaml`; originals kept in `<column>_raw`), maps
diagnosis labels (`configs/diagnosis.yaml`) and stacks the datasets with
`Databank` and `uid = <Databank>_<ID>`. The expansion places every measure in
6-month windows around the first-diagnosis date (`pre…`, `lt6 … lt240`,
`gt240`).

### 4.4 Cohorts

All tools read the expanded master, accept exact stems (`MRI_T1w`, `CDR_SB`,
`Plasma_pTau217_C2N`) or groups (`PET_ABETA`, `PET_TAU`, `PLASMA_PTAU217`,
`PLASMA_NFL`, `CSF_AB42`, ... from `configs/imaging.yaml`), and return a
`CohortResult` (`table`, `images`, `summary`, `files`). Each writes the cohort
table, a long image list (`*_loniIDs.csv`) and the unique image IDs
(`*_loniIDs.txt`, comma-separated for the IDA search).

**Cross-sectional** - one row per subject:

```python
# anchor pairing, one event per subject
ida.create_cross_sectional_cohort(X, method="anchor", anchor="MRI_T1w",
                                  required=["CDR", "CDR_SB", "BMI"], optional=["MRI_T2w"],
                                  margin_months=6, margin_months_by_feature={"BMI": 12},
                                  anchor_events="most_optional")
# best multimodal set
ida.create_cross_sectional_cohort(X, method="multimodal", window_months=12, modalities=[
    {"label": "MRI_T1w", "prefixes": ["MRI_T1w"], "required": True},
    {"label": "PET_FDG", "prefixes": ["PET_FDG"], "required": True},
    {"label": "PET_ABETA", "prefixes": ["PET_ABETA"], "required": False}])
```

`create_anchor_cohort(..., anchor_events="all")` keeps every anchor event
(several rows per subject).

**Longitudinal** - one modality on a schedule:

```python
ida.create_longitudinal_cohort(X, modality_label="MRI_T1w", visit_targets_months=[0, 12, 24],
                               tolerance_months=6, min_visits=2)
```

**Longitudinal multimodal** - every visit must carry a full set:

```python
ida.create_longitudinal_multimodal_cohort(X, schedule_modality="MRI_T1w",
                                          required=["CDR", "PET_ABETA"], optional=["PET_FDG"],
                                          margin_months=12, visit_targets_months=[0, 24, 48],
                                          tolerance_months=9, min_visits=2)
```

A visit counts when the schedule modality has every required measure within
the margin; visits are then matched to the targets as in the longitudinal
tool. Output: one row per subject x matched visit. (This definition is new in
the package; see the report.)

Ready-made YAML specs for the command line are in `configs/cohorts/`.

## 5. Adding or changing a variable

Every scored variable lives in `configs/variables.yaml`. The entry name is the
output stem; each dataset lists where the variable is:

```yaml
variables:
  RAVLT_forgetting:
    label: RAVLT forgetting
    domain: memory
    expand: true                      # also place it in the time windows
    sources:
      ADNI:      {table: neurobat, column: "Ravlt foregetting"}
      PREVENTAD: {table: ravlt,    column: "Ravlt esq"}
```

After pre-processing both datasets write `all_RAVLT_forgetting_date`,
`all_RAVLT_forgetting_score`, ...; the merge puts them in the same columns and
the expansion creates `RAVLT_forgetting_<window>`.

What a source can say:

| Key | Meaning |
|---|---|
| `table` | key of a table in `configs/datasets/<dataset>.yaml` (omit for NIFD, MCSA, NACC, 4RTNI: one main table) |
| `column` | raw column, or a list of candidates (first present wins) |
| `date_column` | when the variable has its own date column (NIFD, MCSA, ADNI ECog) |
| `transform`, `items`, `params` | computed variables: `sum`, `mean`, `max`, `coalesce`, `difference`, `ratio`, `bmi`, `bmi_reported_or_computed`; `params.item_missing_values` drops item codes such as 9 = "don't know" |
| `missing_values` | codes meaning missing, e.g. `[-4, 88, 99]` |
| `valid_range`, `range_inclusive` | `[min, max]`, `null` = open |
| `valid_codes` | keep only these codes (UDS health history: `[0, 1, 2]`) |
| `row_filter` | `{column: PACCNAME, pattern: "composite"}` |
| `case_insensitive` | match the column name regardless of case |
| `instrument`, `suffix`, `date_suffix` | 4RTNI wide columns `<instrument>_<timepoint>_<suffix>` |

Examples:

```yaml
  # derived score: CDR sum of boxes recomputed from the six boxes (all six required)
  CDR_SB_boxes:
    sources:
      ABVIB:
        table: cdr
        transform: sum
        items: [MEMORY, ORIENT, JUDGMENT, COMMUN, HOMEHOBB, PERSCARE]
        params: {min_count: 6}

  # 4RTNI: MMSE_0_MMSETOT, MMSE_6MO_MMSETOT, MMSE_12MO_MMSETOT, dated by MMSE_<tp>_DATE
  MMSE:
    sources:
      4RTNI: {instrument: MMSE, suffix: MMSETOT, date_suffix: DATE}

  # NACC: UDS sentinels out
  ANIMALS_NACC:
    sources:
      NACC: {column: ANIMALS, valid_range: [0, 77]}
```

**A table the pipeline does not read yet.** Declare it in the dataset file
with `generic: true` and its ID / date columns; every variable pointing to it
is built without code:

```yaml
# configs/datasets/adni.yaml   (hypothetical table name)
tables:
  bnt_extra:
    pattern: All_Subjects_BNT_{date}.csv
    required: false
    generic: true
    id_column: PTID
    date_column: VISDATE
```

**A new transform.** Register a function and refer to it by name:

```python
from inbrain_dementia_data_analyser.transforms import register

@register("zscore")
def zscore(df, src):
    x = df[src.items[0]].astype(float)
    return (x - x.mean()) / x.std()
```

After editing:

```bash
inbrain-dementia check-config          # stops with the entry and the problem if something is wrong
inbrain-dementia list-variables --dataset ADNI
python tools/make_docs.py              # refresh docs/variables.md and docs/datasets.md
```

Column renames between datasets (for measures written under different stems by
the existing pipelines) are in `configs/harmonization.yaml`; demographic
codebooks in `configs/demographics.yaml`.

## 6. Configuration files

| File | Content |
|---|---|
| `configs/datasets/<name>.yaml` | folders, output file, subject ID, time axis, raw tables (pattern, required, ID/date columns), cohort settings (e.g. ADNI plasma/CSF table specs) |
| `configs/variables.yaml` | scored variables and their sources ([docs/variables.md](docs/variables.md)) |
| `configs/harmonization.yaml` | renames at merge time, baseline-age priority |
| `configs/demographics.yaml` | sex, race, ethnicity, handedness, education and APOE codebooks |
| `configs/diagnosis.yaml` | diagnosis label map |
| `configs/imaging.yaml` | imaging types, measure groups, fluid-biomarker groups |
| `configs/integration.yaml` | merge order, window size, time-axis units |
| `configs/cohorts/*.yaml` | example cohort definitions for the CLI |

Details: [docs/configuration.md](docs/configuration.md).

## 7. Reproducibility

* Every output gets `<file>.provenance.json`: package version, SHA-256 of the
  configuration, input files (name, size, modification time), parameters and
  output shape. No subject-level data.
* `python tools/compare_outputs.py <root> <reference_folder>` compares a run
  with saved reference outputs (for example the notebooks') value by value.
* `pytest` runs the tests (synthetic data only).

## 8. Known limitations

* **ADNI** - ECog totals (`ECOG_PT`, `ECOG_SP`) are sums of the answered
  items of the 12-item forms (code 9 "don't know" excluded); the usual ECog
  score is the mean, available with `transform: mean`. Plasma tables are found
  by table name anywhere under `Tables_ADNI/`. Ethnicity is not extracted.
* **A4** - no calendar dates: time is `EXAMDAY` (days since consent). The MRI
  day of LEARN subjects is imputed as the median day of the MRI visit (63).
  Visits whose VISCODE has no registry day get the placeholder 2970 days
  (99 x 30), except ADLPQ and PACC (missing).
* **AIBL** - clinical tables end in 2018; 601 subjects exist only in the IDA
  export (imaging, no diagnosis). 1.5 T T1w is a separate family `MRI_T1w_1.5T`.
* **MCSA** - dates are shifted per subject (first visit = 1950-01-01);
  intervals are real, calendar dates are not.
* **NACC** - the main file is read three times in chunks. The full NACC run was
  not repeated during the refactor (only a 30,000-row sample, identical to the
  notebook). Some SCAN images have no date and cannot be paired in time.
* **NIFD** - race is decoded with UDS codes (the file holds 50 / 99); the
  `MoCA` source is the column `MODREY` as in the notebook, which the NIFD
  dictionary describes as the modified Rey figure copy - check before use.
* **OASIS-3** - time is days since entry. The UDS A5 health-history table is
  optional. Centiloid = `Centiloid_fSUVR_rsf_TOT_CORTMEAN`.
* **PREVENT-AD** - dates at month precision; BMI uses the eligibility height.
* **4RTNI** - three timepoints (0, 6, 12 months); a missing date falls back to
  the month offset.
* `RAVLT_I` mixes trial 1 (ADNI, NACC) and the sum of trials 1-5 (PREVENT-AD),
  as in the notebooks.

## 9. Repository layout

```
configs/                YAML configuration (see section 6)
src/inbrain_dementia_data_analyser/
  config.py paths.py io.py validation.py provenance.py preprocess.py cli.py
  series.py variables.py transforms.py clinical.py fluid.py
  imaging/              ida.py (LONI IDA export), sessions.py, scanner.py
  datasets/             one pipeline per dataset (adni.py ... rtni4.py) + base.py
  integration/          merge.py, expand.py, common.py
  cohorts/              anchor.py, multimodal.py, longitudinal.py, longitudinal_multimodal.py
tests/                  pytest, synthetic data only
tools/                  compare_outputs.py (regression), make_docs.py
docs/                   datasets, variables, output schema, configuration, refactor report
examples/               quickstart notebook
```

## 10. Citing and data use

See `CITATION.cff`. Publications must also cite and acknowledge each dataset
as required by its data use agreement (ADNI, AIBL, NACC, OASIS-3, A4,
PREVENT-AD, MCSA, NIFD, 4RTNI, ABVIB). The licence of the code is still to be
chosen (see `LICENSE`).
