# Diagnosis model (v0.4.0)

Every clinical visit is described on two separate axes plus a covariate:

| Axis | Values | Meaning |
| :--- | :--- | :--- |
| **stage** | `CN`, `IMPNOMCI`, `MCI`, `DEM`, `Unknown` | syndromic severity at the visit |
| **etiology** | `AD`, `LBD`, `FTD`, `VaD`, `PSP`, `CBS`, `MSA`, `HD`, `PD`, `Other`, `Epilepsy`, syndromes of NIFD (`bvFTD`, `nfvPPA`, `svPPA`); `none` (CN); `unknown` (impaired, no etiology) | **primary** presumed etiology; several primaries are joined with `+` |
| **epilepsy** | `True` / `False` | epilepsy recorded at the visit (covariate, never part of the label) |

The **harmonised label** (`all_diag`, `1st_diag`, `Lst_diag`, `diag_<w>`,
`diag_at_anchor`) is derived from the two axes:

| stage | label |
| :--- | :--- |
| `CN`, `IMPNOMCI`, `MCI` | the stage (MCI due to AD is `MCI`; its etiology stays in `etiology`) |
| `DEM` | the etiology (`AD`, `LBD`, `AD+VaD` …); `Dementia-NOS` without etiology; `Dementia-Other` for `Other` / `Epilepsy` |
| `Unknown` | the etiology / syndrome (NIFD, 4RTNI, ADNI PD-only); `OtherDx` for `Other`; `Unknown` otherwise |

Trajectories are computed on the **stage** only: `nb_of_change`,
`*_time_of_change` and the conversions (`converted`, `months_to_conversion`,
`conversion_to_stage`, `conversion_direction` = progression / reversion) ignore
visits with an `Unknown` stage, a new etiology and epilepsy. `nb_of_label_change`
counts the changes of the harmonised label. The label used before v0.4.0 is
kept in `1st_diag_legacy`, `Lst_diag_legacy` and `all_diag_legacy`.

## Rules per cohort

| Cohort | Stage | Primary etiology |
| :--- | :--- | :--- |
| NACC | `NACCUDSD` (1 CN, 2 IMPNOMCI, 3 MCI, 4 DEM) | `NACCETPR` (1 AD, 2 LBD, 3 MSA, 4 PSP, 5 CBS, 6–7 FTD, 8 VaD, 11 HD, 15 Epilepsy, others `Other`; 88 / 99 → `none` / `unknown`); raw UDS D1 when both are missing |
| OASIS-3, ABVIB | UDS D1: `NORMCOG` → CN, `DEMENTED` → DEM, any MCI subtype → MCI, `IMPNOMCI` → IMPNOMCI | UDS D1 flags with `…IF = 1` (primary); a flag without its `…IF` column (UDS v1–v2 forms) also counts; contributing (`IF = 2`) and non-contributing (`IF = 3`) never enter the label. UDS v3 names (`alzdis`, `lbdis`, `cvd`, `ftldmo`, `ftldnos`, `msa`, `epilep` …) are read in upper or lower case |
| ADNI | `DIAGNOSIS` (DXSUM, ADNI 2/3/4), else `DXCURREN`, else the legacy flags | MCI: `DXMDUE` (1 AD, 2 Other); dementia: `DXAD` / `DXAPP` → AD, `DXODES` (1, 4, 5 FTD; 9 VaD), `DXPCOG = 3` → LBD, `DXOTHDEM` → Other; Parkinson's only → stage Unknown, etiology PD |
| AIBL | `DXCURREN` (1 CN, 2 MCI, 3 DEM-AD, 5 DEM-Other; 4 PD → Unknown/PD) | from `DXCURREN` |
| MCSA | `FCOGIMPR` (0 CN, 1 MCI, 3 DEM; code 4 is undocumented → Unknown) | not recorded (`unknown` when impaired) |
| PREVENT-AD | `RC_MCI` (research MCI), clinical diagnosis text (SCD / subjective → CN, MCI, Alzheimer → AD, dementia) | from the clinical text |
| NIFD | controls (`CON`) → CN; patients → Unknown (no severity stage) | syndrome (`bvFTD`, `nfvPPA`, `svPPA`, `Other`) |
| 4RTNI | Unknown (static) | `PSP`, `CBS`, `Other` |
| A4 | CN (static, cognitively unimpaired at entry) | `none` |

## Practical consequences

* Cohorts selected on `diag_at_anchor = MCI` now include MCI due to AD; use
  `etiology_at_anchor` to split them.
* `IMPNOMCI` (impaired, not MCI) is a stage of its own and is not grouped with
  MCI in the default figures.
* A CN → AD "conversion" now needs a change of stage (CN → DEM), not only a new
  AD flag at an MCI visit.
* Per-dataset tables written by an older version still integrate: the
  expansion falls back to the label when `all_stage` is missing.
