# Datasets and raw files

Generated from `configs/datasets/*.yaml` by `tools/make_docs.py`; do not edit by hand.

File patterns: `{date}` = any download date (`08Sep2025`, `01-Jun-2018`, `2025-03-21`, `20250321`); `*` = any text. Matching ignores case and accepts the `Copy of ` prefix added by Google Drive. When several files match, the newest date in the name wins.

## A4 - Anti-Amyloid Treatment in Asymptomatic Alzheimer's (A4) / LEARN

* Portal: https://www.a4studydata.org (and LONI IDA image search)
* Raw folder: `<root>/Tables_A4/`
* Output: `<root>/Tables_A4/Integrated_A4/A4_longitudinal.csv`
* Subject ID: `BID`; time axis: days
* LONI IDA search export: `idaSearch_*.csv` in the raw folder (optional; the dataset tables are the fallback)

| Table | File pattern | Required | ID column | Date column | Notes |
|---|---|---|---|---|---|
| `registry` | `A4_REGISTRY_PRV2_{date}.csv` | yes | BID | EXAMDAY | Visit registry (EXAMDAY = days since consent) |
| `ptdemog` | `A4_PTDEMOG_PRV2_{date}.csv` | yes | BID |  |  |
| `subjinfo` | `A4_SUBJINFO_PRV2_{date}.csv` | yes | BID |  |  |
| `cdr` | `A4_CDR_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `mmse` | `A4_MMSE_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `cogdigit` | `A4_COGDIGIT_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `coglogic` | `A4_COGLOGIC_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `cogfcsr16` | `A4_COGFCSR16_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `psychwell` | `A4_PSYCHWELL_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `cfi` | `A4_CFI_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `cfisp` | `A4_CFISP_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `adlpq` | `A4_ADLPQ_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `vitals` | `A4_VITALS_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `habits` | `A4_HABITS_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `ies` | `A4_IES_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `ftpscale` | `A4_FTPSCALE_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `concerns` | `A4_CONCERNS_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `sppacc` | `A4_SPPACC_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `inithealth` | `A4_INITHEALTH_PRV2_{date}.csv` | yes | BID |  |  |
| `petsuvr` | `A4_PETSUVR_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `vmri` | `A4_VMRI_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `tausuvr` | `TAUSUVR_{date}.csv` | yes | ID / BID |  |  |
| `tausuvr_ps` | `TAUSUVR_PETSURFER_{date}.csv` | yes | ID / BID |  |  |
| `dwi_params` | `A4_DWI_PARAMS_PRV2_{date}.csv` | yes | BID | VISCODE |  |
| `fmri_params` | `A4_FMRI_PARAMS_PRV2_{date}.csv` | yes | BID | VISCODE |  |

## ABVIB - Aging Brain: Vasculature, Ischemia, and Behavior

* Portal: LONI IDA - https://ida.loni.usc.edu
* Raw folder: `<root>/Tables_ABVIB/`
* Output: `<root>/Tables_ABVIB/Integrated_ABVIB/ABVIB_longitudinal.csv`
* Subject ID: `PTID`; time axis: date
* LONI IDA search export: `idaSearch_*.csv` in the raw folder (optional; the dataset tables are the fallback)

| Table | File pattern | Required | ID column | Date column | Notes |
|---|---|---|---|---|---|
| `dx` | `UDS_D1_CLINICIAN_DIAGNOSIS_{date}.csv` | yes | PTID / SUBJECT_ID | VISITYR |  |
| `cdr` | `UDS_B4_CDR_{date}.csv` | yes | PTID / SUBJECT_ID | VISITYR |  |
| `gds` | `UDS_B6_GDS_{date}.csv` | yes | PTID / SUBJECT_ID | VISITYR |  |
| `npsych` | `UDS_C1_MMSE_NPSYCH_BATTERY_{date}.csv` | yes | PTID / SUBJECT_ID | VISITYR |  |
| `cogscor` | `COGSCOR_{date}*.csv` | yes | PTID / SUBJECT_ID | TestDate |  |
| `demog` | `UDS_A1_SUBJECT_DEMOGRAPHICS_{date}.csv` | yes | PTID / SUBJECT_ID |  |  |
| `subcat` | `SUBCAT_{date}.csv` | yes | PTID / SUBJECT_ID | DATE |  |
| `apoe` | `IVD_APOE_{date}.csv` | yes | PTID / SUBJECT_ID |  |  |
| `lab` | `LABORATORY_{date}.csv` | yes | PTID / SUBJECT_ID | CDATE |  |
| `imaging` | `UDS_E1_IMAGING_LABS_{date}.csv` | no | PTID / SUBJECT_ID | VISITYR | UDS E1 imaging flags (fallback to IDA) |
| `meds` | `UDS_A4_SUBJECT_MEDICINE_{date}.csv` | yes | PTID / SUBJECT_ID | VISITYR |  |

## ADNI - Alzheimer's Disease Neuroimaging Initiative

* Portal: LONI IDA - https://ida.loni.usc.edu (study data: Download > Study Files)
* Raw folder: `<root>/Tables_ADNI/`
* Output: `<root>/Tables_ADNI/Integrated_ADNI/ADNI_longitudinal.csv`
* Subject ID: `PTID`; time axis: date


| Table | File pattern | Required | ID column | Date column | Notes |
|---|---|---|---|---|---|
| `dx` | `All_Subjects_DXSUM_{date}.csv` | yes | PTID / RID | EXAMDATE | Diagnostic summary |
| `demog` | `PTDEMOG_{date}.csv` | yes | PTID / RID |  | Participant demographics |
| `study_entry` | `All_Subjects_Study_Entry_{date}.csv` | yes | subject_id |  | Study entry (age at entry) |
| `apoe` | `All_Subjects_APOERES_{date}.csv` | yes | PTID / RID |  | APOE genotype |
| `mri` | `All_Subjects_Key_MRI_{date}.csv` | yes | subject_id | image_date | Key MRI image list |
| `pet` | `All_Subjects_PET_Images_{date}.csv` | yes | subject_id | pet_date | PET image list |
| `neurobat` | `All_Subjects_NEUROBAT_{date}.csv` | yes | PTID / RID | EXAMDATE / VISDATE | Neuropsychological battery |
| `mmse` | `All_Subjects_MMSE_{date}.csv` | yes | PTID / RID | EXAMDATE / VISDATE |  |
| `moca` | `All_Subjects_MOCA_{date}.csv` | yes | PTID / RID | EXAMDATE / VISDATE |  |
| `cdr` | `All_Subjects_CDR_{date}.csv` | yes | PTID / RID | EXAMDATE / VISDATE |  |
| `faq` | `All_Subjects_FAQ_{date}.csv` | yes | PTID / RID | EXAMDATE / VISDATE |  |
| `npiq` | `All_Subjects_NPIQ_{date}.csv` | yes | PTID / RID | EXAMDATE / VISDATE |  |
| `gds` | `All_Subjects_GDSCALE_{date}.csv` | yes | PTID / RID | EXAMDATE / VISDATE |  |
| `stai` | `All_Subjects_STAIAD_{date}.csv` | yes | PTID / RID | EXAMDATE / VISDATE |  |
| `ecog_pt` | `All_Subjects_ECOG12PT_{date}.csv` | yes | PTID / RID | COLDATE | Everyday Cognition 12-item, participant (items ECOGPT1-12) |
| `ecog_sp` | `All_Subjects_ECOG12SP_{date}.csv` | yes | PTID / RID |  | Everyday Cognition 12-item, study partner |
| `vitals` | `VITALS_{date}.csv` | yes | PTID / RID | EXAMDATE / VISDATE |  |
| `csf` | `All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_{date}.csv` | yes | PTID / RID | EXAMDATE | CSF Roche Elecsys (UPenn) |

## AIBL - Australian Imaging, Biomarker and Lifestyle study

* Portal: LONI IDA - https://ida.loni.usc.edu (AIBL study data + image collection search)
* Raw folder: `<root>/Tables_AIBL/`
* Output: `<root>/Tables_AIBL/Integrated_AIBL/AIBL_longitudinal.csv`
* Subject ID: `RID`; time axis: date
* LONI IDA search export: `idaSearch_*.csv` in the raw folder (optional; the dataset tables are the fallback)

| Table | File pattern | Required | ID column | Date column | Notes |
|---|---|---|---|---|---|
| `dx` | `aibl_pdxconv_{date}.csv` | yes | RID |  | Diagnosis per visit (no exam date) |
| `dx_date_src` | `aibl_mmse_{date}.csv` | yes | RID | EXAMDATE | Exam dates borrowed for the diagnosis |
| `demog` | `aibl_ptdemog_{date}.csv` | yes | RID |  |  |
| `apoe` | `aibl_apoeres_{date}.csv` | yes | RID |  |  |
| `mri_15t` | `aibl_mrimeta_{date}.csv` | no | RID | EXAMDATE | 1.5 T MRI metadata (fallback to IDA) |
| `mri_3t` | `aibl_mri3meta_{date}.csv` | no | RID | EXAMDATE | 3 T MRI metadata (fallback to IDA) |
| `pet_pib` | `aibl_pibmeta_{date}.csv` | no | RID | EXAMDATE |  |
| `pet_av45` | `aibl_av45meta_{date}.csv` | no | RID | EXAMDATE |  |
| `pet_flute` | `aibl_flutemeta_{date}.csv` | no | RID | EXAMDATE |  |
| `pet_nav` | `aibl_navmeta_{date}.csv` | no | RID | EXAMDATE |  |
| `neurobat` | `aibl_neurobat_{date}.csv` | yes | RID | EXAMDATE |  |
| `mmse` | `aibl_mmse_{date}.csv` | yes | RID | EXAMDATE |  |
| `cdr` | `aibl_cdr_{date}.csv` | yes | RID | EXAMDATE |  |
| `labdata` | `aibl_labdata_{date}.csv` | yes | RID |  |  |
| `medhist` | `aibl_medhist_{date}.csv` | yes | RID |  |  |

## MCSA - Mayo Clinic Study of Aging

* Portal: LONI IDA - https://ida.loni.usc.edu
* Raw folder: `<root>/Tables_MCSA/`
* Output: `<root>/Tables_MCSA/Integrated_MCSA/MCSA_longitudinal.csv`
* Subject ID: `mcsa_id`; time axis: date
* LONI IDA search export: `idaSearch_*.csv` in the raw folder (optional; the dataset tables are the fallback)

| Table | File pattern | Required | ID column | Date column | Notes |
|---|---|---|---|---|---|
| `data` | `MCSA_Data_{date}.csv` | yes | mcsa_id / MCSA_ID | VisitDate | One wide file, one row per visit (dates shifted per subject) |

## NACC - National Alzheimer's Coordinating Center Uniform Data Set

* Portal: https://naccdata.org (investigator data request)
* Raw folder: `<root>/Tables_NACC/`
* Output: `<root>/Tables_NACC/Integrated_NACC/NACC_longitudinal.csv`
* Subject ID: `NACCID`; time axis: date


| Table | File pattern | Required | ID column | Date column | Notes |
|---|---|---|---|---|---|
| `main` | `investigator_ftldlbd_nacc*.csv` | yes | NACCID | VISITYR | UDS visits (large; read in chunks) |
| `csf` | `investigator_fcsf_nacc*.csv` | yes | NACCID |  |  |
| `mri` | `investigator_scan_mriqc_nacc*.csv` | yes | NACCID |  |  |
| `pet_gaain` | `investigator_scan_amyloidpetgaain_nacc*.csv` | yes | NACCID | SCANDATE |  |
| `pet_npdka` | `investigator_scan_amyloidpetnpdka_nacc*.csv` | yes | NACCID | SCANDATE |  |
| `pet_fdg` | `investigator_scan_fdgpetnpdka_nacc*.csv` | yes | NACCID | SCANDATE |  |
| `pet_qc` | `investigator_scan_petqc_nacc*.csv` | yes | NACCID |  |  |
| `pet_tau` | `investigator_scan_taupetnpdka_nacc*.csv` | yes | NACCID | SCANDATE |  |

## NIFD - Frontotemporal Lobar Degeneration Neuroimaging Initiative (NIFD/FTLDNI)

* Portal: LONI IDA - https://ida.loni.usc.edu
* Raw folder: `<root>/Tables_NIFD/`
* Output: `<root>/Tables_NIFD/Integrated_NIFD/NIFD_longitudinal.csv`
* Subject ID: `LONI_ID`; time axis: date
* LONI IDA search export: `idaSearch_*.csv` in the raw folder (required)

| Table | File pattern | Required | ID column | Date column | Notes |
|---|---|---|---|---|---|
| `clinical` | `NIFD_Clinical_Data_{date}*.xlsx` | yes | LONI_ID | CLINICAL_LINKDATE | Clinical spreadsheet, one row per visit |

## OASIS - Open Access Series of Imaging Studies (OASIS-3)

* Portal: https://sites.wustl.edu/oasisbrains (NITRC-IR / XNAT Central)
* Raw folder: `<root>/Tables_OASIS/`
* Output: `<root>/Tables_OASIS/Integrated_OASIS/OASIS_longitudinal.csv`
* Subject ID: `OASISID`; time axis: days


| Table | File pattern | Required | ID column | Date column | Notes |
|---|---|---|---|---|---|
| `dx` | `OASIS3_UDSd1_diagnoses.csv` | yes | OASISID | days_to_visit |  |
| `cdr` | `OASIS3_UDSb4_cdr.csv` | yes | OASISID | days_to_visit |  |
| `gds` | `OASIS3_UDSb6_gds.csv` | yes | OASISID | days_to_visit |  |
| `neuro` | `OASIS3_UDSc1_cognitive_assessments.csv` | yes | OASISID | days_to_visit |  |
| `demo_a1` | `OASIS3_UDSa1_participant_demo.csv` | yes | OASISID | days_to_visit |  |
| `vitals` | `OASIS3_UDSb1_physical_eval.csv` | yes | OASISID | days_to_visit |  |
| `health` | `OASIS3_UDSa5*.csv` | no | OASISID | days_to_visit | UDS A5 health history (seizures, diabetes) |
| `motor` | `OASIS3_UDSb3.csv` | yes | OASISID | days_to_visit |  |
| `demog` | `OASIS3_demographics.csv` | yes | OASISID |  |  |
| `mri` | `OASIS3_MR_json.csv` | yes | subject_id |  |  |
| `pet` | `OASIS3_PET_json.csv` | yes | subject_id |  |  |
| `centiloid` | `OASIS3_amyloid_centiloid.csv` | yes | subject_id |  |  |

## PREVENTAD - PRe-symptomatic EValuation of Experimental or Novel Treatments for AD (PREVENT-AD)

* Portal: https://registeredpreventad.loris.ca (registered data release)
* Raw folder: `<root>/Tables_PREVENT-AD/`
* Output: `<root>/Tables_PREVENT-AD/Integrated_PREVENTAD/PREVENTAD_longitudinal.csv`
* Subject ID: `CONP_CandID`; time axis: date


| Table | File pattern | Required | ID column | Date column | Notes |
|---|---|---|---|---|---|
| `mci` | `MCI_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date |  |
| `clinical_dx` | `Clinical_diagnosis_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date |  |
| `demographics` | `Demographics_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID |  |  |
| `genetics` | `Genetics_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID |  |  |
| `el_cdr_moca` | `EL_CDR_MoCA_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `cdr_fu` | `CDR_FU_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `ravlt` | `RAVLT_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `rbans` | `RBANS_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `tmt` | `TMT_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `dkefs` | `DKEFS-CWIT_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `ad8` | `AD8_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `aps` | `APS_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID |  |  |
| `smell` | `Smell_identification_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `bp_weight` | `BP_Pulse_Weight_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `el_medical` | `EL_Medical_history_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `el_caide` | `EL_CAIDE_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `med_use` | `Med_use_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `lab` | `Lab_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `pet_info` | `PET_info_participants_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID |  |  |
| `csf` | `CSF_proteins_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `plasma_4plex` | `Plasma_4plex_IPMS_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |
| `plasma_ptau217` | `Plasma_p-tau217_Registered_PREVENTAD_{date}.csv` | yes | CONP_CandID | Date_taken |  |

## 4RTNI - 4-Repeat Tauopathy Neuroimaging Initiative

* Portal: LONI IDA - https://ida.loni.usc.edu
* Raw folder: `<root>/Tables_4TRNI/`
* Output: `<root>/Tables_4TRNI/Integrated_4RTNI/4RTNI_longitudinal.csv`
* Subject ID: `SUBID`; time axis: date
* LONI IDA search export: `idaSearch_*.csv` in the raw folder (optional; the dataset tables are the fallback)

| Table | File pattern | Required | ID column | Date column | Notes |
|---|---|---|---|---|---|
| `data` | `4RTNI_DATA_{date}.csv` | yes | SUBID |  | One wide file, one row per subject; timepoints as column suffixes _0, _6MO, _12MO |

