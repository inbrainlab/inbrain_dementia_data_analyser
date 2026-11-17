# Variables

Generated from `configs/variables.yaml` by `tools/make_variable_docs.py`; do not edit by hand.

207 variables. *expand* = placed into time windows in the expanded master.

| Variable | Label | Domain | Expand | Sources |
|---|---|---|---|---|
| `TMT_A` | Trail Making Test part A, time (s) | executive_attention | yes | **ABVIB** npsych: `TRAILA`<br>**ADNI** neurobat: `TRAASCOR`<br>**NACC** main: coalesce(`TRAILA`, `OTRAILA`)<br>**OASIS** neuro: `tma`<br>**PREVENTAD** tmt: `TMT_A_time`<br>**4RTNI** data: `UDSTRAILS_<tp>_TRAILA` |
| `TMT_B` | Trail Making Test part B, time (s) | executive_attention | yes | **ABVIB** npsych: `TRAILB`<br>**ADNI** neurobat: `TRABSCOR`<br>**NACC** main: coalesce(`TRAILB`, `OTRAILB`)<br>**OASIS** neuro: `tmb`<br>**PREVENTAD** tmt: `TMT_B_time`<br>**4RTNI** data: `UDSTRAILS_<tp>_TRAILB` |
| `BNT` | Boston Naming Test, total correct | language | yes | **ABVIB** npsych: `BOSTON`<br>**ADNI** neurobat: `BNTTOTAL`<br>**NACC** main: `BOSTON`<br>**NIFD** clinical: `BNTCORR` (date `CLINICAL_LINKDATE`)<br>**OASIS** neuro: `bnt` |
| `MINT` | Multilingual Naming Test, total | language |  | **ADNI** neurobat: `MINTTOTAL`<br>**NACC** main: `MINTTOTS`<br>**OASIS** neuro: `mintscnc` |
| `RAVLT_T1` | RAVLT trial 1 (words recalled, 0-15) | memory |  | **ADNI** neurobat: `AVTOT1`<br>**NACC** main: `REY1REC`<br>**PREVENTAD** ravlt: `ListA_Recall1_Correct_Answers` |
| `RAVLT_I` | RAVLT immediate recall, total of trials 1-5 (0-75) | memory |  | **ADNI** neurobat: sum(`AVTOT1`, `AVTOT2`, `AVTOT3`, `AVTOT4`, ...)<br>**NACC** main: sum(`REY1REC`, `REY2REC`, `REY3REC`, `REY4REC`, ...)<br>**PREVENTAD** ravlt: `List_A_Sum_Trials_1_5` |
| `RAVLT_D` | RAVLT delayed recall after ~30 min (words recalled, 0-15) | memory |  | **ADNI** neurobat: `AVDEL30MIN`<br>**NACC** main: `REYDREC`<br>**PREVENTAD** ravlt: `ListA_Delayed_Recall_Correct_Answer` |
| `LM_I` | Logical Memory, immediate recall | memory | yes | **A4** coglogic: `LIMMTOTAL`<br>**ABVIB** npsych: `LOGIMEM`<br>**ADNI** neurobat: `LIMMTOTAL`<br>**AIBL** neurobat: `LIMMTOTAL`<br>**NACC** main: `LOGIMEM`<br>**OASIS** neuro: `LOGIMEM` |
| `LM_II` | Logical Memory, delayed recall | memory | yes | **A4** coglogic: `LDELTOTAL`<br>**ABVIB** npsych: `MEMUNITS`<br>**ADNI** neurobat: `LDELTOTAL`<br>**AIBL** neurobat: `LDELTOTAL`<br>**NACC** main: `MEMUNITS`<br>**OASIS** neuro: `MEMUNITS` |
| `Clock_Command` | Clock drawing to command | visuospatial |  | **ADNI** neurobat: `CLOCKSCOR` |
| `Clock_Copy` | Clock copy | visuospatial |  | **ADNI** neurobat: `COPYSCOR` |
| `DS_Fwd` | Digit span forward | attention |  | **ABVIB** npsych: `DIGIF`<br>**ADNI** neurobat: `DSPANFOR`<br>**NACC** main: `DIGIF`<br>**NIFD** clinical: `DIGITFW` (date `CLINICAL_LINKDATE`)<br>**OASIS** neuro: `DIGIF`<br>**4RTNI** data: `DSPAN_<tp>_DIGITFW` |
| `DS_Back` | Digit span backward | attention |  | **ABVIB** npsych: `DIGIB`<br>**ADNI** neurobat: `DSPANBAC`<br>**NACC** main: `DIGIB`<br>**NIFD** clinical: `DIGITBW` (date `CLINICAL_LINKDATE`)<br>**OASIS** neuro: `DIGIB`<br>**4RTNI** data: `DSPAN_<tp>_DIGITBW` |
| `DigSym` | Digit Symbol Substitution / WAIS coding | processing_speed | yes | **A4** cogdigit: `DIGITTOTAL`<br>**ABVIB** npsych: `WAIS`<br>**ADNI** neurobat: `DIGITSCOR`<br>**NACC** main: `WAIS`<br>**OASIS** neuro: `digsym`<br>**4RTNI** data: `DSYM_<tp>` |
| `VF_anim` | Category fluency, animals | language | yes | **ABVIB** npsych: `ANIMALS`<br>**ADNI** neurobat: `CATANIMSC`<br>**NACC** main: `ANIMALS`<br>**NIFD** clinical: `ANCORR` (date `CLINICAL_LINKDATE`)<br>**OASIS** neuro: `ANIMALS`<br>**4RTNI** data: `VERBFLU_AN_<tp>_ANCORR` |
| `VF_veg` | Category fluency, vegetables | language |  | **ABVIB** npsych: `VEG`<br>**ADNI** neurobat: `CATVEGESC`<br>**OASIS** neuro: `VEG` |
| `MMSE` | Mini-Mental State Examination, total | global_cognition | yes | **A4** mmse: `MMSCORE`<br>**ABVIB** npsych: `MMSE`<br>**ADNI** mmse: `MMSCORE`<br>**AIBL** mmse: `MMSCORE`<br>**MCSA** data: `MMSEcalc`<br>**NACC** main: `NACCMMSE`<br>**NIFD** clinical: `MMSE_TOT` (date `MMSE_DCDATE`)<br>**OASIS** cdr: `MMSE`<br>**4RTNI** data: `MMSE_<tp>_MMSETOT` |
| `MoCA` | Montreal Cognitive Assessment, total | global_cognition | yes | **ADNI** moca: `MOCA`<br>**NACC** main: `NACCMOCA`<br>**NIFD** clinical: `MOCA_TotWithEduc` (date `MOCA_DCDATE`)<br>**OASIS** neuro: `mocatots`<br>**PREVENTAD** el_cdr_moca: `MOCA_total_score`<br>**4RTNI** data: `MOCA_<tp>_MOCATOTWITHEDUC` |
| `DS_Fwd_MoCA` |  | attention |  | **ADNI** moca: `DIGFOR` |
| `DS_Back_MoCA` |  | attention |  | **ADNI** moca: `DIGBACK` |
| `CDR` | Clinical Dementia Rating, global | global_staging | yes | **A4** cdr: `CDGLOBAL`<br>**ABVIB** cdr: `CDRGLOB`<br>**ADNI** cdr: `CDGLOBAL`<br>**AIBL** cdr: `CDGLOBAL`<br>**MCSA** data: `CDRGLOB`<br>**NACC** main: `CDRGLOB`<br>**NIFD** clinical: `CDR_TOT` (date `CDR_DCDATE`)<br>**OASIS** cdr: `CDRTOT`<br>**PREVENTAD** el_cdr_moca: `CDR_total_score`<br>**4RTNI** data: `CDR_<tp>_CDRTOT` |
| `CDR_SB` | Clinical Dementia Rating, sum of boxes | global_staging | yes | **A4** cdr: `CDSOB`<br>**ABVIB** cdr: `CDRSUM`<br>**ADNI** cdr: `CDRSB`<br>**MCSA** data: `CDRSUM`<br>**NACC** main: `CDRSUM`<br>**NIFD** clinical: `CDR_BOX_SCORE` (date `CDR_DCDATE`)<br>**OASIS** cdr: `CDRSUM`<br>**4RTNI** data: `CDR_<tp>_BOXSCORE` |
| `GDS` | Geriatric Depression Scale | mood | yes | **A4** psychwell: `GDTOTAL`<br>**ABVIB** gds: `GDS`<br>**ADNI** gds: `GDTOTAL`<br>**NACC** main: `NACCGDS`<br>**OASIS** gds: `GDS` |
| `FAQ` | Functional Activities Questionnaire, total | function | yes | **ADNI** faq: `FAQTOTAL`<br>**MCSA** data: `faq_total`<br>**NIFD** clinical: `FAQ_TOTAL` (date `FAQ_DCDATE`)<br>**4RTNI** data: `FAQ_<tp>_FAQTOT` |
| `NPI` | Neuropsychiatric Inventory (Questionnaire), total | behaviour | yes | **ADNI** npiq: `NPISCORE`<br>**NIFD** clinical: sum(`NPI_DELSEV`, `NPI_HALLSEV`, `NPI_AGITSEV`, `NPI_DEPDSEV`, ...) (date `NPI_DCDATE`)<br>**4RTNI** data: sum_timepoint_items of `NPIQ_<tp>_<item>` (12 items) |
| `STAI` | State-Trait Anxiety Inventory (ADNI 6-item sum) | mood |  | **A4** psychwell: `STAITOTAL`<br>**ADNI** stai: sum(`STAIAD1`, `STAIAD2`, `STAIAD3`, `STAIAD4`, ...) |
| `ECOG_PT` | Everyday Cognition 12-item, participant (sum of answered items; 9 = don't know excluded) | function |  | **ADNI** ecog_pt: sum(`ECOGPT1`, `ECOGPT2`, `ECOGPT3`, `ECOGPT4`, ...) (date `COLDATE`) |
| `ECOG_SP` | Everyday Cognition 12-item, study partner (sum of answered items; 9 = don't know excluded) | function |  | **ADNI** ecog_sp: sum(`ECOGSP1`, `ECOGSP2`, `ECOGSP3`, `ECOGSP4`, ...) (date `COLDATE`) |
| `Height` | Height (units of the source) | vitals | yes | **A4** vitals: `VSHEIGHT`<br>**ADNI** vitals: `VSHEIGHT`<br>**MCSA** data: `HT`<br>**NACC** main: `HEIGHT` (sentinels removed)<br>**OASIS** vitals: `HEIGHT` |
| `Weight` | Weight (units of the source) | vitals | yes | **A4** vitals: `VSWEIGHT`<br>**ADNI** vitals: `VSWEIGHT`<br>**MCSA** data: `WT`<br>**NACC** main: `WEIGHT` (sentinels removed)<br>**OASIS** vitals: `WEIGHT`<br>**PREVENTAD** bp_weight: `Weight` |
| `BP_Sys` | Systolic blood pressure (mmHg) | vitals |  | **A4** vitals: `VSBPSYS`<br>**ADNI** vitals: `VSBPSYS`<br>**MCSA** data: `SysBPAvg`<br>**NACC** main: `BPSYS`<br>**OASIS** vitals: `BPSYS`<br>**PREVENTAD** bp_weight: `Systolic_blood_pressure` |
| `BP_Dias` | Diastolic blood pressure (mmHg) | vitals |  | **A4** vitals: `VSBPDIA`<br>**ADNI** vitals: `VSBPDIA` or `VSBPDIAS`<br>**MCSA** data: `DiaBPAvg`<br>**NACC** main: `BPDIAS`<br>**OASIS** vitals: `BPDIAS`<br>**PREVENTAD** bp_weight: `Diastolic_blood_pressure` |
| `BMI` | Body-mass index (kg/m2) | vitals | yes | **A4** vitals: bmi of `VSWEIGHT` / `VSHEIGHT`<br>**ADNI** vitals: bmi of `VSWEIGHT` / `VSHEIGHT`<br>**MCSA** data: `BMI`<br>**NACC** main: bmi_reported_or_computed of `WEIGHT` / `HEIGHT`<br>**OASIS** vitals: bmi of `WEIGHT` / `HEIGHT` |
| `VF_D` |  | language |  | **NIFD** clinical: `DCORR` (date `CLINICAL_LINKDATE`)<br>**4RTNI** data: `VERBFLU_DWORDS_<tp>_DCORR` |
| `CVLT_total` |  | memory |  | **NIFD** clinical: `TRCOTOT` (date `CLINICAL_LINKDATE`)<br>**4RTNI** data: `CVLT_<tp>_TRCOTOT` |
| `CVLT_30s` |  | memory |  | **NIFD** clinical: `CORR30` (date `CLINICAL_LINKDATE`)<br>**4RTNI** data: `CVLT_<tp>_CORR30` |
| `CVLT_10m` |  | memory |  | **NIFD** clinical: `CORR10` (date `CLINICAL_LINKDATE`)<br>**4RTNI** data: `CVLT_<tp>_CORR10` |
| `CVLT_recog` |  | memory |  | **NIFD** clinical: `RECOG` (date `CLINICAL_LINKDATE`)<br>**4RTNI** data: `CVLT_<tp>_RECOG` |
| `ModTrails_T` |  | executive_attention |  | **NIFD** clinical: `MTTTIME` (date `CLINICAL_LINKDATE`)<br>**4RTNI** data: `MTRAILS_<tp>_MTTIME` |
| `ModTrails_C` |  | executive_attention |  | **NIFD** clinical: `MTCORR` (date `CLINICAL_LINKDATE`)<br>**4RTNI** data: `MTRAILS_<tp>_MTCORR` |
| `GDS_30` |  | mood |  | **NIFD** clinical: `GDSTOT` (date `CLINICAL_LINKDATE`)<br>**4RTNI** data: `GDS_<tp>_GDSTOT` |
| `GDS_15` |  | mood |  | **NIFD** clinical: `GDS15TO` (date `CLINICAL_LINKDATE`)<br>**4RTNI** data: `GDS_<tp>_GDS15TO` |
| `SEADL` |  | function |  | **NIFD** clinical: `SEADL` (date `SEADL_DCDATE`)<br>**4RTNI** data: `SEADL_<tp>` |
| `CGI_sev` |  | global_cognition |  | **NIFD** clinical: `CGI_SEVERITY` (date `CGI_DCDATE`)<br>**4RTNI** data: `CGI_S_<tp>` |
| `CGI_chg` |  | global_cognition |  | **NIFD** clinical: `CGI_CHANGE` (date `CGI_DCDATE`)<br>**4RTNI** data: `CGI_C_<tp>` |
| `PSPRS` |  | motor |  | **NIFD** clinical: `PSPRS_TOTAL` (date `PSPRS_DC_DATE`)<br>**4RTNI** data: `PSPRS_<tp>_TOTAL` |
| `BIS` |  | behaviour |  | **NIFD** clinical: `BIS_TOT` (date `BISBAS_DCDATE`)<br>**4RTNI** data: `BIS_<tp>_BIS_TOT` |
| `BAS` |  | behaviour |  | **NIFD** clinical: `BAS_TOT` (date `BISBAS_DCDATE`) |
| `IRI_EC` |  | behaviour |  | **NIFD** clinical: `IRI_FTDIRIEC` (date `IRI_DCDATE`)<br>**4RTNI** data: `IRI_<tp>_IRIEC` |
| `IRI_PT` |  | behaviour |  | **NIFD** clinical: `IRI_FTDIRIPT` (date `IRI_DCDATE`)<br>**4RTNI** data: `IRI_<tp>_IRIPT` |
| `UPDRS` |  | motor |  | **MCSA** data: `TOT_MOD_UPDRS`<br>**NIFD** clinical: `MDEXAM_UPDRS` (date `MDEXAM_DCDATE`)<br>**4RTNI** data: `UPDRS_<tp>_TOTAL` |
| `UHDRS` |  | motor |  | **NIFD** clinical: `MDEXAM_UHDRS` (date `MDEXAM_DCDATE`) |
| `LM_I_cog` |  | memory |  | **ABVIB** cogscor: `LMIAI` |
| `VF_anim_cog` |  | language |  | **ABVIB** cogscor: `CFTOTANI` |
| `VF_F_cog` |  | language |  | **ABVIB** cogscor: `CFTOTALF` |
| `DigSym_cog` |  | attention |  | **ABVIB** cogscor: `WAISRTOT` |
| `BNT_30` |  | language |  | **ABVIB** cogscor: `BNTCORR30` |
| `BNT_60` |  | language |  | **ABVIB** cogscor: `BNTCORR60` |
| `TMT_A_cog` |  | executive_attention |  | **ABVIB** cogscor: `TRAILSA` |
| `TMT_B_cog` |  | executive_attention |  | **ABVIB** cogscor: `TRAILSB` |
| `GDS_cog` |  | mood |  | **ABVIB** cogscor: `GDS` |
| `VF_A_cog` |  | language |  | **ABVIB** cogscor: `VFTOTALA` |
| `VF_S_cog` |  | language |  | **ABVIB** cogscor: `VFTOTALS` |
| `Glucose` | Blood glucose | lab | yes | **ABVIB** lab: `GLUCOSE` |
| `Cholesterol` | Total cholesterol | lab | yes | **ABVIB** lab: `CHOL`<br>**PREVENTAD** lab: `total_cholesterol_value` |
| `LDL` | LDL cholesterol | lab |  | **ABVIB** lab: `LDL`<br>**PREVENTAD** lab: `LDL_value` |
| `HDL` | HDL cholesterol | lab |  | **ABVIB** lab: `HDL`<br>**PREVENTAD** lab: `HDL_value` |
| `hsCRP` |  | lab |  | **ABVIB** lab: `HSCRP` |
| `HbA1c` | Glycated haemoglobin | lab | yes | **ABVIB** lab: `HGBA1C`<br>**PREVENTAD** lab: `hba1c_value` |
| `Fasting` |  | lab | yes | **ABVIB** subcat: `FASTING` |
| `Diabetes` |  | health_lifestyle | yes | **ABVIB** subcat: `DIABETS`<br>**MCSA** data: `DiabetesGS` |
| `STMS` |  | global_cognition |  | **MCSA** data: `SHTOTSC` |
| `pz_memory` |  | global_cognition |  | **MCSA** data: `pzmemory` |
| `pz_language` |  | global_cognition |  | **MCSA** data: `pzlanguage` |
| `pz_attention` |  | global_cognition |  | **MCSA** data: `pzattention` |
| `pz_visuospatial` |  | global_cognition |  | **MCSA** data: `pzvisualsp` |
| `pz_global` |  | global_cognition |  | **MCSA** data: `pzglobal` |
| `SmellTest` |  | other |  | **MCSA** data: `SMELLTOT` |
| `ECOG` |  | function |  | **MCSA** data: `ECOG_Global` |
| `ECOG_12` |  | function |  | **MCSA** data: `ECOG_12item` |
| `BAI` |  | mood |  | **MCSA** data: `BAItot` |
| `BDI` |  | mood |  | **MCSA** data: `BDITOT` |
| `NPI_count` |  | behaviour |  | **MCSA** data: `n_npsych_symp` |
| `Blessed` |  | global_cognition |  | **MCSA** data: `BLESSTOT` |
| `Hachinski` |  | other |  | **MCSA** data: `HACHTOT` |
| `ESS_informant` |  | other |  | **MCSA** data: `ESSscore` |
| `ESS_subject` |  | other |  | **MCSA** data: `ESSscore_Subject` |
| `SmokingStatus` |  | health_lifestyle |  | **MCSA** data: `SmokeStat` |
| `Exercise_PastYear` |  | health_lifestyle |  | **MCSA** data: `EXWkEq12score` |
| `Exercise_Midlife` |  | health_lifestyle |  | **MCSA** data: `EXWkEq50score` |
| `CognSoc_PastYear` |  | health_lifestyle |  | **MCSA** data: `CSWkEq12score` |
| `CognSoc_Midlife` |  | health_lifestyle |  | **MCSA** data: `CSWkEq50score` |
| `Diet_Calories` |  | health_lifestyle |  | **MCSA** data: `Cals_kcal_` |
| `Diet_Fat_Cal` |  | health_lifestyle |  | **MCSA** data: `Cals_from_fat_kcal_` |
| `Hippo_Vol` |  | imaging_derived |  | **MCSA** data: `HPVOL` (date `ImagingDate`) |
| `TIV` |  | imaging_derived |  | **MCSA** data: `TIV` (date `ImagingDate`) |
| `HVA` |  | imaging_derived |  | **MCSA** data: `HVA` (date `ImagingDate`) |
| `PIB_SUVR` |  | imaging_derived |  | **MCSA** data: `PIB_SUVR` (date `ImagingDate`) |
| `PIB_Centiloid` |  | imaging_derived |  | **MCSA** data: `PIB_CENTILOID` (date `ImagingDate`) |
| `PIB_Abnormal` |  | imaging_derived |  | **MCSA** data: `ABNORMAL_PIB` (date `ImagingDate`) |
| `Hypertension` |  | health_lifestyle |  | **MCSA** data: `HypertensionTxGS` |
| `Dyslipidemia` |  | health_lifestyle |  | **MCSA** data: `DyslipidemiaTxGS` |
| `MI` |  | health_lifestyle |  | **MCSA** data: `MIGS` |
| `Stroke` |  | health_lifestyle |  | **MCSA** data: `StrokeGS` |
| `Diabetes_med_Insulin` |  | health_lifestyle | yes | **MCSA** data: `InsulinTx` |
| `Diabetes_med_Metformin` |  | health_lifestyle | yes | **MCSA** data: `Metformin` |
| `Diabetes_med_OralAgent` |  | health_lifestyle | yes | **MCSA** data: `OralAgent` |
| `Diabetes_med_OralOther` |  | health_lifestyle | yes | **MCSA** data: `OralAgentOTMF` |
| `Diabetes_med` | Current diabetes medication | medication | yes | **MCSA** data: max(`InsulinTx`, `Metformin`, `OralAgent`, `OralAgentOTMF`)<br>**NACC** main: `NACCDBMD` (codes [0, 1, 2]) |
| `PSPRS_History` |  | motor |  | **4RTNI** data: `PSPRS_<tp>_SUBSCORE_HISTORY` |
| `PSPRS_Mentation` |  | motor |  | **4RTNI** data: `PSPRS_<tp>_SUBSCORE_MENTATION` |
| `PSPRS_Bulbar` |  | motor |  | **4RTNI** data: `PSPRS_<tp>_SUBSCORE_BULBAR` |
| `PSPRS_Ocular` |  | motor |  | **4RTNI** data: `PSPRS_<tp>_SUBSCORE_OCULARMOTOR` |
| `PSPRS_Limb` |  | motor |  | **4RTNI** data: `PSPRS_<tp>_SUBSCORE_LIMBMOTOR` |
| `PSPRS_Gait` |  | motor |  | **4RTNI** data: `PSPRS_<tp>_SUBSCORE_GAITMIDLINE` |
| `UPDRS_PDNormal` |  | motor |  | **OASIS** motor: `PDNORMAL`<br>**4RTNI** data: `UPDRS_<tp>_PDNORMAL` |
| `UPDRS_Imputed` |  | motor |  | **4RTNI** data: `UPDRS_<tp>_IMPUTED` |
| `CVLT_cued` |  | memory |  | **4RTNI** data: `CVLT_<tp>_CUEDCOR` |
| `CVLT_long` |  | memory |  | **4RTNI** data: `CVLT_<tp>_CORRLONG` |
| `Benson_copy` | Benson complex figure, copy (modified Rey, 0-17) | visuospatial |  | **NIFD** clinical: `MODREY` (date `CLINICAL_LINKDATE`)<br>**4RTNI** data: `BENSON_<tp>_MODREY` |
| `Benson_10m` | Benson complex figure, 10-minute delayed recall (0-17) | memory |  | **NIFD** clinical: `REY10M` (date `CLINICAL_LINKDATE`)<br>**4RTNI** data: `BENSON_<tp>_REY10M` |
| `Benson_recog` |  | memory |  | **4RTNI** data: `BENSON_<tp>_REYRECG` |
| `BNT_total` |  | language |  | **4RTNI** data: `BNT_<tp>_BNTTOT` |
| `BNT_corr` |  | language |  | **4RTNI** data: `BNT_<tp>_BNTCORR` |
| `Stroop_CW` |  | executive_attention |  | **4RTNI** data: `STROOP_<tp>_STRPCNCOR` |
| `Stroop_C` |  | executive_attention |  | **4RTNI** data: `STROOP_<tp>_STRPCOR` |
| `VOSP` |  | other |  | **4RTNI** data: `VOSP_<tp>_NUMBLOC` |
| `WRAT4` |  | other |  | **4RTNI** data: `WRAT4_<tp>_WRATTOT` |
| `CDR_Memory` |  | global_staging |  | **A4** cdr: `MEMORY`<br>**OASIS** cdr: `memory`<br>**4RTNI** data: `CDR_<tp>_MEMORY` |
| `CDR_Behavior` |  | global_staging |  | **4RTNI** data: `CDR_<tp>_BEHAV` |
| `CDR_Motor` |  | global_staging |  | **4RTNI** data: `CDR_<tp>_MOTOR` |
| `SNQ22` |  | behaviour |  | **4RTNI** data: `SNQ22_<tp>_QTOTC` |
| `RSMS` |  | behaviour |  | **4RTNI** data: `RSMS_<tp>_RSMSTOTI` |
| `RSMS_SelfPres` |  | behaviour |  | **4RTNI** data: `RSMS_<tp>_RSMSSUB1I` |
| `RSMS_Expressive` |  | behaviour |  | **4RTNI** data: `RSMS_<tp>_RSMSSUB2I` |
| `CDR_Orient` |  | global_staging |  | **A4** cdr: `ORIENT`<br>**OASIS** cdr: `orient` |
| `FCSR_FreeRecall` |  | memory |  | **A4** cogfcsr16: `FCTOTF` |
| `FCSR_CuedRecall` |  | memory |  | **A4** cogfcsr16: `FCTOTC` |
| `FCSR_Total` |  | memory |  | **A4** cogfcsr16: `FCTOTAL` |
| `CFI_patient` |  | function |  | **A4** cfi: `CFIPTTOTAL` |
| `CFI_partner` |  | function |  | **A4** cfisp: `CFSPTTOTAL` |
| `IES` |  | mood |  | **A4** ies: `IESCORE` |
| `Pulse` | Pulse (beats/min) | vitals |  | **A4** vitals: `VSPULSE`<br>**NACC** main: `HRATE`<br>**OASIS** vitals: `HRATE`<br>**PREVENTAD** bp_weight: `Pulse` |
| `Sleep_hrs` |  | health_lifestyle |  | **A4** habits: `SLEEP` |
| `Exercise_aerobic` |  | health_lifestyle |  | **A4** habits: `AEROBIC` |
| `Exercise_walking` |  | health_lifestyle |  | **A4** habits: `WALKING` |
| `Alcohol` |  | health_lifestyle |  | **A4** habits: `ALCOHOL` |
| `Caffeine` |  | health_lifestyle |  | **A4** habits: `CAFFEINE` |
| `Smoking` |  | health_lifestyle |  | **A4** habits: `SMOKE` |
| `ADLPQ` | ADL Prevention Questionnaire, total | function |  | **A4** adlpq: sum(`ASCBOOK`, `ASTRAVL`, `ASAPPLI`, `ASLAUN`, ...) |
| `PACC` | Preclinical Alzheimer Cognitive Composite | global_cognition |  | **A4** sppacc: `PACCRN` (rows where `PACCNAME` ~ `composite|PACC|total`) |
| `MoCA_DelayedRecall` |  | global_cognition |  | **PREVENTAD** el_cdr_moca: `delayed_recall_score` |
| `MoCA_Orientation` |  | global_cognition |  | **PREVENTAD** el_cdr_moca: `orientation_score` |
| `CDR_FU` |  | global_staging |  | **PREVENTAD** cdr_fu: `CDR_total_score` |
| `RAVLT_Recog` | RAVLT delayed recognition (ADNI AVDELTOT / NACC REYTCOR; PREVENT-AD words recognised) | memory |  | **ADNI** neurobat: `AVDELTOT`<br>**NACC** main: `REYTCOR`<br>**PREVENTAD** ravlt: `List_A_Recognition_Recognized` |
| `RBANS_total` |  | global_cognition |  | **PREVENTAD** rbans: `total_scale_index_score` |
| `RBANS_ImmMem` |  | global_cognition |  | **PREVENTAD** rbans: `immediate_memory_index_score` |
| `RBANS_DelMem` |  | global_cognition |  | **PREVENTAD** rbans: `delayed_memory_index_score` |
| `RBANS_Attn` |  | global_cognition |  | **PREVENTAD** rbans: `attention_index_score` |
| `DS_RBANS` |  | attention |  | **PREVENTAD** rbans: `digit_span_score` |
| `VF_anim_RBANS` |  | language |  | **PREVENTAD** rbans: `semantic_fluency_score` |
| `DigSym_RBANS` |  | attention |  | **PREVENTAD** rbans: `coding_score` |
| `Stroop_CN` |  | executive_attention |  | **PREVENTAD** dkefs: `Color_Naming_Time` |
| `Stroop_Inhib` |  | executive_attention |  | **PREVENTAD** dkefs: `Inhibition_Time` |
| `Stroop_Inhib_scaled` |  | executive_attention |  | **PREVENTAD** dkefs: `Inhibition_Scaled_Score` |
| `Stroop_Switch_scaled` |  | executive_attention |  | **PREVENTAD** dkefs: `Switching_Scaled_Score` |
| `AD8` |  | global_cognition |  | **PREVENTAD** ad8: `AD8_total_score` |
| `APS` |  | other |  | **PREVENTAD** aps: `APS_score` |
| `UPSIT` |  | other |  | **PREVENTAD** smell: `UPSIT_total_score` |
| `CAIDE` |  | other |  | **PREVENTAD** el_caide: `CAIDE_total_score` |
| `Diabetes_treatment` |  | health_lifestyle | yes | **PREVENTAD** el_medical: `treatment_diabetes` |
| `CAIDE_BMI` |  | other | yes | **PREVENTAD** el_caide: `CAIDE_BMI_subscore` |
| `CDR_Judgment` |  | global_staging |  | **OASIS** cdr: `judgment` |
| `CDR_Community` |  | global_staging |  | **OASIS** cdr: `commun` |
| `CDR_Home` |  | global_staging |  | **OASIS** cdr: `homehobb` |
| `CDR_PersonalCare` |  | global_staging |  | **OASIS** cdr: `perscare` |
| `TMT_A_UDS` |  | executive_attention |  | **OASIS** neuro: `TRAILA` |
| `LM_Delay` |  | memory |  | **OASIS** neuro: `lmdelay` |
| `Craft_I` |  | memory |  | **NACC** main: `CRAFTVRS`<br>**OASIS** neuro: `craftvrs` |
| `Craft_II` |  | memory |  | **OASIS** neuro: `craftdvr` |
| `UPDRS_Brady` |  | motor |  | **OASIS** motor: `BRADYKIN` |
| `SEIZURES` | UDS A5 seizures (0 absent, 1 recent/active, 2 remote) | health_history | yes | **NACC** main: `SEIZURES` (codes [0, 1, 2])<br>**OASIS** health: `SEIZURES` (codes [0, 1, 2]) |
| `DIABETES` | UDS A5 diabetes (0 absent, 1 recent/active, 2 remote) | health_history | yes | **NACC** main: `DIABETES` (codes [0, 1, 2])<br>**OASIS** health: `DIABETES` (codes [0, 1, 2]) |
| `DIABTYPE` | UDS A5 diabetes type (1, 2, 3 other) | health_history | yes | **NACC** main: `DIABTYPE` (codes [1, 2, 3])<br>**OASIS** health: `DIABTYPE` (codes [1, 2, 3]) |
| `EPILEP` | UDS D1 epilepsy diagnosis (0/1) | health_history | yes | **NACC** main: `EPILEP` (codes [0, 1])<br>**OASIS** dx: `EPILEP` (codes [0, 1]) |
| `EPILEPIF` | UDS D1 epilepsy contribution (1 primary, 2 contributing, 3 non-contributing) | health_history | yes | **NACC** main: `EPILEPIF` (codes [1, 2, 3])<br>**OASIS** dx: `EPILEPIF` (codes [1, 2, 3]) |
| `marital_status` | Marital status (UDS code) | demographics |  | **OASIS** demo_a1: `MARISTAT` |
| `living_situation_v2` | Living situation (UDS v2 code) | demographics |  | **OASIS** demo_a1: `LIVSIT` |
| `living_situation_v3` | Living situation (UDS v3 code) | demographics |  | **OASIS** demo_a1: `LIVSITUA` |
| `independence_level` | Level of independence (UDS code) | demographics |  | **OASIS** demo_a1: `INDEPEND` |
| `residence_type` | Type of residence (UDS code) | demographics |  | **OASIS** demo_a1: `RESIDENC` |
| `LM_time` |  | memory |  | **NACC** main: `MEMTIME` |
| `Craft_D` |  | memory |  | **NACC** main: `CRAFTDVR` |
| `Antiepileptic_med` | Current antiepileptic medication (NACC derived) | medication | yes | **NACC** main: `NACCEPMD` (codes [0, 1, 2]) |
| `DIABET` | UDS D2 clinician-assessed diabetes (0 no, 1 type 1, 2 type 2, 3 other) | health_history | yes | **NACC** main: `DIABET` (codes [0, 1, 2, 3]) |
| `Neuro_history` | Neurologic history other than AD (1/0; ADNI MEDHIST, AIBL MH2NEURL, A4 INITHEALTH) | health_history | yes | *built by dataset code* |
| `Endo_history` | Endocrine-metabolic history (1/0; AIBL MH9ENDO, A4 INITHEALTH) | health_history | yes | *built by dataset code* |
| `ABETA42` | CSF amyloid-beta 1-42 (main CSF panel of the cohort) | csf | yes | *built by dataset code* |
| `ABETA40` | CSF amyloid-beta 1-40 | csf | yes | *built by dataset code* |
| `PTAU` | CSF phospho-tau 181 | csf | yes | *built by dataset code* |
| `TAU` | CSF total tau | csf | yes | *built by dataset code* |
| `ABETA42_lumipulse` | CSF amyloid-beta 1-42, Lumipulse (PREVENT-AD) | csf | yes | *built by dataset code* |
| `ABETA40_lumipulse` | CSF amyloid-beta 1-40, Lumipulse (PREVENT-AD) | csf | yes | *built by dataset code* |
| `PTAU217_CSF` | CSF phospho-tau 217 (PREVENT-AD) | csf | yes | *built by dataset code* |
| `Plasma_pTau217_simoa` | Plasma p-tau217, Simoa (PREVENT-AD) | plasma | yes | *built by dataset code* |
| `Plasma_pTau181_simoa` | Plasma p-tau181, Simoa (PREVENT-AD) | plasma | yes | *built by dataset code* |
| `Plasma_pTau231_simoa` | Plasma p-tau231, Simoa (PREVENT-AD) | plasma | yes | *built by dataset code* |

## Imaging types

| Type | Modality | Description | Expand |
|---|---|---|---|
| `MRI_T1w` | MRI | 3D T1-weighted MRI | yes |
| `MRI_T2w` | MRI | T2-weighted MRI | yes |
| `MRI_T2star` | MRI | T2*-weighted / SWI MRI | yes |
| `MRI_T2w_FLAIR` | MRI | T2 FLAIR MRI | yes |
| `MRI_ASL` | MRI | Arterial spin labelling | yes |
| `MRI_fMRI` | MRI | Resting-state / task fMRI | yes |
| `MRI_DTI` | MRI | Diffusion MRI | yes |
| `PET_FDG` | PET | FDG PET | yes |
| `ABETA_AV45` | PET | Amyloid PET, florbetapir (AV-45) | yes |
| `ABETA_FBB` | PET | Amyloid PET, florbetaben | yes |
| `ABETA_Flutemetamol` | PET | Amyloid PET, flutemetamol | yes |
| `ABETA_NAV4694` | PET | Amyloid PET, NAV4694 | yes |
| `ABETA_PIB` | PET | Amyloid PET, Pittsburgh compound B | yes |
| `TAU_Flortaucipir` | PET | Tau PET, flortaucipir (AV-1451) | yes |
| `TAU_MK6240` | PET | Tau PET, MK-6240 | yes |
| `TAU_PI2620` | PET | Tau PET, PI-2620 | yes |
| `MRI_T1w_1.5T` | MRI | AIBL 1.5 T T1-weighted MRI (kept apart from MRI_T1w) |  |
