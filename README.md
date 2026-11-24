# NRSC_510: IPV-BI Biomarker and Neurocognitive Analysis

**Author:** Courtney McCabe  
**Course:** NRSC 510  
**Affiliated lab:** Wellington Lab  
**Last updated:** November 23, 2026  
**Related Dataverse:** https://demo.borealisdata.ca/dataverse/courtmccabe (outcome data and project proposal information)

> **Disclaimer:** This repository analyzes a simulated dataset built to mirror the structure of the ASCEND-IPV Biomarkers study. It does not contain real participant data. No clinical interpretations should be drawn from any results produced by this code.

This README follows the GitHub guidance on repository READMEs (see [References](#references)). It documents what each script does, the order to run them in, the software needed, and the checkpoints that confirm a script ran correctly. It deliberately does not report analysis results; those are produced by running the code.

## Project overview

This project investigates whether plasma biomarker profiles and neuropsychological performance indicate elevated risk of neurological and psychological impairment in survivors of intimate partner violence-caused brain injury (IPV-BI) and whether this risk differs between acute and chronic timepoints.

The analysis compares three study groups:
- Acute/Subacute IPV-BI
- Acute/Subacute non-IPV-BI
- Chronic IPV-BI

### Project goals (all milestones complete)

- Build a reproducible tidyverse-based data cleaning and analysis workflow in R (Milestone 1)
- Summarize intake characteristics across neurological and non-neurological symptom domains (Milestone 1)
- Compare demographics, IPV-BI exposure history, injury mechanisms, symptoms, substance use, and psychiatric diagnoses across study groups using covariate-adjusted models with multiple comparison correction (Milestone 2)
- Test for group differences across the plasma biomarker panel with limma, cross-checked with a nonparametric test, with quality control checks and visualization (Milestone 3)
- Compare neuropsychological performance across groups, then relate the biomarker panel to neuropsychological outcomes using correlation and linear regression, first across all proteins (hypothesis generating) and then for targeted proteins and outcomes (hypothesis driven) (Milestone 4)

## Run order at a glance

The simulated v7 dataset (see [How the data were simulated](#how-the-data-were-simulated)) exists before any script runs and is the input to the first script. Run the scripts in this order. Restart R between scripts so each one starts from a clean session.

| Order | Script | Milestone | Needs | Produces |
|---|---|---|---|---|
| 1 | `Scripts/cleaning_summary_v3.R` | 1: cleaning and intake summaries | `IPV_cohort_simulated_v7.xlsx` | `ipv_cohort_clean.rds`, intake summary tables |
| 2 | `Scripts/group_comparisons_milestone_2.R` | 2: group comparisons | `ipv_cohort_clean.rds` | Comparison tables and age boxplot |
| 3 | `Scripts/heatmap_volcano_milestone_3_v2.R` | 3: biomarker group differences | `ipv_cohort_clean.rds`, `Neuro_220_Panel.csv` | `pairwise_results.rds`, updated `ipv_cohort_clean.rds`, figures and CSVs |
| 4 | `Scripts/volcano_plots_group_comparison.R` | 3 (supplement, optional) | `pairwise_results.rds` | Group-colored volcano plots |
| 5 | `Scripts/neuropsych_biomarker_milestone_4.R` | 4: neuropsychology and biomarker association | `ipv_cohort_clean.rds`, `pairwise_results.rds` | `milestone4_results.xlsx`, figures, updated `ipv_cohort_clean.rds` |

Milestone 2 only needs the cleaned file from Milestone 1, so it does not depend on Milestones 3 and 4. Milestone 3 writes additional columns back into `ipv_cohort_clean.rds`, so re-running Milestone 1 resets that file, and Milestone 3 should then be re-run before Milestone 4. Milestone 4 needs the `pairwise_results.rds` file made by Milestone 3.

Step numbers continue across the project (Steps 1 to 15, then Steps 16 to 35, then Steps 41 to 55). Steps 36 to 40 are not used in the final scripts.

## Data

### How the data were simulated

The simulation was done before any analysis. The v7 workbook it produced is the fixed input for every script, so results do not depend on re-running the simulation.

**Why the data were simulated**

- The project builds and tests an analysis workflow for plasma biomarkers and neuropsychological outcomes in IPV-BI. Participant-level data from survivors of intimate partner violence are sensitive, and the real study data are not used in this project. A simulated dataset lets the full workflow be built, shared openly, and reproduced with no risk to participants.
- A simulation with a known built-in structure lets the workflow be checked: group differences that were put into the data can be looked for, and the code can be tested from start to finish before it is pointed at real data.
- The structure follows a real study design, so the code is directly relevant to it.

**Structure and types of data simulated**

The structure follows the ASCEND-IPV study protocol (BMJ Open 2025, doi:10.1136/bmjopen-2024-098025). There are three study groups: Acute/Subacute IPV-BI, Acute/Subacute non-IPV-BI, and Chronic IPV-BI. There is no chronic non-IPV-BI group. Every participant has intake and exposure history, neuropsychological testing, and a plasma biomarker panel (Alamar NULISAseq Neuro 220). Only the structure is taken from the study; all values are simulated.

| Type | Examples |
|---|---|
| Continuous | Age, height, weight, neuropsychological scores, log2 NPQ biomarker values, months since first exposure, days to symptom resolution, plate quality control values |
| Ordinal | Education category, symptom presence (No, Unknown, Yes), substance use frequency, historic exposure count bins (0, 1-5, 6-10, 11-15, 15+) |
| Binary (Yes/No) | Incident type columns (non-exclusive), head impact and strangulation sign columns, psychiatric diagnoses, symptom persistence, plate quality control pass flags |
| Nominal categories | Sex, gender, ethnicity, sexual orientation, assailant relationship |
| Counts and totals | BISA and HELPS exposure measures, questionnaire totals (PCL-5, PHQ-9, NSI) |
| Missing values | Neuropsychological measures only, missing completely at random by design |

**Parameters used.** Each parameter is recorded in the workbook's `Parameters_Sources` sheet with its source. Values with no literature source are labelled SIMULATION CHOICE. Simulated values and ranges were approximated based on consulting with Wellington lab members and referencing mTBI literature; these prompts were given to Claude Sonnet 5.5 model.

| Parameter | Value used | Source or status |
|---|---|---|
| Group sizes and timing | Acute/subacute: 46 participants within 30 days (23 IPV-BI, 23 non-IPV-BI). Chronic: 20 IPV-BI, 6 months to more than 10 years after the incident. | Study design; the chronic non-IPV-BI group was removed per feedback. Chronic upper bound of 15 years is a SIMULATION CHOICE. |
| Age | Mean 40.99, SD 13.97, clipped to 18 to 75 | Makovec Knight et al. 2025 |
| Sex | 90% female | ASCEND-IPV protocol (predominantly women); 90% is a SIMULATION CHOICE |
| Sexual orientation | 79% heterosexual | Herbenick et al. 2025; split of the remainder is a SIMULATION CHOICE |
| Ethnicity, education, assailant relationship | Category probabilities set in the simulation script | SIMULATION CHOICE; levels follow the project intake form |
| IPV-BI mechanism mix | 38% head impact, 38% strangulation, 24% both | Cimino et al. 2019 |
| Exposure burden | 41% with more than 6 lifetime events, 59% with 1 to 6, split across head impact and strangulation | Makovec Knight et al. 2025; allocation between mechanisms is a SIMULATION CHOICE |
| Substance use | Probability of any use ordered alcohol, then marijuana, then painkillers | Ordering from Mehr et al. 2023; rates are SIMULATION CHOICES |
| Neuropsychological tests | Score scales follow each instrument's scoring (raw, T-score, standard score, z-score). IPV-BI lowers RAVLT and TMT-B scores, with a larger effect at more than 6 lifetime events. | Direction from Makovec Knight et al. 2025; group means and SDs are SIMULATION CHOICES |
| PCL-5, PHQ-9, NSI, TOMM | Symptom scales elevated in both groups and higher with injury burden; most participants pass TOMM Trial 2 | SIMULATION CHOICE |
| Biomarker panel | 222 targets from the Alamar NULISAseq Neuro 220 panel | Alamar panel annotation |
| Biomarker scale | Log2 NPQ. Each target has its own baseline mean drawn between 10.5 and 15.5 log2 units, with a within-target SD of about 1.0 to 1.5. | Output format of the NULISAseq assay; distributions are SIMULATION CHOICES |
| Biomarker group effects | 21 of the 222 targets carry a built-in IPV-BI shift, larger in the acute phase than in the chronic phase. The other 201 have no built-in group difference. Values are listed in the `Panel_Targets` sheet. | Direction set by the project team and informed by Esopenko et al. 2024 (acute GFAP and UCH-L1, prolonged NfL); magnitudes are SIMULATION CHOICES |
| Plate quality control | Eight columns for the negative control (NC), internal control (IC), and inter-plate control (IPC), built around Alamar's documented criteria: IC reads of at least 1000, IC reads within 60 to 140% of the plate median, IPC CV of at most 10%, NC mean below 3 | Alamar NULISAseqR User Guide |
| Missing data | About 4% per neuropsychological variable, missing completely at random. The biomarker panel has no missing values. | SIMULATION CHOICE |

**What the simulation can and cannot support.** It supports building, testing, demonstrating, and sharing the workflow: cleaning, quality control, group comparisons, limma modeling, and correlation and regression analyses. It does not support clinical or biological conclusions. 

### Current dataset: `IPV_cohort_simulated_v7.xlsx`

**Use version 7 (v7) for everything in this repository.** All scripts were written and tested against v7. Earlier versions are not for use; their history is under [Previous versions of the simulated data](#previous-versions-of-the-simulated-data).

| Sheet | Contents |
|---|---|
| `README` | Notes from the simulation. |
| `Acute_Subacute` | 46 participants within 30 days of the incident: 23 IPV-BI and 23 non-IPV-BI. 376 columns. |
| `Chronic` | 20 IPV-BI participants, from 6 months to more than 10 years after the incident. 376 columns. |
| `Data_Dictionary` | Definitions for the variables. |
| `Parameters_Sources` | Simulation parameters, values used, and source or status. |
| `Group_Summary` | Group summary table. |
| `Panel_Targets` | Annotation for each biomarker column: column ID, target name, UniProt ID, hallmark of neurodegenerative disease, and a column recording the effect built into the simulation for each target. No Milestone 3 or 4 script reads that column; it can be used afterward to check what the simulation contains. |

**Key attributes of the data**

- **Design:** three groups (23, 23, and 20 participants, 66 in total). Group (IPV-BI or non-IPV-BI) is assigned at clinical intake.
- **Demographics:** age, sex, gender, education, ethnicity, sexual orientation, height, weight, assailant relationship.
- **IPV-BI exposure history:** incident type (`incident_*`), head impact mechanisms (`head_mech_*`), strangulation signs (`strang_sign_*`), historic head impact and strangulation counts, and months since the first event.
- **Symptoms:** neurological (`neuro_sx_*`) and non-neurological (`non_neuro_sx_*`) symptoms, each with presence (No, Unknown, Yes), days to resolve, and a persistence flag.
- **Substance use and diagnoses:** `substance_*` frequency columns and `dx_*` psychiatric diagnosis columns.
- **Neuropsychological measures used in Milestone 4:** RAVLT total and delayed recall, TMT-A and TMT-B (T-scores), TOPF (standard score), COWA, SDMT oral and written (z-scores), WAIS PSI and WMI (standard scores), and the PCL-5, PHQ-9, and NSI totals. TOMM trials are included as a performance validity indicator. BISA and HELPS are used as exposure measures.
- **Biomarkers:** 222 `NPQ_<target>` columns from the NULISAseq Neuro 220 panel, stored as **log2-scale NPQ** (the Milestone 3 script halts if the values do not look log2-scale). The target MSLN is excluded from modeling in Milestones 3 and 4.
- **Plate quality control columns:** `nc_*`, `ic_*`, and `ipc_*` (negative control, internal control, and inter-plate control).
- **Missing data:** limited to the neuropsychological measures: missing completely at random by design, about 4% per variable. The biomarker panel has no missing values.

### Previous versions of the simulated data

Earlier versions of the simulated dataset were modified in later versions and **are not the most recent files**. Do not use them to run the current scripts.

| Version | Status |
|---|---|
| v1 to v4 | Did not pass the initial quality check. Not uploaded. |
| v5 | Run for the initial analysis, then replaced because the biomarker scores were not adequately distributed. |
| v6 | Did not include quality checks. Replaced by v7. |
| v7 | Current version. Use this for every script. |

If `IPV_cohort_simulated_v5.xlsx` or `IPV_cohort_simulated_v6.xlsx` are present in your copy of the project, they are previously simulated datasets and not the most recent files. Earlier versions of this README referred to v5 as the dataset to use; that is no longer correct.

### Other input

`Neuro_220_Panel.csv` is the panel annotation table (target, UniProt ID, hallmark of neurodegenerative disease, full target name). Milestone 3 reads it to attach hallmark categories to each protein and halts if any protein has no match.

### Derived data files

| File | Created by | Description |
|---|---|---|
| `ipv_cohort_clean.rds` | Milestone 1 | Cleaned and recoded cohort. Milestone 3 adds QC and adjusted p-value columns; Milestone 4 adds the `tomm_flag` validity column. |
| `pairwise_results.rds` | Milestone 3 | List of three results tables, one per group contrast, with limma and Mann-Whitney statistics. |
| `intake_summary_tables.xlsx` | Milestone 1 | Intake summary tables by study group. |
| `results_<contrast>_vN.csv` | Milestone 3 | Per-contrast biomarker results. |
| `milestone4_results.xlsx` | Milestone 4 | One workbook with a sheet per result set (see [Milestone 4](#milestone-4-neuropsychology-and-biomarker-association)). |

## Repository structure and file naming

The scripts expect this layout, with the data files in the project root and scripts in `Scripts/`:

NRSC_510/
├── README.md
├── NRSC_510.Rproj
├── IPV_cohort_simulated_v7.xlsx          # CURRENT simulated dataset
├── Neuro_220_Panel.csv                   # Panel annotation (hallmarks)
├── ipv_cohort_clean.rds                  # Cleaned cohort (created by Milestone 1)
├── pairwise_results.rds                  # Biomarker contrast results (created by Milestone 3)
├── IPV_cohort_simulated.xlsx             # Earlier simulation, not current
├── IPV_cohort_simulated_v5.xlsx          # Earlier simulation, not current
├── IPV_cohort_simulated_v6.xlsx          # Earlier simulation, not current
├── Scripts/
│   ├── cleaning_summary_v3.R                     # Milestone 1 (current)
│   ├── group_comparisons_milestone_2.R           # Milestone 2
│   ├── heatmap_volcano_milestone_3_v2.R          # Milestone 3 (current)
│   ├── volcano_plots_group_comparison.R          # Milestone 3 supplement (standalone)
│   ├── neuropsych_biomarker_milestone_4.R        # Milestone 4
│   ├── cleaning_summary.R, cleaning_summary_v2.R        # Earlier Milestone 1 versions
│   ├── heatmap_volcano_milestone_3.R                    # Earlier Milestone 3 version
│   └── class_demo_pca_script.r, quiz_october_*.R,       # Course exercise scripts,
│       test workshop 3.R, test_10.23.R,                 # not part of the pipeline
│       optional visualization code.R
└── Outputs/                              # All figures, tables, and exported results. Note that multiple versions of figures may exist while troubleshooting visualization occurred. Please use the highest number version for most accurate output
```

### File naming conventions

- **Contrast names** used in file names: `Acute_IPV_vs_NonIPV` (primary), `Chronic_vs_Acute_IPV` (secondary), `Chronic_IPV_vs_NonIPV` (tertiary).
- **`_vN` suffix (Milestone 3 and the volcano supplement):** a run label that is increased by hand each run so earlier outputs are not overwritten. For any figure family, the highest `_vN` is the most recent run. Files without a suffix are the earliest runs.
- **Milestone 3 outputs:**
  - `pval_raw_<contrast>_vN.png` and `pval_adj_<contrast>_vN.png`: p-value histograms (raw and BH-adjusted). `diagnostic_pval_*` files are from an earlier version of this step.
  - `volcano_primary_<contrast>_vN.png`: volcano plots at BH-adjusted p < 0.05.
  - `volcano_exploratory_<contrast>_vN.png`: volcano plots at the looser BH-adjusted p < 0.20 threshold.
  - `volcano_group_<contrast>_vN.png`: volcano plots colored by study group (supplement script).
  - `heatmap_A_*`, `heatmap_B_*`, `heatmap_C_*`, `heatmap_D_*`: four heatmap variants (see [Milestone 3](#milestone-3-biomarker-group-differences)).
  - `results_<contrast>_vN.csv`: result tables.
- **Milestone 4 outputs** have no run label and are overwritten on re-run: `violin_cognitive.png`, `violin_psychiatric.png`, `heatmap_biomarker_cognitive.png`, `heatmap_biomarker_psychiatric.png`, `heatmap_regression_all_proteins.png`, `volcano_regression_all.png`, `volcano_regression_<group>.png`, `scatter_targeted.png`, `scatter_NPQ<protein>_<outcome>.png` (for example `scatter_NPQCST3_WAISWMI.png`), `scatter_panel_<domain>.pdf`, and `milestone4_results.xlsx`.
- **Milestone 1 and 2 tables:** `Table_1_demo.html` to `Table_7_psych.html` (intake tables) and `Categorical_demo.html`, `mechanism_comparison.html`, `non_neuro_vs_neuro.html`, `substance_use_psych.html` (group comparisons).
- Files named `test*` and `Rplot*` are scratch outputs from development.
- Group colors are the same in every figure: Acute IPV-BI `#E41A1C`, Acute Non-IPV-BI `#377EB8`, Chronic IPV-BI `#984EA3`.

---

## Software and package requirements

**Language:** R, version 4.6.1 (2026-06-24, ucrt build) (https://www.r-project.org/)  
**Platform:** Windows 11 x64  
**IDE:** RStudio, version 2026.8.2.200 (https://posit.co/products/open-source/rstudio/)

Versions below are those used to run the scripts. To reproduce your own record, run all scripts in one session and then run `sessionInfo()` (R and packages) and `RStudio.Version()$version` (RStudio).

| Package | Version | Used in | Purpose | Link |
|---|---|---|---|---|
| tidyverse | 2.0.0 | All | Data import, wrangling, piping, plotting | https://www.tidyverse.org/ |
| dplyr | 1.2.1 | All (via tidyverse) | Data manipulation | https://dplyr.tidyverse.org/ |
| tidyr | 1.3.2 | All (via tidyverse) | Reshaping data | https://tidyr.tidyverse.org/ |
| readr | 2.2.0 | All (via tidyverse) | Reading and writing delimited files | https://readr.tidyverse.org/ |
| purrr | 1.2.2 | All (via tidyverse) | Iterating functions over lists | https://purrr.tidyverse.org/ |
| stringr | 1.6.0 | All (via tidyverse) | String handling | https://stringr.tidyverse.org/ |
| tibble | 3.3.1 | All (via tidyverse) | Data frames | https://tibble.tidyverse.org/ |
| forcats | 1.0.1 | All (via tidyverse) | Factor handling | https://forcats.tidyverse.org/ |
| lubridate | 1.9.5 | All (via tidyverse) | Dates | https://lubridate.tidyverse.org/ |
| ggplot2 | 4.0.3 | All (via tidyverse) | Visualization | https://ggplot2.tidyverse.org/ |
| readxl | 1.6.0 | 1 | Reading the Excel input file | https://readxl.tidyverse.org/ |
| gtsummary | 2.6.1 | 1 | Intake summary tables | https://www.danieldsjoberg.com/gtsummary/ |
| naniar | 1.1.0 | 1 | Missing data visualization | https://naniar.njtierney.com/ |
| huxtable | 6.0.2 | 1 | Table export | https://hughjonesd.github.io/huxtable/ |
| openxlsx | 4.2.9 | 1, 4 | Excel workbook export | https://ycphs.github.io/openxlsx/ |
| gt | 1.3.0 | 2 | Formatted comparison tables | https://gt.rstudio.com/ |
| emmeans | 2.0.4 | 2 | Covariate-adjusted marginal means and pairwise contrasts | https://rvlenth.github.io/emmeans/ |
| MASS | 7.3-65 | 2 | Ordinal logistic regression (`polr`) | https://cran.r-project.org/package=MASS |
| broom | 1.0.13 | 2 | Tidying model output | https://broom.tidymodels.org/ |
| ggpubr | 1.0.0 | 2, 4 | Significance brackets, multi-panel figures | https://rpkgs.datanovia.com/ggpubr/ |
| limma | 3.68.5 | 3 | Linear models for high-dimensional biomarker data | https://bioconductor.org/packages/limma/ |
| pheatmap | 1.0.13 | 3, 4 | Heatmaps | https://cran.r-project.org/package=pheatmap |
| ggrepel | 0.9.8 | 3, 4 | Non-overlapping point labels | https://ggrepel.slowkow.com/ |
| RColorBrewer | 1.1-3 | 3, 4 | Color palettes | https://cran.r-project.org/package=RColorBrewer |

## Steps to replicate

1. **Install** R and RStudio (versions above).
2. **Get the project.** Clone the repository (`git clone https://github.com/court-mccabe/NRSC_510.git`) or download it as a ZIP, so the folder layout matches [the structure above](#repository-structure-and-file-naming).
3. **Add the data.** Put `IPV_cohort_simulated_v7.xlsx` and `Neuro_220_Panel.csv` in the project root. The Dataverse linked above also holds the outcome data and the project proposal. Use v7 only.
4. **Open `NRSC_510.Rproj`** in RStudio so the working directory is the project root. Each script also sets the path near the top (`setwd("C:/NRSC_510")` and an `output_dir` line). Change those lines to your own project path, or remove `setwd()` when using the `.Rproj`.
5. **Install the packages:**

   ```r
   install.packages(c("tidyverse", "readxl", "gtsummary", "naniar", "huxtable",
                      "openxlsx", "gt", "emmeans", "MASS", "broom", "ggpubr",
                      "ggrepel", "pheatmap", "RColorBrewer"))
   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   BiocManager::install("limma")
   ```

   To match a listed version exactly, use `remotes::install_version("dplyr", "1.2.1")`.
6. **Run the scripts in the [run order](#run-order-at-a-glance)**, restarting R between scripts (see the [known issues for Milestone 2](#known-issues-milestone-2) before running it). Source each script top to bottom (Source button, or `source("Scripts/<name>.R")`).
7. **Check the checkpoints** for each milestone (tables in the walkthroughs below). A script that works prints its checkpoint tables, raises no unexpected warnings, and writes the files listed under "Produces".
8. **Before re-running Milestone 3 or the volcano supplement,** increase `run_label` so earlier figures are kept.
9. **Record your environment** by running `sessionInfo()` after the last script.

The v7 dataset is itself a fully synthetic test dataset: it contains no real participants, every script can be run on it end to end, and the checkpoints below confirm that the structure and scale of the data are what the code expects.

**How checkpoints are labelled:** **HALT** means the script stops with an error message. **WARN** means it continues and prints a warning. **PRINT** means a small table or message is shown in the console for you to confirm.

## Reproducibility notes

- **Inputs are fixed.** All analyses use the same simulated file (`IPV_cohort_simulated_v7.xlsx`), so every user starts from identical data.
- **Tests are deterministic with one exception.** Milestone 2, Step 12 computes Fisher's exact p-values by Monte Carlo simulation (`simulate.p.value = TRUE`). No random seed is set in the scripts, so these p-values can differ slightly between runs. Run `set.seed()` with a fixed number before Step 12 if you need identical values.
- **Figure layout can vary.** Point jitter (Milestone 4 violin plots) and label placement (`ggrepel` in the volcano plots) are randomized, so figures can differ slightly in layout between runs while showing the same data.
- **Environment.** The scripts were run on R 4.6.1 under Windows 11 with the package versions listed above. Other versions can change formatting and warnings. The full session record is in [Session information](#session-information).

### Milestone 1: Cleaning and intake summaries

**Script:** `Scripts/cleaning_summary_v3.R`  
**Input:** `IPV_cohort_simulated_v7.xlsx`  
**Purpose:** Import the two cohort sheets, recode variables, and summarize intake characteristics by study group. This milestone is descriptive: no group significance tests are run.

**Produces:** `ipv_cohort_clean.rds`, `intake_summary_tables.xlsx`, `Table_1_demo.html` to `Table_7_psych.html`, and a missing-data plot of the neuropsychological variables (`NP_missing_data.png`).

### Milestone 2: Group comparisons

**Script:** `Scripts/group_comparisons_milestone_2.R`  
**Input:** `ipv_cohort_clean.rds`  
**Purpose:** Compare demographics, IPV-BI exposure, injury mechanisms, symptoms, substance use, and psychiatric diagnoses across the three study groups.

All multiple-comparison corrections use Benjamini-Hochberg false discovery rate (BH-FDR), chosen because these are exploratory, hypothesis-generating comparisons with many simultaneous tests. Models in Steps 13 to 15 adjust for age and sex. Demographics in Steps 11 and 12 are unadjusted because they are the covariates used later.

**Checkpoints**

| Step | Type | What to confirm |
|---|---|---|
| Start | Setup | `select <- dplyr::select` is set so `MASS` loading does not mask `dplyr::select()`. |
| 12 | PRINT | The results table lists which test (chi-square or Fisher's exact) was used for each variable. Fisher's exact is expected for most variables at this sample size. |
| 12 | WARN | "Chi-squared approximation may be incorrect" warnings are expected. The script calls `chisq.test()` only to read the expected cell counts before choosing Fisher's exact test. |
| 13 to 15 | WARN | "glm.fit: fitted probabilities numerically 0 or 1 occurred" can appear in the logistic regression steps for columns with very few events in a group. This reflects the documented low-event-count limitation and does not stop the script. |
| 13 to 15 | Built-in | Each model call is wrapped in `tryCatch()`. A column whose model fails returns an `error` row instead of stopping the script, and these rows are removed before correction (`filter(!is.na(p.value))`). |
| 14 | PRINT | `sapply(class)` on the symptom presence columns shows ordered factors (No < Unknown < Yes) before the ordinal models are run. |

**Produces:** formatted `gt` tables (HTML copies in `Outputs/`: `Categorical_demo.html`, `mechanism_comparison.html`, `non_neuro_vs_neuro.html`, `substance_use_psych.html`) and `Age_Comparison_R.png`.

#### Known issues (Milestone 2)

The script is kept as originally written. These issues are documented so they are easy to work around.

1. **`ordinal_emmeans_by_group()` is not defined in this script.** Steps 14 and 15 call it for symptom presence and substance use. The function was defined in an earlier version of Step 13, and Step 13 now defines a binary logistic function instead. The script works only if the function already exists in your R session; in a fresh session, it stops at Step 14 with "object 'ordinal_emmeans_by_group' not found". **Workaround:** run this definition once before Step 14.

   ```r
   ordinal_emmeans_by_group <- function(outcome, data, covariates = c("age_years", "sex")) {
     formula <- as.formula(
       paste(outcome, "~", paste(c("study_group", covariates), collapse = " + "))
     )
     tryCatch({
       model     <- polr(formula, data = data, Hess = TRUE)
       emm       <- emmeans(model, pairwise ~ study_group, mode = "linear.predictor")
       contrasts <- as.data.frame(emm$contrasts) %>%
         mutate(outcome = outcome)
       contrasts
     }, error = function(e) {
       data.frame(outcome = outcome, error = conditionMessage(e))
     })
   }
   ```

2. **One diagnostic line in Step 14 can stop the script.** The block commented `# check to see what's failing` ends with `pull(error)`. It only works when some models fail. When every model succeeds, there is no `error` column, and it stops with "object 'error' not found". It prints information only, so skip that block and continue with the next one; the results are not affected.
3. **Two comments no longer match the code.** The Step 13 header describes ordinal regression (`polr`) for the exposure columns, but the code fits binary logistic regression (`glm`, binomial). The Step 15 comment says the ordinal function comes from Section 13, but it is called from Step 14 onward.

### Milestone 3: Biomarker group differences

**Script:** `Scripts/heatmap_volcano_milestone_3_v2.R` (Steps 16 to 35)  
**Inputs:** `ipv_cohort_clean.rds`, `Neuro_220_Panel.csv`  
**Purpose:** Test every protein for differences between study groups, with quality control, a nonparametric cross-check, and visualization. Biomarker values are log2-scaled NPQ and are used as they are.

Three contrasts are tested: (1) Acute IPV-BI vs Acute Non-IPV-BI (primary), (2) Chronic IPV-BI vs Acute IPV-BI (secondary), and (3) Chronic IPV-BI vs Acute Non-IPV-BI (tertiary).

**Checkpoints**

| Step | Type | Condition | Expected on a successful run |
|---|---|---|---|
| 18 | HALT | NPQ values do not look log2-scale (mean, spread, or minimum outside the expected range). The message says to rebuild `ipv_cohort_clean.rds` from the v7 workbook. | Passes silently. |
| 19 | WARN | Negative control (NC), internal control (IC), or inter-plate control (IPC) thresholds exceeded; the IC warning names the failing samples. Also warns if the QC columns are missing. | Any warning describes the check that failed; the script continues. |
| 21 | HALT | A protein has no hallmark match in `Neuro_220_Panel.csv`. | Passes silently. |
| 24 | HALT | BH-adjusted p-values for the primary contrast are almost all 1.0, which indicates a data format problem. | Passes silently. |

**Produces:** `pval_*`, `volcano_primary_*`, `volcano_exploratory_*`, and `heatmap_A` to `heatmap_D` figures, `results_<contrast>_vN.csv` for each contrast, `pairwise_results.rds`, and the updated `ipv_cohort_clean.rds`.

#### Milestone 3 supplement: volcano plots colored by study group

**Script:** `Scripts/volcano_plots_group_comparison.R` (standalone, optional)  
**Input:** `pairwise_results.rds`
Draws one volcano plot per contrast. Significant proteins are colored with the color of the group in which they are higher; non-significant proteins are grey. Colors match the group colors above. Increase `run_label` before each run. **Produces:** `volcano_group_<contrast>_vN.png`.

### Milestone 4: Neuropsychology and biomarker association

**Script:** `Scripts/neuropsych_biomarker_milestone_4.R` (Steps 41 to 55)  
**Inputs:** `ipv_cohort_clean.rds`, `pairwise_results.rds`  
**Purpose:** Compare neuropsychological scores across groups, then relate the biomarker panel to those scores. 

Scores are put on one scale for comparison: T-scores as (x - 50) / 10, standard scores as (x - 100) / 15, existing z-scores unchanged, and raw scores as within-sample z-scores. The 13 outcomes are split into a cognitive domain (10 scores) and a psychiatric domain (PCL-5, PHQ-9, NSI).

**Checkpoints**

| Step | Type | What to confirm |
|---|---|---|
| 42 | PRINT | A TOMM message appears only if a participant scored below the validity threshold. No message means none were flagged. |
| 44 | PRINT | A two-row summary (cognitive, psychiatric) showing the number of comparisons at FDR < 0.05 and the number with a medium or larger effect size. |
| 48 | PRINT | A summary table with one row per group and family showing the number of correlations at FDR < 0.05 and FDR < 0.20. |
| 49 | PRINT | The ten-row top correlations table with rank, protein, outcome, method, r, n, p, and FDR. |
| 50 | PRINT | A message such as "No proteins at FDR < 0.20 for psychiatric so heatmap skipped" appears when a heatmap has nothing to show. This is expected behaviour, not an error. |
| 51 | PRINT | One-row summary: number of models, and the number at FDR < 0.05 and < 0.20. The model count normally equals the number of proteins times 13 outcomes. |
| 53 | HALT / PRINT | Halts if `pairwise_results.rds` is missing (run Milestone 3 first). Prints the three selected proteins with their log fold change and limma adjusted p-value, then the nine targeted results. |
| 55 | PRINT | A closing message gives the output folder. |

**Produces (Outputs/):** `violin_cognitive.png`, `violin_psychiatric.png`, `heatmap_biomarker_cognitive.png` (and `_psychiatric` when proteins qualify), `heatmap_regression_all_proteins.png`, `volcano_regression_all.png`, `volcano_regression_AcuteIPVBI.png`, `volcano_regression_AcuteNonIPVBI.png`, `volcano_regression_ChronicIPVBI.png`, `scatter_targeted.png`, one `scatter_NPQ<protein>_<outcome>.png` per plotted pair with `scatter_panel_<domain>.pdf`, and `milestone4_results.xlsx`.

**`milestone4_results.xlsx` sheets:** Notes, Top 10 Correlations, Cognitive Group Tests, Psychiatric Group Tests, Biomarker All Neuropsych, Biomarker Cognitive, Biomarker Psychiatric, Biomarker Exposure, Within Group Correlations, Regression Cognitive, Regression Psychiatric, Within Group Regression, Targeted Regression. The Notes sheet explains how to read each one. Every sheet is sorted by FDR with a frozen header and filters.

## Documented limitations

- **Small samples and low event counts (Milestone 2).** Some persistent symptom columns have very few events per group, so non-significant results there likely reflect low power, not a confirmed absence of a difference. Firth's penalized logistic regression (`logistf`) is a more robust option if rare events persist at larger sample sizes.
- **Future extension (Milestone 2).** At larger sample sizes, revisit the correction method: Tukey's HSD for the demographic comparisons, with BH-FDR remaining appropriate for the larger exploratory sets. Reassess whether age and sex are enough covariates.
- **Covariate adjustment (Milestone 3).** The Mann-Whitney cross-check cannot adjust for age and sex; limma can.
- **Milestone 4.**
  1. Power is low for 221 proteins across 13 outcomes at n = 66, so the all-protein analyses are hypothesis generating.
  2. The all-protein correlations and regressions pool the three groups and do not adjust for study group; the targeted step shows models with and without it.
  3. Linear regression assumes straight-line relationships, roughly normal residuals, and no single dominating point; residuals were not checked model by model.
  4. Several outcomes are related to each other, so the 13 outcomes are not independent tests.
  5. PCL-5, PHQ-9, and NSI are summed questionnaire scores treated as continuous, and raw scores are standardized within this sample, not against population norms.
  6. Education, kidney function, and other possible confounders are not modeled.
  7. Data come from a single time point, so associations cannot show direction of effect.
  8. The data are simulated, so results are not biological findings.

## Code documentation conventions

- Variable names describe the test or purpose they serve (for example `age_height_weight_anova_cols`, `ipv_exposure_mechanism_contrasts`) rather than generic names.
- Reusable modeling functions (for example `binary_emmeans_by_group`, `linear_emmeans_by_group`, `logistic_emmeans_by_group` in Milestone 2, and `run_correlations`, `run_lm_per_protein`, `adjusted_line` in Milestone 4) are defined once and reused instead of duplicated.
- In-line comments explain why an approach was chosen, particularly for clinically important decisions (for example, why "Unknown" is kept as its own category).
- Console output is limited to checkpoints and tables that feed later steps.
- Future extension notes are embedded at the relevant place in each script.

## Session information

Output of `sessionInfo()` from the R session used for these scripts (Milestones 2 and 4 packages were loaded).

<details>
<summary>Show full session information</summary>

```
R version 4.6.1 (2026-06-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C
[5] LC_TIME=English_United States.utf8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
 [1] gt_1.3.0           broom_1.0.13       MASS_7.3-65        emmeans_2.0.4
 [5] ggrepel_0.9.8      openxlsx_4.2.9     RColorBrewer_1.1-3 ggpubr_1.0.0
 [9] lubridate_1.9.5    forcats_1.0.1      stringr_1.6.0      dplyr_1.2.1
[13] purrr_1.2.2        readr_2.2.0        tidyr_1.3.2        tibble_3.3.1
[17] ggplot2_4.0.3      tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       rstatix_1.1.0      lattice_0.22-9     tzdb_0.5.0
 [5] vctrs_0.7.3        tools_4.6.1        generics_0.1.4     pkgconfig_2.0.3
 [9] pheatmap_1.0.13    Matrix_1.7-5       S7_0.2.2           lifecycle_1.0.5
[13] compiler_4.6.1     farver_2.1.2       textshaping_1.0.5  carData_3.0-6
[17] sass_0.4.10        htmltools_0.5.9    Formula_1.2-6      pillar_1.11.1
[21] car_3.1-5          abind_1.4-8        nlme_3.1-169       digest_0.6.39
[25] tidyselect_1.2.1   zip_3.0.2          mvtnorm_1.4-2      stringi_1.8.9
[29] labeling_0.4.3     splines_4.6.1      fastmap_1.2.0      cowplot_1.2.0
[33] grid_4.6.1         cli_3.6.6          magrittr_2.0.5     utf8_1.2.6
[37] withr_3.0.3        scales_1.4.0       backports_1.5.1    timechange_0.4.0
[41] estimability_2.0.0 gridExtra_2.3.1    ggsignif_0.6.4     ragg_1.5.2
[45] hms_1.1.4          mgcv_1.9-4         rlang_1.3.0        Rcpp_1.1.2
[49] xtable_1.8-8       glue_1.8.1         xml2_1.6.0         rstudioapi_0.19.0
[53] R6_2.6.1           fs_2.1.0           systemfonts_1.3.2
```

</details>

---

## Acknowledgments

I am grateful to the Wellington Lab for its guidance and support throughout this project.

Thank you to Margo, the peer support tutor for NRSC 510, and to the course instructors and guest lecturers for the material covered throughout the course.

---

## References

- GitHub Docs. *About READMEs.* https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-readmes
- Advocates, Academics, Survivors and Clinicians to END Intimate Partner Violence (ASCEND-IPV) initiative: a prospective observational case-control study protocol to identify plasma biomarkers of intimate partner violence (IPV)-caused brain injury (BI). (2025). *BMJ Open.* doi:10.1136/bmjopen-2024-098025. Reference study design for the simulated cohort.
- Makovec Knight et al. (2025). *J Neurotrauma.* doi:10.1177/08977151251388070. Source of the age distribution, the injury burden split, and the direction of the cognitive effects.
- Cimino et al. (2019). *J Aggress Maltreat Trauma.* Source of the head impact and strangulation mechanism mix (as recorded in the `Parameters_Sources` sheet).
- Herbenick et al. (2025). Source of the sexual orientation distribution and the strangulation mechanism categories (as recorded in the `Parameters_Sources` sheet).
- Mehr et al. (2023). Source of the ordering of substance use rates (as recorded in the `Parameters_Sources` sheet).
- Esopenko et al. (2024). Source of the direction of the built-in biomarker shifts (as recorded in the `Parameters_Sources` sheet).
- Alamar Biosciences. *NULISAseq Neuro 220 panel.* https://alamarbio.com/products-and-services/nulisaseq-neuro-220-panel/
- Alamar Biosciences. *NULISAseqR User Guide* (Chapter 2, quality control). https://nulisaseqr.alamarbio.com/user-guide/quality-control.html. Source of the plate quality control thresholds used in the simulated quality control columns and in Milestone 3, Step 19.
- Ritchie et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7), e47.
- Benjamini & Hochberg (1995). Controlling the false discovery rate. *Journal of the Royal Statistical Society B*, 57(1), 289-300.
- Cohen (1988). *Statistical Power Analysis for the Behavioral Sciences* (2nd ed.). Effect size benchmarks used in Milestone 4.
- Ferreira & Castro (2022). Linear and logistic regression models: when to use and how to interpret them. *J Bras Pneumol.* PMC9747134.

## Citation

If reusing this code or its generating logic, please cite:

McCabe, C. (2026). *Identifying Biomarkers and Neurocognitive Function to Improve Diagnostic Accuracy for Brain Injury in Survivors of Intimate Partner Violence* [Source code]. GitHub. https://github.com/court-mccabe/NRSC_510

This code operates on **simulated data** and was developed for coursework (NRSC 510) under the guidance of the Wellington Lab.
