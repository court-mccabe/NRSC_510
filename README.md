NRSC_510: IPV-BI Biomarker and Neurocognitive Analysis
Author: Courtney McCabe Course: NRSC 510 
Affiliated Lab: Wellington Lab 
Related Dataverse: https://demo.borealisdata.ca/dataverse/courtmccabe/

Disclaimer: This repository analyzes a simulated dataset built to mirror the structure of the ASCEND-IPV Biomarkers study. It does not contain real participant data. No clinical interpretations should be drawn from any results produced by this code.

Project Overview

This project investigates whether plasma biomarker profiles and neuropsychological performance indicate elevated risk of long-term neurological and psychological impairment in survivors of intimate partner violence–caused brain injury (IPV-BI), and whether this risk differs between acute and chronic timepoints.

The analysis compares three study groups:

Acute/Subacute IPV-BI
Acute/Subacute non-IPV-BI
Chronic IPV-BI

Project Goals

Build a reproducible tidyverse-based data cleaning and analysis workflow in R
Summarize intake characteristics across neurological and non-neurological symptom domains
Compare demographics, IPV-BI exposure history, and injury mechanisms across study groups
Apply covariate-adjusted statistical models (ordinal logistic, linear, logistic regression via emmeans) with multiple comparison correction
Extend into high-dimensional biomarker analysis (limma) and biomarker–cognition correlation in later project stages

Repository Structure
NRSC_510/
├── README.md                          # This file
├── data/
│   ├── IPV_cohort_simulated_v5_xlsx       # Raw simulated dataset
│   └── ipv_cohort_clean.rds               # Cleaned dataset 
├── scripts/
│   └── cleaning_summary.R              # summary variables and organize data
└── outputs/
    ├── intake_summary_tables.xlsx     # Section 7 exported tables
    └── group_comparison_results.xlsx  # Sections 11–15 exported results


Software and Package Requirements

Software/Package	Version Used	Purpose	Link
R	(4.6.1)	Base language	https://www.r-project.org/
tidyverse	(2.0.0)	Data import, wrangling, piping	https://www.tidyverse.org/
dplyr	(e.g. 1.2.1)	Data manipulation (select, mutate, filter, across)	https://dplyr.tidyverse.org/
readxl	(1.5.0.1)	Reading Excel input files	https://readxl.tidyverse.org/
ggplot2	(4.0.3)	Data visualization	https://ggplot2.tidyverse.org/
gtsummary	(2.6.1)	Demographic and intake summary tables	https://www.danieldsjoberg.com/gtsummary/
naniar	(1.1.0)	Missing data visualization and handling	https://naniar.njtierney.com/
emmeans	Covariate-adjusted estimated marginal means and pairwise contrasts	https://rvlenth.github.io/emmeans/
MASS	Ordinal logistic regression (polr)	https://cran.r-project.org/package=MASS
broom	(1.0.13)	Tidying model output	https://broom.tidymodels.org/
ggpubr	(1.0.0)	Significance brackets on boxplots	https://rpkgs.datanovia.com/ggpubr/
huxtable	(6.0.2)	Table export to Excel	https://hughjonesd.github.io/huxtable/
openxlsx	(4.2.9)	Excel workbook export	https://ycphs.github.io/openxlsx/
limma	(planned — Weeks 5–6)	High-dimensional biomarker modeling with FDR correction	https://bioconductor.org/packages/limma/
ComplexHeatmap	(planned — Weeks 5–6)	Biomarker heatmap visualization	https://bioconductor.org/packages/ComplexHeatmap/
EnhancedVolcano	(planned — Weeks 5–6)	Biomarker volcano plot visualization	https://bioconductor.org/packages/EnhancedVolcano/

Script Section Overview

Code is fully annotated in-script; this is a broad overview of what each section does and the key functions it relies on.

Section	Purpose	Key functions
1–4	Data import and setup	read_csv(), bind_rows()
5	Variable recoding: logical, ordered factor, and resolution-column splits	mutate(), across(), case_when(), factor(), str_remove()
7	Intake summary tables by study group	tbl_summary(), add_overall(), bold_labels(), modify_caption()
8	Export summary tables to Excel	createWorkbook(), writeData(), saveWorkbook()
11	Continuous demographic comparisons (age, height, weight)	aov(), emmeans(), p.adjust(method = "BH")
12	Categorical demographic comparisons	table(), chisq.test(), fisher.test()
13	IPV-BI exposure mechanism comparisons	polr(), emmeans()
14	Neuro vs. non-neuro symptom comparisons (presence, days to resolve, persistence)	polr(), lm(), glm(family = binomial), emmeans()
15	Substance use and psychiatric diagnosis comparisons	Reuses model functions from Sections 13–14

Sections 11–15 adjust for age and sex and use BH-FDR correction throughout. In-script notes flag documented limitations (e.g. low event rates in some persistent symptom columns) and future-extension points (e.g. revisiting correction method at larger sample sizes).

Planned: Weeks 5–8
limma-based biomarker modeling by group and phase
Biomarker–cognition correlation
Final report write-up

Code Documentation Conventions

Variable names are descriptive of the statistical test or purpose they serve (e.g. age_height_weight_anova_cols, ipv_exposure_mechanism_contrasts) rather than generic names, to support readability when revisiting the script later.
Reusable modeling functions (ordinal_emmeans_by_group, linear_emmeans_by_group, logistic_emmeans_by_group) are defined once and reused across sections rather than duplicated.
In-line comments explain not just what each line does but why that approach was chosen, particularly around clinically important decisions (e.g. why "Unknown" is preserved as its own category rather than collapsed to NA for certain variables).

Citation

If reusing this code or its generating logic, please cite:

McCabe, C. (2026). Identifying Biomarkers and Neurocognitive Function to Improve Diagnostic Accuracy for Brain Injury in Survivors of Intimate Partner Violence [Source code]. GitHub. https://github.com/court-mccabe/NRSC_510

Please also note that this code operates on simulated data and was developed for coursework (NRSC 510) under the guidance of the Wellington Lab.
Future extension notes are embedded directly in the script at the relevant section, flagging what to revisit if the dataset or sample size changes.
