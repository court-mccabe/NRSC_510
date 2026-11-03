Identifying Biomarkers and Neurocognitive Function to Improve Diagnostic Accuracy for Brain Injury in Survivors of Intimate Partner Violence

Author: Courtney McCabe 
Course: NRSC 510 Affiliated 
Lab: Wellington Lab 
GitHub Repository: https://github.com/court-mccabe/NRSC_510

Project Overview

Intimate partner violence (IPV) is a leading cause of traumatic brain injury (TBI), yet IPV-caused brain injury (IPV-BI) remains underrepresented in TBI and Alzheimer's Disease and Related Dementias (ADRD) research. This project investigates whether plasma biomarker profiles and neuropsychological performance indicate elevated risk of long-term neurological and psychological impairment in IPV-BI survivors, and whether this risk differs between acute and chronic timepoints.

This Dataverse repository supports a course project (NRSC 510) analyzing a simulated dataset that mirrors the structure and variables of the ASCEND-IPV Biomarkers study. The real cohort data is not available for public use, as it is being analyzed separately by another member of the Wellington Lab. The simulated dataset was built using Alamar datasets, normative neuropsychological scores, and predicted cohort outcomes to allow independent analysis and skill-building without compromising participant confidentiality.

Disclaimer: This is simulated data. It was generated to mirror the structure, variable types, and approximate distributions of a real prospective cohort study (ASCEND-IPV), but it does not represent real participants or real clinical outcomes.

No clinical interpretations, diagnostic conclusions, or treatment decisions should be made based on this dataset. It is intended solely for methodological development, statistical training, and reproducibility demonstration purposes.

Repository Contents
File / Folder	Description

Project Proposal: Identifying Biomarkers and Neurocognitive Function to Improve Diagnostic Accuracy for Brain Injury in Survivors of Intimate Partner Violence --	Original NRSC 510 project proposal outlining background, research question, dataset description, and analysis plan
Update_Presentation_1:	Project update presentation summarizing progress through Weeks 1–4 of the milestone timeline
Simulated_Dataset_NRSC_510.xlsx:	Raw simulated dataset prior to cleaning
IPV_RSD_clean: Cleaned and recoded dataset following Section 5 variable recoding (logical, ordered factor, and resolution column transformations)
intake_summary_tables.xlsx	Exported intake summary tables (Section 7) — demographics, IPV-BI exposure history, symptomatology, substance use, and psychiatric diagnoses by study group
README_Dataverse.md	This file

The code used to generate ipv_cohort_clean.rds, the summary tables, and the statistical comparisons is maintained separately on GitHub in scripts/cleaning_summary.R (see link above) and is not duplicated in this Dataverse. The GitHub repository contains its own README documenting the code specifically.

Key Attributes of the Data
Sample size: 66 simulated participants
Acute (≤30 days since incident): n = 46, split 50/50 into IPV-BI (n = 23) and non-IPV-BI (n = 23)
Chronic (>6 months since incident): n = 20, IPV-BI only
Study groups: Acute/Subacute IPV-BI, Acute/Subacute non-IPV-BI, Chronic IPV-BI
Variable categories:
Demographics (age, sex, gender, ethnicity, education, height, weight, sexual orientation)
Days since injury / incident timing
Historic and cumulative lifetime IPV-BI exposure (strangulation, head/face impacts)
Assailant relationship to participant
Substance use (5-level frequency scale)
Strangulation mechanism and signs
Head impact mechanism and neurological signs
Neurological and non-neurological symptomatology at intake, including symptom resolution status (resolved/persistent, days to resolution)
Psychiatric diagnoses and related clinical variables (suicide ideation, substance coping, sleep disturbance)
Neuropsychological test scores (RAVLT, TMT-A/B, TOPF, COWA, SDMT, WAIS indices, TOMM, PCL-5, PHQ-9, NSI, HELPS)
Alamar NULISAseq Neuro 220 biomarker z-scores
Missingness: Approximately 3–4% of cognitive data are missing by design, modeled on realistic data collection patterns and parameterized using published normative work. 

Dataset Access

This dataset is stored in this Dataverse repository. No special access permissions are required to use the simulated dataset for research or educational purposes — see the License section below.

Note: This simulated dataset is distinct from the real ASCEND-IPV cohort data, which is restricted and analyzed separately under lab-specific data use agreements. Requests for access to the real cohort data should be directed to the Wellington Lab directly and are outside the scope of this Dataverse.

License

License: CC BY 4.0 (Attribution 4.0 International)

Justification

This dataset is shared under a CC BY 4.0 license, which permits others to copy, redistribute, and adapt the data for any purpose, including commercial use, provided appropriate credit is given.

CC BY 4.0 was selected over a more permissive license such as CC0 for two reasons:

Attribution requirement: Although this is simulated data, it was constructed to reflect the structure and design of ongoing, original research conducted by the Wellington Lab. Requiring attribution ensures the lab and this project receive appropriate credit if the simulated dataset or its generating methodology is reused, cited, or built upon by others.
Academic norms: CC BY 4.0 is a standard and widely recognized license in academic data sharing, balancing open reuse with proper scholarly credit, which is particularly important for a dataset tied to an active area of IPV/TBI research.

When reusing this dataset, please credit: Courtney McCabe and the Wellington Lab, and clearly note that the data is simulated and not suitable for clinical interpretation.