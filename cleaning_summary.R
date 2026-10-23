# Week 1–2: Data Cleaning and Intake Summary
# Goal: Import the simulated ASCEND-IPV dataset, clean and recode variables,
# assess missingness, and produce an intake summary table describing
# recent incident characteristics, demographics, and substance use by group.

# Groups:
#   - Acute IPV-BI       (≤30 days post-incident, intimate partner violence-caused brain injury [IPV-BI])
#   - Acute Non-IPV-BI   (≤30 days post-incident, no IPV-BI)
#   - Chronic IPV-BI     (>6 months post-incident, IPV-BI)
# Note: in a larger cohort, chronic IPV-BI would be compared to healthy normative population

# Step 1: Load packages

library(tidyverse)   # data wrangling (dplyr, tidyr, stringr) and plotting
library(readxl)      # read .xlsx files (tidyverse-adjacent, loaded separately)
library(naniar)      # missing data assessment and visualization
library(gtsummary)   # formatted summary tables
install.packages("huxtable")
install.packages("openxlsx")  # huxtable uses this as its xlsx engine
library(huxtable)
library(openxlsx)

# Step 2: Import data

# The dataset has two sheets: one for the acute/subacute phase and one for
# chronic. We read them separately then stack them into one dataframe so all
# downstream cleaning and summaries operate on the full cohort at once.

# Update file path
file_path <- "C:/NRSC_510/IPV_cohort_simulated_v5.xlsx"

acute_raw   <- read_excel(file_path, sheet = "Acute_Subacute")
chronic_raw <- read_excel(file_path, sheet = "Chronic")

# bind_rows() stacks two dataframes that share the same column names.
# for example, age_years will all be stacked in the same column not separated by tabs
# Any column present in one sheet but not the other will be filled with NA.
# This prevents any missing data from throwing an error in the code 
cohort_raw <- bind_rows(acute_raw, chronic_raw)

# Quick check: confirm expected row counts (46 acute + 20 chronic = 66 total)
# Use cat function -- Concatenate and Print
cat("Total rows:", nrow(cohort_raw), "\n")
cat("Columns:", ncol(cohort_raw), "\n")

# Step 3: Standardize "N/A" strings to true NA 

# Identify which of the true_missing_cols are character type
# (these may contain "N/A" strings that need converting)
# Numeric columns already have true NA where data is missing — no conversion needed

neuro_test_cols <- c(
  "RAVLT_total", "RAVLT_delayed_recall",
  "TMT_A_Tscore", "TMT_B_Tscore",
  "TOPF_standard_score", "COWA_total_words",
  "SDMT_oral_z", "SDMT_written_z",
  "WAIS_PSI", "WAIS_WMI",
  "TOMM_trial1", "TOMM_trial2",
  "PCL5_total", "PHQ9_total", "NSI_total",
  "BISA_BI_count", "HELPS_total"
)

char_missing_cols <- neuro_test_cols[
  sapply(cohort_raw[neuro_test_cols], is.character)
]

# Apply na_if() only to character columns within true_missing_cols
cohort <- cohort_raw %>%
  mutate(across(all_of(char_missing_cols), ~ na_if(.x, "N/A")))

# Should show true NA count in neuro tests, not zero
miss_var_summary(cohort %>% select(all_of(neuro_test_cols)))

# Step 4: Create a clean group label 

# The dataset has separate 'phase' and 'group' columns. We combine them into
# one three-level factor for use in all summaries and plots.
# use case_when fxn to assign new study_group column if both phase and group are 
# true 
# factor() with explicit levels, easier visual display of acute groups next to each
# other and then chronic group

cohort <- cohort %>%
  mutate(
    study_group = case_when(
      phase == "Acute/subacute" & group == "IPV-BI"     ~ "Acute IPV-BI",
      phase == "Acute/subacute" & group == "Non-IPV-BI" ~ "Acute Non-IPV-BI",
      phase == "Chronic"        & group == "IPV-BI"     ~ "Chronic IPV-BI"
    ),
    study_group = factor(study_group,
                         levels = c("Acute IPV-BI",
                                    "Acute Non-IPV-BI",
                                    "Chronic IPV-BI"))
  )

# Confirm group sizes match the study design
table(cohort$study_group)

# Step 5: Recode variables
# Raw data arrives with Yes/No strings, mixed-type resolution columns, and
# frequency scales stored as character. This section converts each column type
# to its appropriate R class before any summaries or analyses are run.
#
#
# 5a. Yes/No columns → logical (TRUE/FALSE)
# Applied only to fields where Unknown is not a clinically meaningful response:
# incident type and psychiatric diagnoses. A participant who does not know their
# anxiety diagnosis is functionally equivalent to missing data for that field.
# For dx_ columns specifically, Unknown was not offered as a response option on
# the intake form — participants could only select Yes or No. Any NA in these
# columns reflects true missing data, not patient uncertainty.
#
# Excluded from this block:
#   substance_     → frequency scale, handled in 5b
#   n_historic_*  → count bins, handled in 5c
#   neuro_sx_      → three-level symptom recall, handled in 5e
#   non_neuro_sx_ → three-level symptom recall, handled in 5e
#   strang_sign_   → amnesia during strangulation makes Unknown clinically meaningful, handled in 5e
#   head_mech_     → amnesia during head impact makes Unknown clinically meaningful, handled in 5e
#   neuro_sign_    → LOC and other signs may not be recalled due to amnesia, handled in 5e
#
#  FUTURE EXTENSION NOTE 
#
# If the real data intake form is updated to offer Unknown as a valid response
# for dx_ columns (e.g. participant is unsure of their diagnosis), move dx_
# out of yes_no_cols here and into the three_level_cols block in 5e so that
# Unknown is preserved as its own ordered factor level (No < Unknown < Yes)
# rather than collapsed to NA.
# This is consistent with the clinical reasoning applied to neuro_sign_,
# strang_sign_, and head_mech_ — where patient recall uncertainty is
# meaningfully distinct from confirmed absence or missing data.
# Update the Data Dictionary and intake form documentation at the same time.


yes_no_cols <- cohort %>%
  select(
    starts_with("incident_"),
    starts_with("dx_")
  ) %>%
  names()

# case_when evaluates conditions in order and stops at the first match.
# The final TRUE ~ NA is the catch-all: any value that is not "Yes" or "No"
# (blanks, NAs, unexpected strings) becomes NA — true missing data.

cohort <- cohort %>%
  mutate(across(
    all_of(yes_no_cols),
    ~ case_when(
      .x == "Yes" ~ TRUE,
      .x == "No"  ~ FALSE,
      TRUE         ~ NA
    )
  ))

# 5b. Substance use columns → ordered factor
# Substance use is a five-level frequency scale representing use in the past month.
# ordered = TRUE tells R these levels have a meaningful sequence (no use → daily use),
# which affects how comparisons, sorting, and some models treat the variable.
# This is qualitative (never/rarely/sometimes/often/daily) rather than a count.

substance_levels <- c(
  "Never",
  "1-2 days",
  "Several days",
  "More than half the days",
  "Nearly every day"
)

cohort <- cohort %>%
  mutate(across(
    starts_with("substance_"),
    ~ factor(.x, levels = substance_levels, ordered = TRUE)
  ))

# 5c. Historic exposure count bins → ordered factor
# These columns record how many times a participant has experienced strangulation
# or head impact across their lifetime, binned into ranges.
# Unlike 5b (qualitative frequency), these are quantitative count bins —
# "6-10" means more exposures than "1-5", which matters for dose-response analyses.
# ordered = TRUE preserves that ascending exposure gradient.

exposure_levels <- c("0", "1-5", "6-10", "11-15", "15+")

cohort <- cohort %>%
  mutate(
    n_historic_strangulations_ipv = factor(n_historic_strangulations_ipv,
                                           levels = exposure_levels,
                                           ordered = TRUE),
    n_historic_head_impacts_ipv   = factor(n_historic_head_impacts_ipv,
                                           levels = exposure_levels,
                                           ordered = TRUE)
  )

# 5d. Symptom resolution columns → split into logical + numeric
# Each neuro_sx_*_days_resolved column contains mixed data:
#   - an integer (days until symptom resolved)
#   - "Persistent" (symptom ongoing at time of intake)
#   - NA (symptom was absent at intake — no resolution to track)
#
# A single mixed-type column cannot be summarized cleanly, so we split each
# into two purpose-built columns:
#   _persistent      → logical (TRUE = still ongoing, FALSE = resolved, NA = symptom absent)
#   _days_to_resolve → numeric (days until resolution, NA if persistent or absent)
#
# The for loop iterates over every resolution column and applies this split.
# str_remove() strips the suffix to get the symptom base name, then paste0()
# constructs the two new column names dynamically.
#
# Inside mutate():
#   !!flag := uses !! (injection operator) to tell R to treat the value stored
#             in `flag` as the column name rather than the word "flag" itself.
#             := is required instead of = when the left side uses !!.
#   .data[[col]] accesses the current column by name stored in a variable —
#             equivalent to cohort[[col]] but safe inside a pipe.
#   suppressWarnings(as.numeric()) converts day integers to numeric;
#             "Persistent" cannot be converted and becomes NA with a warning —
#             suppressWarnings keeps the console clean since this is expected.

symptom_resolution <- cohort %>%
  select(ends_with("_days_resolved")) %>%
  names()

for (col in symptom_resolution) {
  base  <- str_remove(col, "_days_resolved")
  flag  <- paste0(base, "_persistent")
  days  <- paste0(base, "_days_to_resolve")
  
  cohort <- cohort %>%
    mutate(
      !!flag := .data[[col]] == "Persistent",
      !!days := suppressWarnings(as.numeric(.data[[col]]))
    )
}

# Drop the original mixed-type resolution columns now that they are split
cohort <- cohort %>%
  select(-all_of(symptom_resolution))

# FUTURE EXTENSION NOTE
# If future data collection adds "Unsure" or unquantified responses
# (e.g. "A little while", "A few days") to any _days_resolved column,
# implement the three-column solution below in place of the two-column
# mutate() block above.
#
# IMPORTANT: This solution preserves the clinical distinction between all
# five response categories. These are NOT interchangeable:
#   - "Absent"     → participant never had the symptom
#   - "Resolved"   → participant had the symptom, gave a quantified number of days
#   - "Vague"      → participant HAD the symptom, gave an unquantified time estimate
#                    (e.g. "A little while", "A few days")
#   - "Unsure"     → participant HAD the symptom but could not quantify resolution
#   - "Persistent" → participant currently has the symptom ongoing
#
# "Unsure" and "Vague" must NOT be treated as NA — they are preserved as
# their own categories because these participants DID experience the symptom.
# Grouping them with "Absent" would misrepresent the data clinically.
#
# ALSO UPDATE 5e WHEN IMPLEMENTING: add !ends_with("_resolution_status") to
# the select() filter in 5e, otherwise the new status columns (which also
# start with neuro_sx_) will be caught by the ordered factor recoding and
# produce all NA. The updated 5e filter would read:
#
#   select(
#     !ends_with("_days_to_resolve") &
#     !ends_with("_persistent") &
#     !ends_with("_resolution_status")   # add this line when extension is active
#   )
#
# To implement, replace the mutate() block inside the for loop with:
#
#   persistent_col <- str_replace(col, "_days_resolved", "_persistent")
#   days_col       <- str_replace(col, "_days_resolved", "_days_to_resolve")
#   status_col     <- str_replace(col, "_days_resolved", "_resolution_status")
#
#   cohort <- cohort %>%
#     mutate(
#       !!persistent_col := case_when(
#         .data[[col]] == "Persistent"     ~ TRUE,
#         .data[[col]] == "Unsure"         ~ FALSE,
#         .data[[col]] == "A little while" ~ FALSE,
#         is.na(.data[[col]])              ~ NA,
#         TRUE                             ~ FALSE
#       ),
#       !!days_col := case_when(
#         .data[[col]] == "Persistent"     ~ NA_real_,
#         .data[[col]] == "Unsure"         ~ NA_real_,
#         .data[[col]] == "A little while" ~ NA_real_,
#         is.na(.data[[col]])              ~ NA_real_,
#         TRUE ~ suppressWarnings(as.numeric(.data[[col]]))
#       ),
#       !!status_col := factor(
#         case_when(
#           .data[[col]] == "Persistent"     ~ "Persistent",
#           .data[[col]] == "Unsure"         ~ "Unsure",
#           .data[[col]] == "A little while" ~ "Vague",
#           is.na(.data[[col]])              ~ "Absent",
#           TRUE                             ~ "Resolved"
#         ),
#         levels = c("Absent", "Resolved", "Vague", "Unsure", "Persistent")
#       )
#     )
#
# Also update the Data Dictionary to document any new response options
# added to the data collection tool before modifying this code.


# 5e. Three-level response columns → ordered factor (No < Unknown < Yes)
# Applied to all fields where Unknown is clinically meaningful and distinct
# from both absence (No) and missing data (NA):
#
#   neuro_sx_      → neurological symptoms (headache, dizziness, memory problems, etc.)
#                    A participant with a head injury may have amnesia and genuinely
#                    not recall whether a symptom occurred — Unknown ≠ No
#   non_neuro_sx_ → non-neurological symptoms (neck pain, vision changes, etc.)
#                    Same recall uncertainty applies
#   strang_sign_   → signs during strangulation (petechiae, voice changes, etc.)
#                    Loss of consciousness during strangulation means signs may not be recalled
#   head_mech_     → mechanism of head impact (slapped, thrown, hit object, etc.)
#                    Retrograde amnesia around the impact event makes mechanism uncertain
#   neuro_sign_    → neurological signs at presentation (LOC, seizure, amnesia itself)
#                    A participant cannot reliably self-report LOC or seizure they were
#                    unconscious for — Unknown is a valid and informative response
#
# The !ends_with() filter excludes _days_to_resolve and _persistent columns created
# in 5d — they also start with neuro_sx_ but are numeric/logical, not Yes/No/Unknown,
# and would produce all NA if passed through this factor recoding.
# NOTE: if the future extension in 5d is implemented, also add
# !ends_with("_resolution_status") to this filter — see the note above.
#
# ordered = TRUE means No < Unknown < Yes, so R treats these as having a
# meaningful sequence rather than arbitrary categories.

three_level_cols <- cohort %>%
  select(
    starts_with("neuro_sx_"),
    starts_with("non_neuro_sx_"),
    starts_with("strang_sign_"),
    starts_with("head_mech_"),
    starts_with("neuro_sign_")
  ) %>%
  select(!ends_with("_days_to_resolve") & !ends_with("_persistent")) %>%
  names()

cohort <- cohort %>%
  mutate(across(
    all_of(three_level_cols),
    ~ factor(.x, levels = c("No", "Unknown", "Yes"), ordered = TRUE)
  ))

# Step 6: Assess missingness 

# Per study design, only neuropsychological test scores should be missing
# (~3–4% Missing Completely at Random [MCAR]). Biomarker (NPQz_*) and demographic columns should be complete.
# This section confirms that pattern before any analysis.
# Note: this data is simulated and models data collected from the ASCEND Clinic 
# and previously published IPV work

neuro_test_cols <- c(
  "RAVLT_total", "RAVLT_delayed_recall",
  "TMT_A_Tscore", "TMT_B_Tscore",
  "TOPF_standard_score", "COWA_total_words",
  "SDMT_oral_z", "SDMT_written_z",
  "WAIS_PSI", "WAIS_WMI",
  "TOMM_trial1", "TOMM_trial2",
  "PCL5_total", "PHQ9_total", "NSI_total",
  "BISA_BI_count", "HELPS_total"
)

# Summary table: % missing per neuro test
miss_var_summary(cohort %>% select(all_of(neuro_test_cols)))

# Confirm biomarkers have no missing data
npq_cols <- cohort %>% select(starts_with("NPQz_")) %>% names()
npq_missing <- cohort %>%
  select(all_of(npq_cols)) %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "biomarker", values_to = "n_missing") %>%
  filter(n_missing > 0)

if (nrow(npq_missing) == 0) {
  cat("✓ No missing values in biomarker panel (as expected)\n")
} else {
  cat("Warning: missing biomarker values found:\n")
  print(npq_missing)
}

# Simple visual: which participants and which tests have missing data
# This plot shows missingness as a grid — rows are participants, columns are tests.
vis_miss(cohort %>% select(all_of(neuro_test_cols))) +
  labs(title = "Missing data in neuropsychological test scores",
       subtitle = "Expected: ~3–4% MCAR; biomarker columns excluded")

# medium complexity visualization
# Missingness bar chart by study group
# Each bar = one neuropsychological test
# Shows % missing separately for each group side by side

gg_miss_var(
  cohort %>% select(study_group, all_of(neuro_test_cols)),
  facet = study_group,
  show_pct = TRUE
) +
  labs(
    title    = "Missing neuropsychological test data by study group",
    subtitle = "Expected: ~3-4% MCAR",
    x        = "Neuropsychological test",
    y        = "% missing"
  ) +
  theme_minimal()
# Step 7: Intake Summary Table

# gtsummary::tbl_summary() builds a publication-ready descriptive table.
# by = study_group splits columns by our three groups.
# statistic = controls the format: median (IQR) for continuous, n (%) for categorical.
# missing = "no" suppresses the missing row for cleaner output — we have already
# documented missingness in section 6.

# 7a. Demographics
tbl_demographics <- cohort %>%
  select(
    study_group,
    age_years,
    sex,
    gender,
    sexual_orientation,
    ethnicity,
    education,
    height_cm,
    weight_kg
  ) %>%
  tbl_summary(
    by = study_group,
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      age_years         ~ "Age (years)",
      sex               ~ "Assigned sex at birth",
      gender            ~ "Gender identity",
      sexual_orientation ~ "Sexual orientation",
      ethnicity         ~ "Ethnic identity",
      education         ~ "Highest education",
      height_cm         ~ "Height (cm)",
      weight_kg         ~ "Weight (kg)"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_caption("**Table 1. Participant demographics by study group**")

tbl_demographics

# 7b. Most recent incident characteristics
tbl_incident <- cohort %>%
  select(
    study_group,
    days_since_index_incident,
    incident_head_impact,
    incident_strangulation,
    incident_physical_assault,
    incident_sexual_assault,
    assailant_relationship,
    mechanism_of_strangulation,
    strang_n_times,
    strang_duration
  ) %>%
  tbl_summary(
    by = study_group,
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      days_since_index_incident  ~ "Days since incident (median, IQR)",
      incident_head_impact       ~ "Head impact",
      incident_strangulation     ~ "Strangulation",
      incident_physical_assault  ~ "Physical assault",
      incident_sexual_assault    ~ "Sexual assault",
      assailant_relationship     ~ "Assailant relationship",
      mechanism_of_strangulation ~ "Strangulation type",
      strang_n_times             ~ "Number of strangulation episodes",
      strang_duration            ~ "Duration of strangulation"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_caption("**Table 2. Most recent incident characteristics by study group**")

tbl_incident

# 7c. Substance use
# The substance columns each represent a single substance 
# Participants can answer how often they have used the respective 
# substance within the past month

tbl_substance <- cohort %>%
  select(study_group, starts_with("substance_")) %>%
  tbl_summary(
    by = study_group,
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    label = as.list(
      setNames(
        str_remove(names(select(cohort, starts_with("substance_"))), "substance_"),
        names(select(cohort, starts_with("substance_")))
      )
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_caption("**Table 3. Substance use by study group**")

tbl_substance

# 7d. Cumulative lifetime IPV-BI exposure history
tbl_history <- cohort %>%
  select(
    study_group,
    n_historic_head_impacts_ipv,
    n_historic_strangulations_ipv,
    diagnosed_concussion,
    diagnosed_concussion_n
  ) %>%
  tbl_summary(
    by = study_group,
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      n_historic_head_impacts_ipv   ~ "Lifetime head/face impacts from IPV",
      n_historic_strangulations_ipv ~ "Lifetime strangulations from IPV",
      diagnosed_concussion          ~ "Prior diagnosed concussion",
      diagnosed_concussion_n        ~ "Number of prior concussions"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_caption("**Table 4. Lifetime IPV-BI exposure history by study group**")

tbl_history

# 7e. Non-neurological symptomatology at intake
# non_neuro_sx_ columns are ordered factors (No < Unknown < Yes) from Section 5e.
# Resolution columns (_persistent, _days_to_resolve) are included because this
# is a subacute cohort — symptoms reported at intake may have already partially
# or fully resolved since the index incident, making resolution trajectory
# part of the intake clinical picture, not outcome data.

tbl_non_neuro_sx <- cohort %>%
  select(
    study_group,
    starts_with("non_neuro_sx_")
  ) %>%
  tbl_summary(
    by = study_group,
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_caption("**Table 5. Non-neurological symptoms at intake by study group**")

tbl_non_neuro_sx


# 7f. Neurological symptomatology, signs, and cognitive screening at intake
# neuro_sx_ and neuro_sign_ columns are ordered factors (No < Unknown < Yes) from 5e.
# No mutate() conversion needed — tbl_summary() reads factor levels directly.
# mini_cog_score is kept as continuous (median/IQR) and mini_cog_normal is
# derived here using the standard clinical threshold of ≥4 = Normal.
# Resolution columns included for the same reason as 7e.

tbl_neuro_sx <- cohort %>%
  select(
    study_group,
    starts_with("neuro_sx_"),
    starts_with("neuro_sign_"),
    mini_cog_score,
    mini_cog_normal
  ) %>%
  mutate(
    mini_cog_normal = case_when(
      mini_cog_score >= 4   ~ "Normal (4+)",
      mini_cog_score <  4   ~ "Abnormal (<4)",
      is.na(mini_cog_score) ~ NA_character_
    )
  ) %>%
  tbl_summary(
    by = study_group,
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_caption("**Table 6. Neurological signs, symptoms, and cognitive screening by study group**")

tbl_neuro_sx


# 7g. Other intake information — psychiatric diagnoses and clinical variables
# dx_ columns are logical (TRUE/FALSE/NA) from Section 5a.
# Converted here to Yes/No/Unknown for clean display in the summary table.
# This conversion is display-only — it occurs inside the pipe and does not
# modify the cohort dataset itself.
# Note: Unknown here reflects true missing data, not patient uncertainty —
# Unknown was not an offered response option for diagnosis fields on the
# intake form. See the future extension note in 5a if this changes.

tbl_other_intake <- cohort %>%
  select(
    study_group,
    starts_with("dx_"),
    suicide_ideation,
    etoh_drugs_to_cope,
    sleep_disturbance
  ) %>%
  mutate(
    across(
      starts_with("dx_"),
      ~ case_when(
        .x == TRUE  ~ "Yes",
        .x == FALSE ~ "No",
        is.na(.x)   ~ "Unknown"
      )
    )
  ) %>%
  tbl_summary(
    by = study_group,
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    label = list(
      dx_ptsd            ~ "PTSD diagnosis",
      dx_adhd            ~ "ADHD diagnosis",
      dx_anxiety         ~ "Anxiety disorder diagnosis",
      dx_depression      ~ "Depression diagnosis",
      dx_bipolar         ~ "Bipolar disorder diagnosis",
      dx_other           ~ "Other psychiatric diagnosis",
      suicide_ideation   ~ "Suicide ideation",
      etoh_drugs_to_cope ~ "Uses alcohol/drugs to cope",
      sleep_disturbance ~ "Sleep disturbances reported"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_caption("**Table 7. Psychiatric diagnoses and other clinical intake information by study group**")

tbl_other_intake

# Step 8: Substance Use Prevalence

# 8a. Frequency of substance use in the past month (one loop per substance)
# picks up any new substances automatically via starts_with("substance_")

substance_cols <- cohort %>% select(starts_with("substance_")) %>% names()

substance_summary <- vector("list", length(substance_cols))

for (i in seq_along(substance_cols)) {
  col <- substance_cols[i]
  
  substance_summary[[i]] <- cohort %>%
    group_by(study_group) %>%
    summarise(
      substance = str_remove(col, "substance_"),
      n_never   = sum(.data[[col]] == "Never", na.rm = TRUE),
      n_low     = sum(.data[[col]] == "1-2 days", na.rm = TRUE),
      n_several = sum(.data[[col]] == "Several days", na.rm = TRUE),
      n_half    = sum(.data[[col]] == "More than half the days", na.rm = TRUE),
      n_daily   = sum(.data[[col]] == "Nearly every day", na.rm = TRUE),
      n_total   = sum(!is.na(.data[[col]])),
      .groups   = "drop"
    )
}

substance_prev <- bind_rows(substance_summary)
print(substance_prev)

#  FUTURE EXTENSION NOTE 
# If future data collection adds a variable for substance use during the
# most recent IPV incident, this should be handled separately from the
# frequency scale above because it is a distinct clinical question —
# not how often a participant uses substances generally, but whether
# substances were present during a specific event.
#
# The new variable would likely be Yes/No/Unknown per substance, so it
# should be coded as a logical column (TRUE/FALSE/NA) following the same
# pattern as Section 5a, NOT as an ordered factor like Section 5b.
#
# To implement, add a new column to the dataset (e.g. substance_during_incident_alcohol)
# and add a second loop below:
#
# 8b. Substance use during most recent IPV incident
#
# incident_substance_cols <- cohort %>%
#   select(starts_with("substance_during_incident_")) %>%
#   names()
#
# incident_substance_summary <- vector("list", length(incident_substance_cols))
#
# for (i in seq_along(incident_substance_cols)) {
#   col <- incident_substance_cols[i]
#
#   incident_substance_summary[[i]] <- cohort %>%
#     group_by(study_group) %>%
#     summarise(
#       substance  = str_remove(col, "substance_during_incident_"),
#       n_yes      = sum(.data[[col]] == TRUE, na.rm = TRUE),
#       n_no       = sum(.data[[col]] == FALSE, na.rm = TRUE),
#       n_unknown  = sum(is.na(.data[[col]])),
#       n_total    = n(),
#       pct_yes    = round(n_yes / n_total * 100, 1),
#       .groups    = "drop"
#     )
# }
#
# incident_substance_prev <- bind_rows(incident_substance_summary)
# print(incident_substance_prev)
#
# Also add substance_during_incident_* columns to the Data Dictionary
# and to the tbl_substance summary table in Section 7c as a separate
# row group so it is visually distinct from the frequency scale above.

# Step 9: Save Clean Dataset

# Save the cleaned dataframe as an .rds file — R's native format.
# .rds preserves factor levels, logical types, and all recoding done above,
# so future scripts can load it without repeating the cleaning steps.

saveRDS(cohort, file = "ipv_cohort_clean.rds")
cat("Cleaned dataset saved to ipv_cohort_clean.rds\n")
cat("Rows:", nrow(cohort), "| Columns:", ncol(cohort), "\n")

# Step 10: Export summary data to excel tables, easier visualization 

getwd()
setwd("C:/NRSC_510")

# Create a new workbook
wb <- createWorkbook()

# Convert each gtsummary table to a dataframe and write to its own sheet
tables_to_export <- list(
  "Demographics"       = tbl_demographics,
  "Incident"           = tbl_incident,
  "Substance Use"      = tbl_substance,
  "IPV History"        = tbl_history,
  "Non-Neuro Symptoms" = tbl_non_neuro_sx,
  "Neuro Symptoms"     = tbl_neuro_sx,
  "Other Intake"       = tbl_other_intake
)

for (sheet_name in names(tables_to_export)) {
  addWorksheet(wb, sheet_name)
  writeData(
    wb,
    sheet  = sheet_name,
    x      = as_tibble(tables_to_export[[sheet_name]]),
    rowNames = FALSE
  )
}

saveWorkbook(wb, "intake_summary_tables.xlsx", overwrite = TRUE)

# In Week 3+, load the cleaned data with:
# cohort <- readRDS("ipv_cohort_clean.rds")


