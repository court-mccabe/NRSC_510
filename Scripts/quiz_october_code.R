# Update file path
file_path <- "C:/NRSC_510/IPV_cohort_simulated_v5.xlsx"

#combine data
acute_raw   <- read_excel(file_path, sheet = "Acute_Subacute")
chronic_raw <- read_excel(file_path, sheet = "Chronic")
cohort_raw <- bind_rows(acute_raw, chronic_raw)

#neuropsych data
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

#identify characters
char_missing_cols <- neuro_test_cols[
  sapply(cohort_raw[neuro_test_cols], is.character)
]

# Apply na.if
cohort <- cohort_raw %>%
  mutate(across(all_of(char_missing_cols), ~ na_if(.x, "N/A")))

#Create clean group labels
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

#identify binary variables
yes_no_cols <- cohort %>%
  select(
    starts_with("incident_"),
    starts_with("dx_")
  ) %>%
  names()

cohort <- cohort %>%
  mutate(across(
    all_of(yes_no_cols),
    ~ case_when(
      .x == "Yes" ~ TRUE,
      .x == "No"  ~ FALSE,
      TRUE         ~ NA
    )
  ))