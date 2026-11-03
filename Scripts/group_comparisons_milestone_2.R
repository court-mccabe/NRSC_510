# Week 3-4: Group Comparisons
# Compares demographics, IPV-BI exposure, injury mechanisms, symptom reports,
# substance use, and psychiatric diagnoses across the three study groups:
#   - Acute/Subacute IPV-BI
#   - Acute/Subacute non-IPV-BI
#   - Chronic IPV-BI

# Multiple comparison correction: Benjamini-Hochberg False Discovery Rate (BH-FDR)
# BH-FDR controls the proportion of false positives among significant findings
# rather than the probability of any single false positive (Tukey/Bonferroni).
# This is appropriate for exploratory analyses with many simultaneous comparisons
# where the goal is hypothesis generation rather than confirmatory inference.

# FUTURE EXTENSION NOTE
# If sample size increases , revisit the correction method used at each stage:
#   - Section 11 (demographic ANOVA pairwise comparisons): switch to Tukey's
#     HSD, which is the standard expectation in clinical neuropsychology
#     publications for Table 1 demographic comparisons. With only three
#     pairwise comparisons the conservativeness penalty is minimal at larger n.
#   - Sections 13-15 (ordinal/linear/logistic models across symptom/exposure
#     columns): BH-FDR remains appropriate even at larger n given the number
#     of simultaneous tests and the exploratory aim. State clearly in the
#     methods section that results are hypothesis-generating.
# Also reassess whether two covariates (age and sex) remain sufficient or
# whether additional covariates (e.g. education, time since incident) are
# warranted at larger n — return to the model specification in Section 13
# when doing so.

# Limitation: 
# Persistent symptom columns (particularly non-neuro) have low event counts
# per study group in this simulated cohort — some symptoms have as few as
# 0-6 persistent cases per group. Logistic regression requires a reasonable
# number of events per group to produce stable, well-powered estimates
# (general guidance: at least ~10 events per predictor level). Non-significant
# results for persistent symptom comparisons likely reflect low statistical
# power rather than a confirmed absence of a true group difference. Report
# as a limitation. At larger n, Firth's penalized logistic regression
# (logistf package) is a more robust alternative if rare events persist.

library(tidyverse)
library(emmeans)
library(MASS)      # polr() for ordinal logistic regression
library(broom)     # tidy() for clean model output
library(ggplot2)
library(ggpubr)
library(gt)

getwd()
setwd("C:/NRSC_510")
cohort <- readRDS("ipv_cohort_clean.rds")

# Explicitly resolve MASS::select conflict with dplyr::select
# MASS masks dplyr's select function when loaded — this ensures
# select() always refers to dplyr throughout this script
select <- dplyr::select

# Step 11: Continuous Variables ANOVA
# One-way ANOVA tests whether any group difference exists for continuous
# demographic variables. emmeans then extracts pairwise comparisons between
# all three study groups with BH-FDR correction applied across all pairs.

# Note: demographics are reported unadjusted (no covariates) because they ARE
# the covariates in later steps — adjusting for age when comparing age across
# groups is circular.

# 11a. Define ANOVA columns
# These are the continuous demographic variables being compared across groups.

age_height_weight_anova_cols <- c("age_years", "height_cm", "weight_kg")

# 11b. Run one-way ANOVA and extract overall p-value
# aov() fits a one-way ANOVA model with study_group as the predictor.
# summary(model)[[1]][["Pr(>F)"]][1] extracts the overall group effect p-value.

age_height_weight_anova_results <- lapply(age_height_weight_anova_cols, function(var) {
  formula <- as.formula(paste(var, "~ study_group"))
  model   <- aov(formula, data = cohort)
  
  anova_p <- summary(model)[[1]][["Pr(>F)"]][1]
  
# 11c. Run emmeans pairwise comparisons with BH-FDR
  # emmeans extracts estimated marginal means for each study group, then
  # computes all pairwise contrasts between them.
  # estimate: mean difference between the two groups in the contrast
  # SE: standard error of that difference
  # df: degrees of freedom used in the t-test for that contrast
  # t.ratio: the test statistic (estimate divided by SE)
  # p.value: raw p-value before correction
  # p.adj.BH: BH-FDR corrected p-value
  
  emm   <- emmeans(model, pairwise ~ study_group)
  pairs <- as.data.frame(emm$contrasts) %>%
    mutate(
      variable = var,
      p.adj.BH = p.adjust(p.value, method = "BH")
    ) %>%
    select(variable, contrast, estimate, SE, df, t.ratio, p.value, p.adj.BH)
  
  list(anova_p = anova_p, pairs = pairs)
})

names(age_height_weight_anova_results) <- age_height_weight_anova_cols

# 11d. Combine and print pairwise results

cat("Overall ANOVA p-values:\n")
for (var in age_height_weight_anova_cols) {
  cat(sprintf("  %-20s p = %.4f\n", var, age_height_weight_anova_results[[var]]$anova_p))
}

age_height_weight_pairwise <- bind_rows(
  lapply(age_height_weight_anova_results, `[[`, "pairs")
)
print(age_height_weight_pairwise)

# 11e. Visualize age comparison across study groups
# A simple boxplot showing the distribution of age in each study group,
# with significance brackets from the emmeans pairwise comparisons in 11c.
# This gives a visual preview of what a significant (or non-significant)
# group difference looks like before moving into the more complex ordinal
# and logistic models in later sections.

age_group_comparisons <- list(
  c("Acute IPV-BI", "Acute Non-IPV-BI"),
  c("Acute IPV-BI", "Chronic IPV-BI"),
  c("Acute Non-IPV-BI", "Chronic IPV-BI")
)

age_boxplot <- ggplot(cohort, aes(x = study_group, y = age_years, fill = study_group)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8, width = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
  stat_compare_means(
    comparisons     = age_group_comparisons,
    method          = "t.test",
    p.adjust.method = "BH",
    label           = "p.signif"
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Age by Study Group",
    x     = "Study Group",
    y     = "Age (years)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 15, hjust = 1),
    plot.title      = element_text(face = "bold", hjust = 0.5)
  )

print(age_boxplot)


# Step 12: Categorical demographic comparisons 
# Chi-square requires expected cell counts of at least 5 to be reliable.
# Fisher's exact test is used instead when any expected cell count falls
# below 5. BH-FDR correction applied across all categorical variables.
# Note: Fisher's exact is expected for most categorical variables in this
# simulated cohort given the small sample size and number of levels per
# variable — this is not an error.

# 12a. Define categorical demographic columns

sex_edu_ethnicity_chi_cols <- c(
  "sex",
  "gender",
  "ethnicity",
  "education",
  "sexual_orientation"
)

# 12b. Run chi-square or Fisher's exact per column
# table() builds a contingency table of each variable against study_group.
# chisq.test(tbl)$expected checks expected cell counts before choosing
# the test. simulate.p.value = TRUE uses Monte Carlo simulation for
# Fisher's exact, which is more reliable for larger tables.

sex_edu_ethnicity_chi_results <- lapply(sex_edu_ethnicity_chi_cols, function(var) {
  tbl      <- table(cohort[[var]], cohort$study_group)
  expected <- chisq.test(tbl)$expected
  
  test <- if (any(expected < 5)) {
    fisher.test(tbl, simulate.p.value = TRUE)
  } else {
    chisq.test(tbl)
  }
  
  # 12c. Record which test was used alongside the raw p-value
  data.frame(
    variable = var,
    test     = if (any(expected < 5)) "Fisher's exact" else "Chi-square",
    p.value  = test$p.value
  )
})

# 12d. Apply BH-FDR correction and print results

sex_edu_ethnicity_chi_results <- bind_rows(sex_edu_ethnicity_chi_results) %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH"))

print(sex_edu_ethnicity_chi_results)

sex_edu_ethnicity_chi_results %>%
  gt() %>%
  tab_header(title = "Categorical Demographic Comparisons") %>%
  cols_label(
    variable = "Variable",
    test     = "Test Used",
    p.value  = "p-value",
    p.adj.BH = "Adjusted p (BH-FDR)"
  ) %>%
  fmt_number(columns = c(p.value, p.adj.BH), decimals = 4) %>%
  tab_style(
    style = cell_fill(color = "#fff3cd"),
    locations = cells_body(rows = p.adj.BH < 0.05)
  )


# Step 13: IPV-BI exposure mechanism comparisons 
# Tests whether reported exposure mechanisms differ across study groups,
# adjusted for age and sex. Outcomes are ordered factors (No < Unknown < Yes)
# from Section 5e, so ordinal logistic regression (MASS::polr) is used.
# emmeans extracts covariate-adjusted pairwise contrasts. BH-FDR applied
# across all exposure columns simultaneously.

# 13a. Define exposure mechanism columns
# Captures head impact mechanisms, strangulation signs, and incident type —
# all recoded as ordered factors in Section 5e.

ipv_exposure_mechanism_cols <- cohort %>%
  select(
    starts_with("head_mech_"),
    starts_with("strang_sign_"),
    starts_with("incident_")
  ) %>%
  names()

# 13b. Define reusable binary logistic emmeans function
binary_emmeans_by_group <- function(outcome, data, covariates = c("age_years", "sex")) {
  formula <- as.formula(
    paste(outcome, "~", paste(c("study_group", covariates), collapse = " + "))
  )
  tryCatch({
    model     <- glm(formula, data = data, family = binomial)
    emm       <- emmeans(model, pairwise ~ study_group, type = "response")
    contrasts <- as.data.frame(emm$contrasts) %>%
      mutate(outcome = outcome)
    contrasts
  }, error = function(e) {
    data.frame(outcome = outcome, error = conditionMessage(e))
  })
}

# 13c. Run across all exposure columns
ipv_exposure_mechanism_contrasts <- bind_rows(
  lapply(ipv_exposure_mechanism_cols, binary_emmeans_by_group, data = cohort)
) %>%
  filter(!is.na(p.value)) %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH"))

print(ipv_exposure_mechanism_contrasts)

# 13d. Produce table of statistical comparisons

ipv_exposure_mechanism_contrasts %>%
  mutate(
    contrast = gsub("study_group", "", contrast),  # clean up label if needed
    sig = case_when(
      p.adj.BH < 0.05  ~ "*",
      TRUE             ~ ""
    )
  ) %>%
  select(outcome, contrast, odds.ratio, SE, z.ratio, p.value, p.adj.BH, sig) %>%
  gt(groupname_col = "outcome") %>%
  tab_header(
    title    = "IPV-BI Exposure Mechanism Comparisons",
    subtitle = "Binary logistic regression, adjusted for age and sex | BH-FDR corrected"
  ) %>%
  cols_label(
    contrast   = "Comparison",
    odds.ratio = "Odds Ratio",
    SE         = "SE",
    z.ratio    = "z",
    p.value    = "p-value",
    p.adj.BH   = "Adjusted p (BH)",
    sig        = ""
  ) %>%
  fmt_number(columns = c(odds.ratio, SE, z.ratio), decimals = 3) %>%
  fmt_number(columns = c(p.value, p.adj.BH), decimals = 3) %>%
  tab_style(
    style     = cell_fill(color = "#fff3cd"),
    locations = cells_body(rows = p.adj.BH < 0.05)
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_footnote(
    footnote  = "* p<0.05 (BH-FDR adjusted)",
    locations = cells_column_labels(columns = sig)
  )

# table shows odds ratio, SE, z, ratio, and p value
# no significant differences found between groups

# Step 14: Symptom Comparisons
# Tests whether symptom reporting differs across study groups, adjusted for
# age and sex. Neurological (neuro_sx_) and non-neurological (non_neuro_sx_)
# symptoms are kept separate throughout — a group difference in neuro
# symptoms that does not appear in non-neuro symptoms is clinically
# meaningful and would be obscured if the two were combined.
#
# Three model types depending on column type:
#   - Symptom presence (No < Unknown < Yes) → ordinal logistic regression
#   - Days to resolve (numeric) → linear regression
#   - Persistent symptoms (TRUE/FALSE) → logistic regression
# BH-FDR correction applied within each symptom type and column type
# separately — combining them would dilute the correction across clinically
# distinct outcome families.

# 14a. Define symptom column types — neuro and non-neuro kept separate
# Chained select() calls: starts_with() narrows to the right prefix first,
# then ends_with() filters by suffix within that subset.

neuro_presence_cols <- cohort %>%
  select(starts_with("neuro_sx_")) %>%
  select(!ends_with("_days_to_resolve") & !ends_with("_persistent")) %>%
  names()

non_neuro_presence_cols <- cohort %>%
  select(starts_with("non_neuro_sx_")) %>%
  select(!ends_with("_days_to_resolve") & !ends_with("_persistent")) %>%
  names()

neuro_days_cols <- cohort %>%
  select(starts_with("neuro_sx_")) %>%
  select(ends_with("_days_to_resolve")) %>%
  names()

non_neuro_days_cols <- cohort %>%
  select(starts_with("non_neuro_sx_")) %>%
  select(ends_with("_days_to_resolve")) %>%
  names()

neuro_persistent_cols <- cohort %>%
  select(starts_with("neuro_sx_")) %>%
  select(ends_with("_persistent")) %>%
  names()

non_neuro_persistent_cols <- cohort %>%
  select(starts_with("non_neuro_sx_")) %>%
  select(ends_with("_persistent")) %>%
  names()

# 14b. Define reusable linear regression emmeans function
# Used for continuous outcomes (days to resolve). lm() fits ordinary least
# squares. emmeans() extracts the adjusted mean days to resolve per group
# holding age and sex constant, then computes pairwise contrasts.

linear_emmeans_by_group <- function(outcome, data, covariates = c("age_years", "sex")) {
  formula <- as.formula(
    paste(outcome, "~", paste(c("study_group", covariates), collapse = " + "))
  )
  tryCatch({
    model     <- lm(formula, data = data)
    emm       <- emmeans(model, pairwise ~ study_group)
    contrasts <- as.data.frame(emm$contrasts) %>%
      mutate(outcome = outcome)
    contrasts
  }, error = function(e) {
    data.frame(outcome = outcome, error = conditionMessage(e))
  })
}

# 14c. Define reusable logistic regression emmeans function
# Used for binary outcomes (persistent symptoms TRUE/FALSE).
# glm(family = binomial) fits logistic regression. type = "response"
# returns probability-scale contrasts rather than log-odds.

logistic_emmeans_by_group <- function(outcome, data, covariates = c("age_years", "sex")) {
  formula <- as.formula(
    paste(outcome, "~", paste(c("study_group", covariates), collapse = " + "))
  )
  tryCatch({
    model     <- glm(formula, data = data, family = binomial)
    emm       <- emmeans(model, pairwise ~ study_group, type = "response")
    contrasts <- as.data.frame(emm$contrasts) %>%
      mutate(outcome = outcome)
    contrasts
  }, error = function(e) {
    data.frame(outcome = outcome, error = conditionMessage(e))
  })
}

# 14d. Run models across all symptom column types
# Neuro and non-neuro run separately into distinct results tables so group
# differences can be interpreted independently. filter removes rows where
# the model failed to converge.

#cohort %>%  getting error with code below, checking what neuro cols are ordered as
 # select(all_of(neuro_presence_cols)) %>%
  # sapply(class)

# Fix: recode neuro_sx columns as ordered factors
cohort <- cohort %>%
  mutate(across(
    all_of(neuro_presence_cols),
    ~ factor(., levels = c("No", "Unknown", "Yes"), ordered = TRUE)
  ))

# Verify
cohort %>%
  select(all_of(neuro_presence_cols)) %>%
  sapply(class)

bind_rows(
  lapply(neuro_presence_cols, ordinal_emmeans_by_group, data = cohort)
) %>% names()

#neuro_presence_contrasts <- bind_rows(
 #pply(neuro_presence_cols, ordinal_emmeans_by_group, data = cohort)
#%>%
 # filter(!is.na(p.value)) # producing error of p.value name

bind_rows(  # check to see what's failing
  lapply(non_neuro_presence_cols, ordinal_emmeans_by_group, data = cohort)
) %>% pull(error) %>% unique()

cohort %>%
  select(all_of(non_neuro_presence_cols)) %>%
  sapply(class)

cohort <- cohort %>%
  mutate(across(
    all_of(non_neuro_presence_cols),
    ~ factor(., levels = c("No", "Unknown", "Yes"), ordered = TRUE)
  ))

neuro_presence_contrasts <- bind_rows(
  lapply(neuro_presence_cols, ordinal_emmeans_by_group, data = cohort)
) %>%
  filter(!is.na(p.value))

non_neuro_presence_contrasts <- bind_rows(
  lapply(non_neuro_presence_cols, ordinal_emmeans_by_group, data = cohort)
) %>%
  filter(!is.na(p.value))

neuro_days_contrasts <- bind_rows(
  lapply(neuro_days_cols, linear_emmeans_by_group, data = cohort)
) %>%
  filter(!is.na(p.value))

non_neuro_days_contrasts <- bind_rows(
  lapply(non_neuro_days_cols, linear_emmeans_by_group, data = cohort)
) %>%
  filter(!is.na(p.value))

neuro_persistent_contrasts <- bind_rows(
  lapply(neuro_persistent_cols, logistic_emmeans_by_group, data = cohort)
) %>%
  filter(!is.na(p.value))

non_neuro_persistent_contrasts <- bind_rows(
  lapply(non_neuro_persistent_cols, logistic_emmeans_by_group, data = cohort)
) %>%
  filter(!is.na(p.value))

# 14e. Apply BH-FDR correction and print results

neuro_presence_contrasts <- neuro_presence_contrasts %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH"))

non_neuro_presence_contrasts <- non_neuro_presence_contrasts %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH"))

neuro_days_contrasts <- neuro_days_contrasts %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH"))

non_neuro_days_contrasts <- non_neuro_days_contrasts %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH"))

neuro_persistent_contrasts <- neuro_persistent_contrasts %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH"))

non_neuro_persistent_contrasts <- non_neuro_persistent_contrasts %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH"))

print(neuro_presence_contrasts)
print(non_neuro_presence_contrasts)
print(neuro_days_contrasts)
print(non_neuro_days_contrasts)
print(neuro_persistent_contrasts)
print(non_neuro_persistent_contrasts)

# 14f. Diagnostic — summarize significant findings in plain language

all_findings <- bind_rows(
  neuro_presence_contrasts       %>% mutate(category = "Neuro symptom presence"),
  non_neuro_presence_contrasts   %>% mutate(category = "Non-neuro symptom presence"),
  neuro_days_contrasts           %>% mutate(category = "Neuro days to resolve"),
  non_neuro_days_contrasts       %>% mutate(category = "Non-neuro days to resolve"),
  neuro_persistent_contrasts     %>% mutate(category = "Neuro persistent"),
  non_neuro_persistent_contrasts %>% mutate(category = "Non-neuro persistent")
) %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH"))

all_findings %>%
  gt(groupname_col = "category") %>%
  tab_header(
    title    = "Symptom Comparisons",
    subtitle = "BH-FDR adjusted | adjusted for age and sex"
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style     = cell_fill(color = "#fff3cd"),
    locations = cells_body(rows = p.adj.BH < 0.05)
  )
# Step 15: Compare Substance Abuste and Psychiatric Diagnoses
# Tests whether substance use frequency and psychiatric diagnosis rates
# differ across study groups, adjusted for age and sex.
# Substance use (ordered factor) → ordinal logistic regression, reusing
# ordinal_emmeans_by_group from Section 13.
# dx_ columns (logical TRUE/FALSE) → logistic regression, reusing
# logistic_emmeans_by_group from Section 14.
# BH-FDR correction applied within each column type separately.

# 15a. Define substance use and diagnosis columns

substance_use_cols <- cohort %>%
  select(starts_with("substance_")) %>%
  names()

dx_diagnosis_cols <- cohort %>%
  select(starts_with("dx_")) %>%
  names()

# 15b. Run models across both column types
# filter removes rows where the model failed to converge.

substance_use_cols # problematic naming convention with spaces and / used

cohort <- cohort %>%
  rename(
    substance_pain_killers  = "substance_Pain killers",
    substance_MDMA_ecstasy  = "substance_MDMA/ecstasy",
    substance_meth          = "substance_Meth",
    substance_alcohol       = "substance_Alcohol",
    substance_stimulants    = "substance_Stimulants",
    substance_sedatives     = "substance_Sedatives",
    substance_marijuana     = "substance_Marijuana",
    substance_cocaine       = "substance_Cocaine",
    substance_hallucinogens = "substance_Hallucinogens",
    substance_heroin        = "substance_Heroin",
    substance_inhalants     = "substance_Inhalants"
  )

# Redefine the column vector with the new names
substance_use_cols <- c(
  "substance_alcohol", "substance_pain_killers", "substance_stimulants",
  "substance_sedatives", "substance_marijuana", "substance_cocaine",
  "substance_MDMA_ecstasy", "substance_hallucinogens", "substance_heroin",
  "substance_inhalants", "substance_meth"
)

substance_use_contrasts <- bind_rows(
  lapply(substance_use_cols, ordinal_emmeans_by_group, data = cohort)
) %>%
  filter(!is.na(p.value))

dx_diagnosis_contrasts <- bind_rows(
  lapply(dx_diagnosis_cols, logistic_emmeans_by_group, data = cohort)
) %>%
  filter(!is.na(p.value))

# 15c. Apply BH-FDR correction and print results

substance_use_contrasts <- substance_use_contrasts %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH"))

dx_diagnosis_contrasts <- dx_diagnosis_contrasts %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH"))

print(substance_use_contrasts)
print(dx_diagnosis_contrasts)

# 15d. Summarize significant findings

bind_rows(
  substance_use_contrasts  %>% mutate(category = "Substance use"),
  dx_diagnosis_contrasts   %>% mutate(category = "Psychiatric diagnosis")
) %>%
  mutate(
    p.adj.BH = p.adjust(p.value, method = "BH"),
    sig = case_when(
      p.adj.BH < 0.05  ~ "*",
      TRUE             ~ ""
    )
  ) %>%
  gt(groupname_col = "category") %>%
  tab_header(
    title    = "Substance Use & Psychiatric Diagnosis Comparisons",
    subtitle = "BH-FDR adjusted | adjusted for age and sex"
  ) %>%
  cols_label(sig = "") %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style     = cell_fill(color = "#fff3cd"),
    locations = cells_body(rows = p.adj.BH < 0.05)
  ) %>%
  tab_footnote(
    footnote  = "* p<0.05  ** p<0.01  *** p<0.001 (BH-FDR adjusted)",
    locations = cells_column_labels(columns = sig)
  )
