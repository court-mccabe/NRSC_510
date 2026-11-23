# Week 7-8: Neuropsychological Impairment and Biomarker-Cognition Correlations
# Goal: Assess group differences in neuropsychological performance across the
# three study groups, then relate the full NULISA biomarker panel to
# neuropsychological outcomes.
#
# Analysis order:
#   Steps 44 and 45   group comparisons of the 13 neuropsychological scores
#   Steps 46 to 52    all proteins against all 13 outcomes (hypothesis generating):
#                     correlations in Steps 47 and 48, a top ten table in Step 49,
#                     heatmaps in Step 50, linear regression in Step 51, and a
#                     regression heatmap with volcano plots in Step 52
#   Step 53           targeted test of the top milestone 3 proteins against
#                     selected outcomes (hypothesis driven)
#
# Groups:
#   - Acute IPV-BI       (<=30 days post-incident, IPV-caused brain injury)
#   - Acute Non-IPV-BI   (<=30 days post-incident, no IPV-BI)
#   - Chronic IPV-BI     (>6 months post-incident, IPV-BI)
#
# Multiple comparison correction: Benjamini-Hochberg FDR throughout.
# Correlation method: Shapiro-Wilk fallback: Pearson if both variables
# pass normality (p > 0.05), Spearman if either fails. This is appropriate
# at n = 66 where parametric assumptions cannot be taken for granted.
#
# Score standardization note (limitation):
# Several neuropsychological tests are scored here as within-sample z-scores
# because population-based normative T-score tables require demographic lookup
# (age, sex, education) that is feasible with real data but not implemented
# here with simulated data. Tests that already provide T-scores or standard
# scores in the dataset (TMT-A, TMT-B, TOPF, WAIS PSI, WAIS WMI) are used
# as-is. For a real dataset, raw scores would arrive pre-converted to T-scores
# at the point of data entry; see pseudocode in Step 42.
#
# Limitations:
#   1. Sample size. With n = 66 and 2873 tests in the all-protein analyses
#      (Steps 47 and 51), power after BH correction is low. Only large
#      associations can be detected, and the strongest hits tend to
#      overstate the true effect. These steps are hypothesis generating.
#   2. Pooled groups. Steps 47 and 51 combine all three study groups and do
#      not adjust for study group, so an association can reflect a group
#      difference rather than a relationship within groups. Step 53 fits
#      models with and without study group to show the difference.
#   3. Model assumptions. Linear regression assumes straight line
#      relationships, roughly normal residuals, and no single point
#      dominating the fit. Residuals were not checked model by model.
#   4. Related outcomes. Several scores are related (for example RAVLT total
#      and delayed recall, or SDMT oral and written), so one true finding can
#      appear as several hits and the 13 outcomes are not independent tests.
#   5. Scale of the scores. PCL-5, PHQ-9 and NSI are summed questionnaire
#      items treated as continuous. Raw scores are standardized within this
#      sample rather than against population norms (see the note above).
#   6. Unmeasured confounding. Education, kidney function (relevant to CST3),
#      and other factors are not in the models.
#   7. Design. Data come from a single time point, so associations cannot show
#      whether a protein affects cognition or the reverse. Each model drops
#      participants missing that score, so n varies from about 60 to 66.
#   8. Simulated data. Results reflect how the data were built and are not
#      biological findings.
#
# References:
#   Ferreira & Castro (2022). Linear and logistic regression models: when to
#     use and how to interpret them. J Bras Pneumol. PMC9747134.
#   Wikipedia: Generalized linear mixed model
#   Wikipedia: Likelihood-ratio test
#   MathWorks: Wilkinson notation


# Step 41: Load packages and data

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(openxlsx)
library(ggrepel)

setwd("C:/NRSC_510")
cohort <- readRDS("ipv_cohort_clean.rds")

output_dir <- "C:/NRSC_510/Outputs"

select <- dplyr::select

qc_meta_cols <- c("nc_mean_npq", "nc_sd_npq", "nc_pass",
                  "ic_reads", "ic_pct_of_median", "ic_pass",
                  "ipc_cv_pct", "ipc_pass")

# Step 42: Define neuropsychological column groupings, score types, and flags
# Score type key:
#   "tscore"   - already a T-score (mean 50, SD 10) in the dataset; used as-is
#   "ss"       - already a standard score (mean 100, SD 15); used as-is
#   "zscore"   - already a z-score in the dataset; used as-is
#   "raw_z"    - raw score; will be within-sample z-scored in Step 43
#
# Pseudocode for real data import:
#   When real data are imported, raw scores for RAVLT, COWA, PCL-5, PHQ-9,
#   NSI, BISA, and HELPS would arrive pre-converted to normative T-scores
#   using published norms (e.g., Schmidt 1996 for RAVLT; Ruff et al. 1996
#   for COWA; Weathers et al. 2013 for PCL-5). The score_type column for
#   those tests would then be changed to "tscore" and Step 43 would skip
#   them. No conversion code is needed here; the T-score column would be
#   the data column at import.

score_meta <- tibble::tribble(
  ~col,                    ~label,                          ~score_type, ~domain,
  "RAVLT_total",           "RAVLT Total",                   "raw_z",     "cognitive",
  "RAVLT_delayed_recall",  "RAVLT Delayed Recall",          "raw_z",     "cognitive",
  "TMT_A_Tscore",          "TMT-A",                         "tscore",    "cognitive",
  "TMT_B_Tscore",          "TMT-B",                         "tscore",    "cognitive",
  "TOPF_standard_score",   "TOPF",                          "ss",        "cognitive",
  "COWA_total_words",      "COWA",                          "raw_z",     "cognitive",
  "SDMT_oral_z",           "SDMT Oral",                     "zscore",    "cognitive",
  "SDMT_written_z",        "SDMT Written",                  "zscore",    "cognitive",
  "WAIS_PSI",              "WAIS PSI",                      "ss",        "cognitive",
  "WAIS_WMI",              "WAIS WMI",                      "ss",        "cognitive",
  "PCL5_total",            "PCL-5",                         "raw_z",     "psychiatric",
  "PHQ9_total",            "PHQ-9",                         "raw_z",     "psychiatric",
  "NSI_total",             "NSI",                           "raw_z",     "psychiatric"
)

cognitive_cols    <- score_meta %>% filter(domain == "cognitive")    %>% pull(col)
psychiatric_cols  <- score_meta %>% filter(domain == "psychiatric")  %>% pull(col)
all_outcome_cols  <- score_meta %>% pull(col)
exposure_cols     <- c("BISA_BI_count", "HELPS_total")

# TOMM validity flag: participants scoring below 45 on Trial 2 are flagged.
# TOMM is not included in group comparisons or correlations; it is a
# performance validity test. A flag is recorded here for transparency and
# could be used as an exclusion covariate in sensitivity analyses.
cohort <- cohort %>%
  mutate(tomm_flag = !is.na(TOMM_trial2) & TOMM_trial2 < 45)

n_flagged <- sum(cohort$tomm_flag, na.rm = TRUE)
if (n_flagged > 0) {
  cat("TOMM validity flag: ", n_flagged,
      "participant(s) scored below threshold (Trial 2 < 45)\n")
}

# Step 43: Standardize scores to a common scale for comparison and plotting
# T-scores (mean 50, SD 10) and standard scores (mean 100, SD 15) are
# rescaled to z-scores so all measures share a mean-0, SD-1 scale on plots.
# Raw columns (raw_z) are within-sample z-scored. Z-score columns are used
# as-is. This produces a standardized dataframe used only for visualization
# and group comparison; the original columns in cohort are not overwritten.
#
# Limitation: within-sample z-scoring anchors performance to the cohort mean
# rather than a normative population mean. This means a score of 0 reflects
# average performance within this IPV-exposed sample, not the general
# population. In a real dataset, normative T-scores would replace raw_z
# columns (see pseudocode in Step 42), making the mean-50 T-score the
# clinically interpretable anchor.

standardize_to_z <- function(x, type) {
  if (type == "tscore") return((x - 50) / 10)
  if (type == "ss")     return((x - 100) / 15)
  if (type == "zscore") return(x)
  if (type == "raw_z")  return(scale(x)[, 1])
  return(x)
}

std_cols <- paste0(score_meta$col, "_std")

cohort_std <- cohort %>%
  select(study_group, all_of(all_outcome_cols)) %>%
  mutate(across(
    all_of(all_outcome_cols),
    ~ standardize_to_z(.x, score_meta$score_type[score_meta$col == cur_column()]),
    .names = "{col}_std"
  )) %>%
  select(study_group, all_of(std_cols))

# Step 44: Pairwise t-tests across neuropsychological columns
# Welch two-sample t-test (var.equal = FALSE) used for all pairwise group
# comparisons. This does not assume equal variance between groups, which
# is appropriate for small and unequal group sizes (n = 23, 23, 20).
# All tests are computed once, then BH-FDR is applied two ways:
#   p_adj_BH      within each domain (cognitive, psychiatric); used for the
#                 significance brackets in Step 45
#   p_adj_BH_all  across the full set of tests
#
# Three pairwise contrasts per column:
#   Acute IPV-BI vs Acute Non-IPV-BI
#   Acute IPV-BI vs Chronic IPV-BI
#   Acute Non-IPV-BI vs Chronic IPV-BI
#
# Cohen's d is computed as a measure of clinical significance for each
# pairwise contrast, using the pooled SD:
#   d = (mean_g1 - mean_g2) / sd_pooled
#   sd_pooled = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1+n2-2))
#
# Effect size benchmarks (Cohen 1988):
#   |d| < 0.2  negligible
#   |d| >= 0.2 small
#   |d| >= 0.5 medium
#   |d| >= 0.8 large
#
# A result may be statistically significant but clinically negligible (small d),
# or miss significance but show a meaningful effect, particularly relevant
# at n = 66 where power is limited.

run_pairwise_ttests <- function(data, cols) {
  contrasts <- list(
    c("Acute IPV-BI", "Acute Non-IPV-BI"),
    c("Acute IPV-BI", "Chronic IPV-BI"),
    c("Acute Non-IPV-BI", "Chronic IPV-BI")
  )

  map_dfr(cols, function(col) {
    map_dfr(contrasts, function(pair) {
      g1 <- data[[col]][data$study_group == pair[1]]; g1 <- g1[!is.na(g1)]
      g2 <- data[[col]][data$study_group == pair[2]]; g2 <- g2[!is.na(g2)]
      if (length(g1) < 3 || length(g2) < 3) return(NULL)

      tt        <- t.test(g1, g2, var.equal = FALSE)
      sd_pooled <- sqrt(((length(g1) - 1) * var(g1) + (length(g2) - 1) * var(g2)) /
                          (length(g1) + length(g2) - 2))
      cohens_d  <- (mean(g1) - mean(g2)) / sd_pooled

      tibble(
        domain    = score_meta$domain[score_meta$col == col],
        col       = col,
        label     = score_meta$label[score_meta$col == col],
        contrast  = paste(pair[1], "vs", pair[2]),
        t_stat    = unname(tt$statistic),
        df        = unname(tt$parameter),
        p_value   = tt$p.value,
        mean_g1   = mean(g1),
        mean_g2   = mean(g2),
        mean_diff = mean(g1) - mean(g2),
        cohens_d  = cohens_d,
        d_label   = as.character(cut(abs(cohens_d), c(0, 0.2, 0.5, 0.8, Inf),
                                     labels = c("negligible", "small", "medium", "large"),
                                     right = FALSE))
      )
    })
  })
}

ttest_full <- run_pairwise_ttests(cohort, all_outcome_cols) %>%
  group_by(domain) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(p_adj_BH_all = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj_BH)

ttest_cognitive   <- filter(ttest_full, domain == "cognitive")
ttest_psychiatric <- filter(ttest_full, domain == "psychiatric")

ttest_full %>%
  group_by(domain) %>%
  summarise(fdr_sig_05 = sum(p_adj_BH < 0.05, na.rm = TRUE),
            effect_d_ge_0.5 = sum(abs(cohens_d) >= 0.5, na.rm = TRUE),
            .groups = "drop") %>%
  print()

# Step 45: Violin plots, neuropsychological scores by study group
# One panel per domain, with one violin per study group within each test.
# Points are individual participants colored by study group and black bars
# mark group means.
# Brackets are drawn only for pairwise comparisons with BH-adjusted p < 0.05
# from Step 44. Comparisons that are not significant have no bracket.

# Colors match group_colors in milestone 3 (Step 29).
group_palette <- c(
  "Acute IPV-BI"     = "#E41A1C",
  "Acute Non-IPV-BI" = "#377EB8",
  "Chronic IPV-BI"   = "#984EA3"
)
groups <- names(group_palette)

make_violin_panel <- function(std_data, cols, ttest_results, panel_title) {
  dodge     <- 0.8
  offsets   <- setNames((seq_along(groups) - (length(groups) + 1) / 2) * dodge / length(groups), groups)
  label_map <- setNames(score_meta$label[match(cols, score_meta$col)], paste0(cols, "_std"))

  plot_long <- std_data %>%
    select(study_group, all_of(names(label_map))) %>%
    pivot_longer(-study_group, names_to = "measure", values_to = "z_score") %>%
    filter(!is.na(z_score)) %>%
    mutate(measure_label = factor(unname(label_map[measure]), levels = unname(label_map)),
           study_group   = factor(study_group, levels = groups))

  brackets <- ttest_results %>%
    filter(p_adj_BH < 0.05) %>%
    mutate(measure_label = factor(unname(label_map[paste0(col, "_std")]), levels = unname(label_map)),
           g1    = str_remove(contrast, " vs .*"),
           g2    = str_remove(contrast, ".* vs "),
           stars = case_when(p_adj_BH < 0.001 ~ "***", p_adj_BH < 0.01 ~ "**", TRUE ~ "*")) %>%
    left_join(plot_long %>% group_by(measure_label) %>%
                summarise(y_top = max(z_score), .groups = "drop"),
              by = "measure_label") %>%
    group_by(measure_label) %>%
    mutate(y = y_top + 0.3 * row_number()) %>%
    ungroup() %>%
    mutate(x = as.integer(measure_label),
           xmin = x + unname(offsets[g1]),
           xmax = x + unname(offsets[g2]))

  segs <- bind_rows(
    brackets %>% transmute(x = xmin, xend = xmax, y = y,        yend = y),
    brackets %>% transmute(x = xmin, xend = xmin, y = y - 0.08, yend = y),
    brackets %>% transmute(x = xmax, xend = xmax, y = y - 0.08, yend = y)
  )

  ggplot(plot_long, aes(x = measure_label, y = z_score,
                        fill = study_group, color = study_group)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_violin(position = position_dodge(dodge), width = 0.75, trim = TRUE,
                alpha = 0.15, color = "grey45", linewidth = 0.35) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = dodge),
               size = 0.9, alpha = 0.85) +
    stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)),
                 geom = "errorbar", position = position_dodge(dodge),
                 width = 0.22, linewidth = 0.8, color = "black") +
    geom_segment(data = segs, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = "grey25", linewidth = 0.4) +
    geom_text(data = brackets, aes(x = (xmin + xmax) / 2, y = y, label = stars),
              inherit.aes = FALSE, vjust = 0.2, size = 4) +
    scale_fill_manual(values = group_palette, name = "Study group") +
    scale_color_manual(values = group_palette, name = "Study group") +
    labs(title = panel_title, x = NULL, y = "Standardized score (z)") +
    theme_classic(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold", hjust = 0.5, size = 10),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      axis.text.x        = element_text(angle = 25, hjust = 1, size = 8),
      legend.title       = element_text(face = "bold", size = 8),
      legend.text        = element_text(size = 7),
      legend.key.size    = unit(0.4, "cm"),
      legend.background  = element_rect(color = "grey80", linewidth = 0.3)
    )
}

violin_specs <- list(
  cognitive   = list(cols = cognitive_cols,   tt = ttest_cognitive,
                     title = "Cognitive Test Scores by Study Group",    width = 14),
  psychiatric = list(cols = psychiatric_cols, tt = ttest_psychiatric,
                     title = "Psychiatric Symptom Scores by Study Group", width = 8)
)

iwalk(violin_specs, function(spec, domain) {
  p <- make_violin_panel(cohort_std, spec$cols, spec$tt, spec$title)
  print(p)
  ggsave(file.path(output_dir, paste0("violin_", domain, ".png")),
         plot = p, width = spec$width, height = 6, dpi = 300)
})


# Step 46: Shapiro-Wilk normality test for correlation method assignment
# One Shapiro-Wilk p-value per variable (proteins, neuropsychological
# outcomes, and exposure counts), computed on the full cohort (n = 66).
# For each protein and outcome pair: Pearson if both variables pass
# (p > 0.05), Spearman if either fails. The same assignment is used in
# Steps 47 and 48, including the within-group runs.
#
# Note: at n = 66, Shapiro-Wilk has moderate power. Spearman is a
# conservative fallback and is appropriate for mixed score types.

npq_cols <- cohort %>% select(starts_with("NPQ_"), -any_of(qc_meta_cols)) %>% names()
npq_cols <- npq_cols[str_remove(npq_cols, "^NPQ_") != "MSLN"]   # exclude MSLN per milestone 3

analysis_meta <- bind_rows(
  score_meta %>% select(col, label, domain),
  tibble(col = exposure_cols, label = c("BISA", "HELPS"), domain = "exposure")
)

label_outcomes <- function(df) {
  left_join(df, analysis_meta %>% rename(outcome = col, outcome_label = label), by = "outcome")
}

sw_p <- sapply(c(npq_cols, analysis_meta$col), function(col) {
  x <- cohort[[col]][!is.na(cohort[[col]])]
  if (length(x) < 3) NA_real_ else shapiro.test(x)$p.value
})

pick_method <- function(a, b) if (isTRUE(sw_p[a] > 0.05 && sw_p[b] > 0.05)) "pearson" else "spearman"

# Step 47: Biomarker correlations, full cohort
# Every protein is correlated with every outcome once (cognitive, psychiatric,
# and the BISA/HELPS exposure counts). BH-FDR is then applied within each
# family of tests, so the same correlations are never recomputed:
#   full_neuropsych  cognitive + psychiatric   221 x 13 = 2873 tests
#   cognitive        cognitive only            221 x 10 = 2210 tests
#   psychiatric      PCL-5, PHQ-9, NSI         221 x 3  =  663 tests
#   exposure         BISA, HELPS               221 x 2  =  442 tests
# Exposure counts are a dose-response question (does lifetime IPV-BI burden
# track protein levels) distinct from the cognition associations.
# Separating cognitive and psychiatric families shows whether biomarker
# associations are domain specific. Psychiatric has only 3 outcomes, so
# its results are exploratory.

run_correlations <- function(data, outcomes) {
  grid <- expand_grid(protein = npq_cols, outcome = outcomes)

  map2_dfr(grid$protein, grid$outcome, function(prot, out) {
    x  <- data[[prot]]
    y  <- data[[out]]
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 5) return(NULL)

    m  <- pick_method(prot, out)
    ct <- tryCatch(suppressWarnings(cor.test(x[ok], y[ok], method = m)),
                   error = function(e) NULL)
    if (is.null(ct)) return(NULL)

    tibble(protein = prot, outcome = out, method = m,
           r = unname(ct$estimate), p_value = ct$p.value, n = sum(ok))
  }) %>%
    label_outcomes()
}

corr_families <- list(
  full_neuropsych = c("cognitive", "psychiatric"),
  cognitive       = "cognitive",
  psychiatric     = "psychiatric",
  exposure        = "exposure"
)

adjust_families <- function(raw, families, group_label) {
  imap_dfr(families, function(domains, family) {
    raw %>%
      filter(domain %in% domains) %>%
      mutate(group = group_label, family = family,
             p_adj_BH = p.adjust(p_value, method = "BH")) %>%
      arrange(p_adj_BH)
  })
}

corr_all <- adjust_families(run_correlations(cohort, analysis_meta$col),
                            corr_families, "full_cohort")

# Step 48: Within-group supplementary correlations
# The neuropsychological families from Step 47 (without exposure) are
# re-run separately within each study group, using the same method
# assignment. Group sizes are n = 23, 23, 20, so power is very limited and
# results are hypothesis generating only. BH-FDR is applied within each
# group and family.

corr_within <- map_dfr(groups, function(g) {
  adjust_families(run_correlations(filter(cohort, study_group == g), all_outcome_cols),
                  corr_families[1:3], g)
})

corr_results <- bind_rows(corr_all, corr_within)

corr_results %>%
  group_by(group, family) %>%
  summarise(fdr_05 = sum(p_adj_BH < 0.05, na.rm = TRUE),
            fdr_20 = sum(p_adj_BH < 0.20, na.rm = TRUE),
            .groups = "drop") %>%
  print(n = Inf)

# Step 49: Top ten biomarker correlations
# The ten strongest protein and outcome pairs from the Step 47 correlations
# across all 13 neuropsychological scores, ranked by BH-adjusted p-value.
# r is Pearson or Spearman depending on the method assigned in Step 46.
# The table is printed here and saved as the Top 10 Correlations sheet in
# Step 55.

top_corr <- corr_all %>%
  filter(family == "full_neuropsych") %>%
  arrange(p_adj_BH, p_value) %>%
  slice_head(n = 10) %>%
  transmute(Rank = row_number(), Protein = str_remove(protein, "^NPQ_"),
            Outcome = outcome_label, Method = str_to_title(method),
            r = round(r, 3), n = n, p = signif(p_value, 3), FDR = signif(p_adj_BH, 3))

print(top_corr)

# Step 50: Correlation heatmaps, cognitive and psychiatric panels
# Rows: proteins with at least one association at FDR < 0.20 in the
# full-cohort family. Columns: outcomes, in the order used in Step 45.
# Cell color: Pearson or Spearman r. Cells marked * where FDR < 0.05.
# Saved only, as in milestone 3. Each heatmap is built silently and drawn to
# its own png device, so it does not depend on which plot device is open.
# to_wide() and save_pheatmap() are also used for the regression heatmap in
# Step 52.

to_wide <- function(df, value) {
  df %>%
    mutate(protein = str_remove(protein, "^NPQ_"),
           outcome_label = factor(outcome_label, levels = score_meta$label)) %>%
    select(protein, outcome_label, all_of(value)) %>%
    pivot_wider(names_from = outcome_label, values_from = all_of(value), names_sort = TRUE) %>%
    column_to_rownames("protein") %>%
    as.matrix()
}

save_pheatmap <- function(mat, marks, panel_title, file_name, limit = 1,
                          row_font = 7, row_height = 0.25) {
  ph <- pheatmap::pheatmap(
    mat,
    display_numbers = marks,
    number_color    = "black",
    fontsize_number = 9,
    color           = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
    breaks          = seq(-limit, limit, length.out = 101),
    cluster_rows    = nrow(mat) > 1,
    cluster_cols    = FALSE,
    fontsize_row    = row_font,
    fontsize_col    = 8,
    fontsize        = 9,
    main            = panel_title,
    border_color    = NA,
    na_col          = "grey90",
    silent          = TRUE
  )

  png(file.path(output_dir, file_name),
      width = max(8, ncol(mat) * 0.9 + 3), height = max(6, nrow(mat) * row_height + 2),
      units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  invisible(dev.off())
}

make_corr_heatmap <- function(fam, panel_title, file_name) {
  plot_data <- corr_all %>% filter(family == fam)
  keep <- plot_data %>% filter(p_adj_BH < 0.20) %>% pull(protein) %>% unique()
  if (length(keep) == 0) {
    cat("No proteins at FDR < 0.20 for", fam, "so heatmap skipped\n")
    return(invisible(NULL))
  }

  plot_data <- filter(plot_data, protein %in% keep)
  p_mat <- to_wide(plot_data, "p_adj_BH")
  save_pheatmap(to_wide(plot_data, "r"), ifelse(!is.na(p_mat) & p_mat < 0.05, "*", ""),
                panel_title, file_name)
}

tribble(
  ~fam,          ~panel_title,                                                                          ~file_name,
  "cognitive",   "Biomarker-Cognitive Correlations (FDR < 0.20 proteins shown; * = FDR < 0.05)",   "heatmap_biomarker_cognitive.png",
  "psychiatric", "Biomarker-Psychiatric Correlations (FDR < 0.20 proteins shown; * = FDR < 0.05)", "heatmap_biomarker_psychiatric.png"
) %>%
  pwalk(make_corr_heatmap)

# Step 51: Linear regression, neuropsychological outcomes predicted by biomarkers
# For every protein and every outcome, a separate linear regression is fit:
#   outcome_z ~ NPQ_protein + age_years + sex
# All 221 proteins are modelled (2873 models), matching the correlation
# scope in Step 47, and BH-FDR is applied across all 2873 models.
#
# Outcomes are the standardized z-scores from Step 43, so betas are on a
# common scale: change in outcome SD per 1 log2 NPQ unit. Because each
# outcome is only rescaled, p-values are the same as for the raw scores.
# Age and sex are included as in milestone 3 (Step 22). Unlike the limma
# design, study_group is not a covariate here, so associations are pooled
# across the three groups and may partly reflect group differences.
# sex is encoded as numeric (Male = 1, Female = 0).
#
# r_squared reflects the full model (protein + age + sex), not the protein
# alone, and is reported for descriptive completeness only.

cohort_lm <- bind_cols(cohort, select(cohort_std, -study_group)) %>%
  mutate(sex_num = ifelse(sex == "Male", 1L, 0L))

run_lm_per_protein <- function(data, outcome_col, proteins = npq_cols,
                               covariates = c("age_years", "sex_num")) {
  map_dfr(proteins, function(prot) {
    df <- data %>%
      select(y = all_of(paste0(outcome_col, "_std")), x = all_of(prot), all_of(covariates)) %>%
      drop_na()
    if (nrow(df) < 10) return(NULL)

    covs <- covariates[map_lgl(covariates, ~ n_distinct(df[[.x]]) > 1)]
    s  <- summary(lm(reformulate(c("x", covs), response = "y"), data = df))
    cf <- coef(s)
    if (!"x" %in% rownames(cf)) return(NULL)

    tibble(protein = prot, outcome = outcome_col,
           beta = cf["x", "Estimate"], se = cf["x", "Std. Error"],
           t_stat = cf["x", "t value"], p_value = cf["x", "Pr(>|t|)"],
           r_squared = s$r.squared, adj_r_sq = s$adj.r.squared, n = nrow(df))
  })
}

# adjusted_line() gives the regression line and 95% confidence band drawn in
# the scatter plots of Steps 53 and 54. Covariates are held at their sample
# means (group proportions for study group), so the line has the same slope
# as the reported beta and runs through the center of the data. Points in the
# plots stay on the raw scale.
adjusted_line <- function(df, covariates) {
  fit  <- lm(reformulate(c("x", covariates), response = "y"), data = df)
  mm   <- model.matrix(fit)
  grid <- seq(min(df$x), max(df$x), length.out = 80)
  X0   <- matrix(colMeans(mm), nrow = length(grid), ncol = ncol(mm), byrow = TRUE,
                 dimnames = list(NULL, colnames(mm)))
  X0[, "x"] <- grid
  est  <- drop(X0 %*% coef(fit))
  se   <- sqrt(rowSums((X0 %*% vcov(fit)) * X0))
  crit <- qt(0.975, df.residual(fit))
  tibble(x = grid, y = est, ymin = est - crit * se, ymax = est + crit * se)
}

lm_results <- map_dfr(all_outcome_cols, ~ run_lm_per_protein(cohort_lm, .x)) %>%
  label_outcomes() %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj_BH)

# Supplementary: the same models run separately inside each study group
# (n = 20 to 23, adjusted for age and sex; a covariate with no variation in a
# group is dropped). BH-FDR is applied across the 2873 models in each group.
# Used for the group volcano plots in Step 52. Hypothesis generating only.
lm_within <- map_dfr(groups, function(g) {
  data_g <- filter(cohort_lm, study_group == g)
  map_dfr(all_outcome_cols, ~ run_lm_per_protein(data_g, .x)) %>%
    label_outcomes() %>%
    mutate(group = g, p_adj_BH = p.adjust(p_value, method = "BH"))
})

lm_results %>%
  summarise(models = n(),
            fdr_05 = sum(p_adj_BH < 0.05, na.rm = TRUE),
            fdr_20 = sum(p_adj_BH < 0.20, na.rm = TRUE)) %>%
  print()

# Step 52: All-protein regression visuals
# The 2873 models from Step 51 are too many for individual scatter plots, so
# they are shown two ways. Scatter plots are used for specific pairs in
# Steps 53 and 54.
#   Heatmap  beta for every protein (rows) and outcome (columns), in SD of the
#            outcome per 1 log2 NPQ unit. Cells are marked * where FDR < 0.05
#            and + where FDR < 0.20.
#   Volcano  one panel per outcome, beta against raw -log10 p, with the two
#            smallest p-values in each panel labelled. Points with raw p < 0.05
#            are colored. The grey dashed line is raw p = 0.05 and the black
#            dashed line is the raw p at the FDR 0.20 cutoff, drawn only if a
#            model reaches it. One volcano is made for all participants (dark
#            grey) and one for each study group using the within-group models
#            from Step 51, in that group's milestone 3 color.
# Raw p is on the y axis because FDR is calculated across all 2873 models
# and would flatten every panel.

beta_mat <- to_wide(lm_results, "beta")
fdr_mat  <- to_wide(lm_results, "p_adj_BH")

save_pheatmap(
  beta_mat,
  marks = ifelse(is.na(fdr_mat), "", ifelse(fdr_mat < 0.05, "*", ifelse(fdr_mat < 0.20, "+", ""))),
  panel_title = "Regression beta, all proteins by outcome (* FDR < 0.05, + FDR < 0.20)",
  file_name   = "heatmap_regression_all_proteins.png",
  limit       = ceiling(max(abs(beta_mat), na.rm = TRUE) * 10) / 10,
  row_font    = 5, row_height = 0.15
)

make_regression_volcano <- function(df, title, hit_color) {
  fdr20_p <- if (any(df$p_adj_BH < 0.20)) max(df$p_value[df$p_adj_BH < 0.20]) else NA_real_
  df <- df %>%
    mutate(protein = str_remove(protein, "^NPQ_"),
           outcome_label = factor(outcome_label, levels = score_meta$label))

  ggplot(df, aes(x = beta, y = -log10(p_value))) +
    geom_point(aes(color = p_value < 0.05), alpha = 0.7, size = 0.9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_hline(data = tibble(y = -log10(fdr20_p)) %>% filter(!is.na(y)), aes(yintercept = y),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_text_repel(data = df %>% group_by(outcome) %>%
                      slice_min(p_value, n = 2, with_ties = FALSE) %>% ungroup(),
                    aes(label = protein), size = 2.2, max.overlaps = 20,
                    segment.color = "grey50", box.padding = 0.3) +
    scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = hit_color), guide = "none") +
    facet_wrap(~ outcome_label, ncol = 4) +
    labs(title    = title,
         subtitle = "Beta adjusted for age and sex; colored points raw p < 0.05; grey line raw p = 0.05; black line raw p at the FDR 0.20 cutoff",
         x = "Beta (SD of outcome per 1 log2 NPQ unit)", y = expression(-log[10](raw~p))) +
    theme_classic(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 10),
          plot.subtitle = element_text(color = "grey40", hjust = 0.5, size = 7),
          strip.text    = element_text(face = "bold", size = 8))
}

c(list(all = lm_results), split(lm_within, lm_within$group)) %>%
  iwalk(function(df, set) {
    is_all <- set == "all"
    p <- make_regression_volcano(
      df,
      title     = if (is_all) "All-Protein Regression by Outcome" else paste0(set, ": Within-Group Regression by Outcome"),
      hit_color = if (is_all) "grey25" else unname(group_palette[set])
    )
    ggsave(file.path(output_dir, paste0("volcano_regression_", str_replace_all(set, "[^A-Za-z0-9]", ""), ".png")),
           plot = p, width = 12, height = 9, dpi = 300)
  })

# Step 53: Targeted regression, top milestone 3 proteins
# The three most significant proteins in the primary milestone 3 contrast
# (Acute IPV-BI vs Acute Non-IPV-BI, by limma BH-adjusted p) are tested
# against three outcomes chosen to span the neuropsychological data:
#   RAVLT Delayed Recall  the test with FDR-significant group differences in Step 44
#   COWA                  verbal fluency, a second cognitive domain
#   PCL-5                 PTSD symptoms, the psychiatric domain
# This is 3 x 3 = 9 tests per model, with BH-FDR across those 9 only.
#
# Two models are fit so the effect of study group is visible:
#   age + sex                  pooled across groups, as in Step 51
#   age + sex + study group    the primary model for this step
# Study group matters here because these proteins differ between groups and
# so do several test scores, so a pooled association could come from the
# group difference alone. The group model assumes one common slope across
# groups. Each panel of the figure shows that adjusted line with a 95%
# confidence ribbon, so the line matches the beta in the panel title.
# To select proteins from another contrast, change target_contrast.

target_contrast <- "Acute_IPV_vs_NonIPV"
if (!file.exists("pairwise_results.rds")) stop("pairwise_results.rds not found, run milestone 3 first.")

target_tbl <- readRDS("pairwise_results.rds")[[target_contrast]] %>%
  slice_min(limma_adj_p, n = 3, with_ties = FALSE) %>%
  select(protein, limma_logFC, limma_adj_p)
print(target_tbl)

target_proteins <- paste0("NPQ_", target_tbl$protein)
target_outcomes <- c("RAVLT_delayed_recall", "COWA_total_words", "PCL5_total")
target_models   <- list("age + sex"               = c("age_years", "sex_num"),
                        "age + sex + study group" = c("age_years", "sex_num", "study_group"))
primary_model   <- names(target_models)[2]

lm_targeted <- imap_dfr(target_models, function(covs, model) {
  map_dfr(target_outcomes, ~ run_lm_per_protein(cohort_lm, .x, target_proteins, covs)) %>%
    mutate(model = model)
}) %>%
  label_outcomes() %>%
  group_by(model) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(model, p_adj_BH)

lm_targeted %>%
  filter(model == primary_model) %>%
  transmute(protein = str_remove(protein, "^NPQ_"), outcome_label,
            beta = round(beta, 3), p = signif(p_value, 3), FDR = signif(p_adj_BH, 3)) %>%
  print()

fig <- lm_targeted %>%
  filter(model == primary_model) %>%
  arrange(match(protein, target_proteins), match(outcome, target_outcomes)) %>%
  mutate(strip = paste0(str_remove(protein, "^NPQ_"), " and ", outcome_label,
                        "\nbeta = ", round(beta, 2), ", FDR = ", signif(p_adj_BH, 2)))

fig_data <- pmap_dfr(list(fig$protein, fig$outcome, fig$strip), function(prot, out, strip) {
  tibble(study_group = cohort_lm$study_group, age_years = cohort_lm$age_years,
         sex_num = cohort_lm$sex_num, x = cohort_lm[[prot]],
         y = cohort_lm[[paste0(out, "_std")]], strip = strip)
}) %>%
  drop_na(x, y, age_years, sex_num) %>%
  mutate(strip = factor(strip, levels = fig$strip))

fig_lines <- fig_data %>%
  group_by(strip) %>%
  group_modify(~ adjusted_line(.x, target_models[[primary_model]])) %>%
  ungroup()

p_targeted <- ggplot(fig_data, aes(x = x, y = y)) +
  geom_ribbon(data = fig_lines, aes(x = x, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = "grey85", alpha = 0.4) +
  geom_line(data = fig_lines, aes(x = x, y = y), inherit.aes = FALSE,
            color = "grey40", linewidth = 0.7) +
  geom_point(aes(color = study_group), alpha = 0.75, size = 1.2) +
  facet_wrap(~ strip, ncol = 3, scales = "free") +
  scale_color_manual(values = group_palette, name = "Study group") +
  labs(title    = "Top milestone 3 proteins and selected neuropsychological scores",
       subtitle = "Regression line adjusted for age, sex and study group with 95% confidence ribbon; FDR across the 9 targeted tests",
       x = "Protein level (NPQ log2)", y = "Standardized score (z)") +
  theme_classic(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 10),
        plot.subtitle = element_text(color = "grey40", hjust = 0.5, size = 7),
        strip.text    = element_text(face = "bold", size = 7),
        legend.title  = element_text(face = "bold", size = 8),
        legend.text   = element_text(size = 7))

ggsave(file.path(output_dir, "scatter_targeted.png"),
       plot = p_targeted, width = 11, height = 9, dpi = 300)

# Step 54: Scatter plots, protein vs neuropsychological outcome
# One scatter plot per protein and outcome pair with FDR < 0.20 in Step 51
# (the all-protein regression; Step 53 has its own figure).
# Each plot shows individual participants colored by study group and the
# regression line from Step 51, adjusted for age and sex, with a 95%
# confidence ribbon. The line has the same slope as the reported beta.
#
# If no pair reaches FDR < 0.20, the 5 pairs with the smallest raw p-value
# are plotted instead and labelled as exploratory with the unadjusted p.
# Plots are saved individually and as one multi-panel PDF per domain.

make_scatter <- function(row, data, exploratory) {
  df <- data %>%
    select(study_group, age_years, sex_num,
           x = all_of(row$protein), y = all_of(paste0(row$outcome, "_std"))) %>%
    drop_na()
  line <- adjusted_line(df, c("age_years", "sex_num"))
  protein_label <- str_remove(row$protein, "^NPQ_")
  p_text <- if (exploratory) paste0("exploratory, unadjusted p = ", signif(row$p_value, 2))
            else paste0("FDR p = ", signif(row$p_adj_BH, 2))

  ggplot(df, aes(x = x, y = y)) +
    geom_ribbon(data = line, aes(x = x, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = "grey85", alpha = 0.4) +
    geom_line(data = line, aes(x = x, y = y), inherit.aes = FALSE,
              color = "grey40", linewidth = 0.7) +
    geom_point(aes(color = study_group), alpha = 0.75, size = 2.0) +
    scale_color_manual(values = group_palette, name = "Study group") +
    labs(title    = paste0(protein_label, " vs ", row$outcome_label),
         subtitle = paste0("beta = ", round(row$beta, 3),
                           ", adj R-squared = ", round(row$adj_r_sq, 3), ", ", p_text),
         x = paste0(protein_label, " (NPQ log2)"),
         y = paste0(row$outcome_label, " (z-score)")) +
    theme_classic(base_size = 10) +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 10),
          plot.subtitle   = element_text(color = "grey40", hjust = 0.5, size = 7),
          legend.text     = element_text(size = 7),
          legend.title    = element_text(face = "bold", size = 8),
          legend.key.size = unit(0.35, "cm"))
}

scatter_hits <- lm_results %>% filter(p_adj_BH < 0.20)
exploratory_scatter <- nrow(scatter_hits) == 0
if (exploratory_scatter) scatter_hits <- slice_min(lm_results, p_value, n = 5, with_ties = FALSE)

scatter_plots <- scatter_hits %>%
  split(seq_len(nrow(scatter_hits))) %>%
  map(make_scatter, data = cohort_lm, exploratory = exploratory_scatter)

iwalk(scatter_plots, function(p, i) {
  hit   <- scatter_hits[as.integer(i), ]
  fname <- paste0("scatter_", str_replace_all(hit$protein, "[^A-Za-z0-9]", ""), "_",
                  str_replace_all(hit$outcome, "[^A-Za-z0-9]", ""), ".png")
  ggsave(file.path(output_dir, fname), plot = p, width = 6, height = 5, dpi = 300)
})

iwalk(split(scatter_plots, scatter_hits$domain), function(plots, dom) {
  n_col <- min(3L, length(plots))
  n_row <- ceiling(length(plots) / n_col)
  ggsave(file.path(output_dir, paste0("scatter_panel_", dom, ".pdf")),
         ggarrange(plotlist = unname(plots), ncol = n_col, nrow = n_row,
                   common.legend = TRUE, legend = "right"),
         width = n_col * 6, height = n_row * 5)
})

# Step 55: Export results and save updated RDS
# One Excel workbook (milestone4_results.xlsx) with one sheet per result set.
# Every sheet is sorted by FDR, has a frozen header row, and has filters on.
#   Notes                      what each sheet contains and how to read it
#   Top 10 Correlations        ten strongest protein and outcome pairs, all-protein correlations
#   Cognitive Group Tests      Welch t-tests and Cohen's d, cognitive scores
#   Psychiatric Group Tests    Welch t-tests and Cohen's d, symptom scales
#   Biomarker All Neuropsych   correlations, cognitive and psychiatric together
#   Biomarker Cognitive        correlations, cognitive scores only
#   Biomarker Psychiatric      correlations, symptom scales only
#   Biomarker Exposure         correlations with BISA and HELPS
#   Within Group Correlations  supplementary correlations in each study group
#   Regression Cognitive       linear regression, cognitive scores
#   Regression Psychiatric     linear regression, symptom scales
#   Within Group Regression    the all-protein regressions run inside each study group
#   Targeted Regression        top 3 milestone 3 proteins against RAVLT Delayed Recall, COWA, PCL-5
#
# The cohort RDS is updated to include the TOMM validity flag column.

tidy_ttests <- function(df) {
  df %>% transmute(
    Test = label,
    `Group 1` = str_remove(contrast, " vs .*"),
    `Group 2` = str_remove(contrast, ".* vs "),
    `Mean group 1` = round(mean_g1, 2),
    `Mean group 2` = round(mean_g2, 2),
    `Mean difference` = round(mean_diff, 2),
    `Cohen's d` = round(cohens_d, 2),
    `Effect size` = d_label,
    t = round(t_stat, 2),
    df = round(df, 1),
    p = signif(p_value, 3),
    `FDR within domain` = signif(p_adj_BH, 3),
    `FDR all tests` = signif(p_adj_BH_all, 3)
  )
}

tidy_corr <- function(df, with_group = FALSE) {
  out <- df %>% transmute(
    Group = group, Family = family,
    Protein = str_remove(protein, "^NPQ_"),
    Outcome = outcome_label,
    Method = str_to_title(method),
    r = round(r, 3), n = n,
    p = signif(p_value, 3),
    FDR = signif(p_adj_BH, 3)
  )
  if (with_group) out else select(out, -Group, -Family)
}

tidy_lm <- function(df) {
  df %>% transmute(
    Protein = str_remove(protein, "^NPQ_"),
    Outcome = outcome_label,
    Beta = round(beta, 3), SE = round(se, 3), t = round(t_stat, 2),
    p = signif(p_value, 3), FDR = signif(p_adj_BH, 3),
    `R-squared` = round(r_squared, 3), `Adj R-squared` = round(adj_r_sq, 3), n = n
  )
}

corr_family <- function(fam) filter(corr_all, family == fam) %>% tidy_corr()

notes <- tribble(
  ~Item, ~Description,
  "Group tests", "Welch t-tests for the three group pairs. Means are in original score units. FDR within domain is BH across the tests in that sheet and is used for the plot brackets. FDR all tests is BH across both domains.",
  "Cohen's d", "Effect size for the group difference. Under 0.2 negligible, 0.2 small, 0.5 medium, 0.8 large.",
  "Top 10 Correlations", "The ten strongest protein and outcome pairs from the all-protein correlations across the 13 scores, ranked by FDR.",
  "Biomarker sheets", "Each protein correlated with each outcome across all 66 participants. Method is Pearson if both variables pass Shapiro-Wilk, otherwise Spearman. FDR is BH within the sheet.",
  "Within Group Correlations", "The same correlations run inside each study group (n = 20 to 23), with FDR applied within each group and family. Hypothesis generating only.",
  "Regression sheets", "One model per protein and outcome: outcome z-score ~ protein + age + sex. Beta is change in outcome SD per 1 log2 NPQ unit. FDR is BH across all models in both regression sheets together. Study group is not a covariate.",
  "Within Group Regression", "The all-protein regressions run separately in each study group (n = 20 to 23), adjusted for age and sex, with FDR across the 2873 models in each group. Hypothesis generating only.",
  "Targeted Regression", "The three most significant proteins from the primary milestone 3 contrast against RAVLT Delayed Recall, COWA and PCL-5. Two models are shown, with and without study group. FDR is BH across the 9 tests within each model.",
  "Direction", "Higher is better for cognitive scores. Higher is worse for PCL-5, PHQ-9 and NSI.",
  "Scores", "Raw scores are standardized within this sample (see Step 42 limitation note)."
)

result_sheets <- list(
  "Notes"                     = notes,
  "Top 10 Correlations"       = top_corr,
  "Cognitive Group Tests"     = tidy_ttests(ttest_cognitive),
  "Psychiatric Group Tests"   = tidy_ttests(ttest_psychiatric),
  "Biomarker All Neuropsych"  = corr_family("full_neuropsych"),
  "Biomarker Cognitive"       = corr_family("cognitive"),
  "Biomarker Psychiatric"     = corr_family("psychiatric"),
  "Biomarker Exposure"        = corr_family("exposure"),
  "Within Group Correlations" = tidy_corr(corr_within, with_group = TRUE),
  "Regression Cognitive"      = tidy_lm(filter(lm_results, domain == "cognitive")),
  "Regression Psychiatric"    = tidy_lm(filter(lm_results, domain == "psychiatric")),
  "Within Group Regression"   = tidy_lm(lm_within) %>% mutate(Group = lm_within$group, .before = 1),
  "Targeted Regression"       = tidy_lm(lm_targeted) %>% mutate(Model = lm_targeted$model, .before = 1)
)

wb <- createWorkbook()
iwalk(result_sheets, function(df, name) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = createStyle(textDecoration = "bold", border = "bottom"))
  freezePane(wb, name, firstRow = TRUE)
  addFilter(wb, name, rows = 1, cols = seq_len(ncol(df)))
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
})
setColWidths(wb, "Notes", cols = 1:2, widths = c(28, 120))
addStyle(wb, "Notes", createStyle(wrapText = TRUE, valign = "top"),
         rows = 2:(nrow(notes) + 1), cols = 1:2, gridExpand = TRUE)
saveWorkbook(wb, file.path(output_dir, "milestone4_results.xlsx"), overwrite = TRUE)

saveRDS(cohort, file = "ipv_cohort_clean.rds")
cat("Milestone 4 outputs written to:", output_dir, "\n")
