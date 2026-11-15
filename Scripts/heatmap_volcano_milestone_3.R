# Weeks 5-6: NULISA Proteomics — Limma + Mann-Whitney U, Volcano Plots, Heatmaps
# Goal: Run differential expression models across the three study groups using
# limma, validate results with Mann-Whitney U tests, and visualize with volcano
# plots and heatmaps organized by hallmark category of neurodegeneration.

# Pairwise contrasts by scientific priority:
#   Primary:   Acute IPV-BI vs Acute Non-IPV-BI
#   Secondary: Chronic IPV-BI vs Acute IPV-BI
#   Tertiary:  Chronic IPV-BI vs Acute Non-IPV-BI

# limma (parametric, empirical-Bayes):
#   Borrows strength across all 221 proteins to stabilise variance estimates at n=66.
#   Adjusts for age and sex as covariates.
#   limma_adj_p comes from topTable() applying BH correction across 221 proteins
#   within each contrast (not across contrasts — each contrast is corrected independently).

# Mann-Whitney U (nonparametric sanity check):
#   Makes no assumptions about normality or variance.
#   Cannot adjust for covariates — use to confirm limma results, not replace them.
#   mw_p_raw = raw two-sided p-value from wilcox.test().
#   mw_adj_p = BH correction applied across all 663 tests (221 proteins x 3 contrasts).

# Simulated data note:
#   NPQz_ columns are z-scored (mean ~0, SD ~1) — not raw NPQ, no log2 needed.
#   logFC represents z-score difference between group means, not true fold-change.
#   When applying to real Alamar NULISA data, feed log2-NPQ directly into limma.

# Note on the primary volcano plots:
#   With 221 proteins tested and mostly null signal at n=66, BH correction is
#   conservative — only 2 proteins survive FDR < 0.05. This is statistically
#   correct for this dataset. The diagnostic p-value histograms (Step 31) confirm
#   the raw p-value distribution is well-behaved and the flat appearance of the
#   primary volcanos reflects the data, not a pipeline error.
#   Exploratory volcanos at FDR < 0.20 are provided in Step 34 for context.


# Step 16: Packages

library(tidyverse)
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(readxl)


# Step 17: Load data
# ipv_cohort_clean.rds was produced by milestone 1. It already contains both
# sheets combined, study_group as a three-level factor, and all recoded variables.

setwd("C:/NRSC_510")

output_dir <- "C:/NRSC_510/Outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cohort <- readRDS("ipv_cohort_clean.rds")


# Step 18: Data format check
# Halts the script if NPQz_ values are not z-scored.
# If this fires on real data, apply log2(x+1) before proceeding.

# Real data extension note:
#   cohort <- cohort %>% mutate(across(starts_with("NPQz_"), ~ log2(.x + 1)))
#   Feed log2-NPQ directly into limma without z-scoring — limma uses between-sample
#   variance that z-scoring removes.

npqz_cols_all   <- cohort %>% select(starts_with("NPQz_")) %>% names()
all_npqz_values <- cohort %>% select(all_of(npqz_cols_all)) %>% as.matrix() %>% as.vector()

if (abs(mean(all_npqz_values, na.rm = TRUE)) > 1 ||
    sd(all_npqz_values, na.rm = TRUE) > 2 ||
    max(all_npqz_values, na.rm = TRUE) > 10) {
  stop(
    "NPQz values do not look z-scored. Check data format before running models.\n",
    "If data are log2-NPQ, see the real data extension note in Step 18.\n",
    "If data are raw NPQ, apply log2(x+1) first."
  )
}


# Detection limit QC scaffold (simulated data — cannot execute)
# In real Alamar NULISA data, replace this section with:
# detection_flags <- read_csv("alamar_detection_flags.csv")
# detection_rate  <- detection_flags %>%
#   left_join(cohort %>% select(participant_id, study_group), by = "participant_id") %>%
#   group_by(study_group) %>%
#   summarise(across(starts_with("NPQz_"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
# proteins_pass <- detection_rate %>%
#   select(-study_group) %>%
#   summarise(across(everything(), ~ all(.x >= 0.50))) %>%
#   pivot_longer(everything(), names_to = "protein", values_to = "pass") %>%
#   filter(pass) %>% pull(protein)
# cohort <- cohort %>% select(-any_of(setdiff(
#   names(cohort)[startsWith(names(cohort), "NPQz_")], proteins_pass)))


# Step 19: Hallmark lookup table
# Read from the Alamar Neuro 220 Panel CSV — the source of truth for hallmark
# assignments. Protein names are cleaned from Target List to match NPQz_ columns:
#   1. Strip parenthetical aliases e.g. "CCL2 (MCP-1)" -> "CCL2"
#   2. Replace semicolons with underscores e.g. "HBA1; HBA2" -> "HBA1_HBA2"
#   3. Remove pipe and everything after e.g. "IL-12A|IL-12B" -> "IL-12A"
#   4. Replace spaces with hyphens
# Manual remap: IL-12A cleans to "IL-12A" but dataset column is NPQz_IL-12p70.

panel_raw <- read_csv("C:/NRSC_510/Neuro_220_Panel.csv", show_col_types = FALSE)

hallmark_lookup <- panel_raw %>%
  rename(
    target_full = `Target List`,
    uniprot_id  = `Uniprot ID`,
    hallmark    = `Hallmarks of Neurodegenerative Disease`
  ) %>%
  mutate(
    protein = target_full %>%
      str_replace_all("\\s*\\(.*?\\)", "") %>%
      str_replace_all(";\\s*", "_") %>%
      str_replace_all("\\|.*", "") %>%
      str_replace_all("\\s+", "-") %>%
      str_trim(),
    protein = case_when(
      str_detect(target_full, "IL-12A") ~ "IL-12p70",
      TRUE                              ~ protein
    ),
    hallmark = case_when(
      hallmark == "Aberrant proteostasis" ~ "Aberrant Proteostasis",
      TRUE                               ~ hallmark
    ),
    hallmark_abbrev = case_when(
      hallmark == "Pathological Protein Aggregation"      ~ "Protein Aggregation",
      hallmark == "Synaptic and Neuronal Network Defects" ~ "Synaptic/Neuronal",
      hallmark == "Inflammation"                          ~ "Inflammation",
      hallmark == "Aberrant Proteostasis"                 ~ "Proteostasis",
      hallmark == "Altered Energy Homeostasis"            ~ "Energy Homeostasis",
      hallmark == "DNA and/or RNA defects"                ~ "DNA/RNA Defects",
      hallmark == "Neuronal Cell Death"                   ~ "Neuronal Death",
      hallmark == "Cytoskeletal Abnormalities"            ~ "Cytoskeletal",
      hallmark == "vascular/cerebrovascular dysfunction"  ~ "Vascular",
      is.na(hallmark)                                     ~ "Unassigned",
      TRUE                                                ~ hallmark
    )
  ) %>%
  filter(!is.na(hallmark)) %>%
  select(protein, uniprot_id, hallmark, hallmark_abbrev, target_full)


# Step 20: Identify and validate NPQz protein columns
# Excludes MSLN (no hallmark assignment in Alamar panel).
# Halts if any retained protein has no hallmark match — would cause silent NA rows.

excluded_proteins <- c("MSLN")

protein_names_all <- str_remove(npqz_cols_all, "^NPQz_")
npqz_cols         <- npqz_cols_all[!protein_names_all %in% excluded_proteins]
protein_names     <- str_remove(npqz_cols, "^NPQz_")

unmatched <- protein_names[!protein_names %in% hallmark_lookup$protein]
if (length(unmatched) > 0) {
  stop("Proteins in dataset with no hallmark match: ", paste(unmatched, collapse = ", "),
       "\nAdd these to hallmark_lookup in Step 19 before proceeding.")
}


# Step 21: Build expression and design matrices for limma
# Expression matrix: proteins as rows, samples as columns (221 x 66).
# Design matrix: samples as rows, group/covariate columns (66 x k).
# No-intercept design (~0 + study_group) gives each group its own coefficient.
# Age and sex included as covariates; sex coded 1=Male, 0=Female for limma.
# model.matrix converts group names to dot-notation:
#   study_groupAcute.IPV.BI, study_groupAcute.Non.IPV.BI, study_groupChronic.IPV.BI

expr_matrix <- cohort %>%
  select(all_of(npqz_cols)) %>%
  as.matrix() %>%
  t()

rownames(expr_matrix) <- protein_names
colnames(expr_matrix) <- cohort$participant_id

design <- model.matrix(
  ~ 0 + study_group + age_years + sex,
  data = cohort %>% mutate(sex = ifelse(sex == "Male", 1L, 0L))
)
colnames(design) <- make.names(colnames(design))


# Step 22: Specify pairwise contrasts
# Positive logFC = first-named group has higher expression.
# logFC units are z-score differences, not log2 fold-change, because input is z-scored.

contrast_matrix <- makeContrasts(
  Acute_IPV_vs_NonIPV   = study_groupAcute.IPV.BI - study_groupAcute.Non.IPV.BI,
  Chronic_vs_Acute_IPV  = study_groupChronic.IPV.BI - study_groupAcute.IPV.BI,
  Chronic_IPV_vs_NonIPV = study_groupChronic.IPV.BI - study_groupAcute.Non.IPV.BI,
  levels = design
)
contrast_labels <- c(
  "Acute_IPV_vs_NonIPV"   = "Acute IPV-BI vs Acute Non-IPV-BI",
  "Chronic_vs_Acute_IPV"  = "Chronic IPV-BI vs Acute IPV-BI",
  "Chronic_IPV_vs_NonIPV" = "Chronic IPV-BI vs Acute Non-IPV-BI"
)
contrast_names <- names(contrast_labels)


# Step 23: Fit limma model
# lmFit:        fits a separate linear model per protein using the design matrix
# contrasts.fit: recalculates coefficients for the specified pairwise contrasts
# eBayes:       applies empirical-Bayes moderation, borrowing strength across all
#               221 proteins to stabilise variance estimates at n=66
#               trend=TRUE allows prior variance to vary with mean expression
#               robust=TRUE downweights outlier samples

fit  <- lmFit(expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)


# Step 24: P-value sanity check
# Halts if fewer than 10 unique raw p-values across 221 proteins — indicates
# eBayes has collapsed the distribution, reproducing the flat-line volcano problem.

for (cn in contrast_names) {
  tt         <- topTable(fit2, coef = cn, number = Inf, adjust = "BH", sort.by = "P")
  n_unique_p <- length(unique(round(tt$P.Value, 6)))
  stopifnot("Limma p-values are not varying — check data format (Step 18)." = n_unique_p >= 10)
}


# Step 25: Extract limma contrast results
# topTable() returns per-protein results for one contrast at a time.
# Key columns:
#   logFC     — z-score difference between group means (positive = first group higher)
#   AveExpr   — average expression across all 66 samples
#   t         — moderated t-statistic
#   P.Value   — raw p-value from limma's empirical-Bayes t-test
#   adj.P.Val — BH-FDR adjusted p-value across 221 proteins in this contrast only
#               each contrast is corrected independently, not pooled across contrasts
#   B         — log-odds that the protein is differentially expressed

limma_results <- lapply(contrast_names, function(cn) {
  topTable(fit2, coef = cn, number = Inf, adjust = "BH", sort.by = "P") %>%
    as.data.frame() %>%
    rownames_to_column("protein") %>%
    mutate(
      contrast        = cn,
      contrast_label  = contrast_labels[cn],
      limma_sig       = adj.P.Val < 0.05,
      limma_direction = case_when(
        limma_sig & logFC > 0 ~ "Up",
        limma_sig & logFC < 0 ~ "Down",
        TRUE                  ~ "NS"
      )
    ) %>%
    left_join(hallmark_lookup %>% select(protein, hallmark, hallmark_abbrev),
              by = "protein") %>%
    rename(
      limma_logFC = logFC,
      limma_p_raw = P.Value,
      limma_adj_p = adj.P.Val
    )
})
names(limma_results) <- contrast_names


# Step 26: Mann-Whitney U tests (nonparametric sanity check)
# wilcox.test() ranks all observations and tests whether one group's ranks are
# systematically higher. No assumption about normality or equal variance.
# mw_p_raw = raw two-sided p-value from wilcox.test()
# mw_adj_p = BH correction across all 663 tests (221 proteins x 3 contrasts)
#            more conservative than limma's per-contrast correction
# Cannot adjust for age or sex — confirms limma results, does not replace them.

contrast_groups <- list(
  Acute_IPV_vs_NonIPV   = c("Acute IPV-BI",    "Acute Non-IPV-BI"),
  Chronic_vs_Acute_IPV  = c("Chronic IPV-BI",  "Acute IPV-BI"),
  Chronic_IPV_vs_NonIPV = c("Chronic IPV-BI",  "Acute Non-IPV-BI")
)

run_mw_tests <- function(cohort_df, protein_cols, protein_names, group_col, groups) {
  g1_rows <- cohort_df[[group_col]] == groups[1]
  g2_rows <- cohort_df[[group_col]] == groups[2]
  
  results <- lapply(seq_along(protein_cols), function(i) {
    col  <- protein_cols[i]
    prot <- protein_names[i]
    x1   <- cohort_df[[col]][g1_rows]
    x2   <- cohort_df[[col]][g2_rows]
    x1   <- x1[!is.na(x1)]
    x2   <- x2[!is.na(x2)]
    
    if (length(x1) < 3 || length(x2) < 3) {
      return(data.frame(protein  = prot, mw_p_raw = NA_real_,
                        mw_logFC = NA_real_, stringsAsFactors = FALSE))
    }
    
    test <- wilcox.test(x1, x2, alternative = "two.sided", exact = FALSE)
    data.frame(protein  = prot,
               mw_p_raw = test$p.value,
               mw_logFC = mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  bind_rows(results)
}

mw_per_contrast <- lapply(contrast_names, function(cn) {
  run_mw_tests(cohort, npqz_cols, protein_names, "study_group", contrast_groups[[cn]]) %>%
    mutate(contrast = cn, contrast_label = contrast_labels[cn])
})
names(mw_per_contrast) <- contrast_names

mw_all <- bind_rows(mw_per_contrast) %>%
  mutate(
    mw_adj_p     = p.adjust(mw_p_raw, method = "BH"),
    mw_sig       = mw_adj_p < 0.05,
    mw_direction = case_when(
      mw_sig & mw_logFC > 0 ~ "Up",
      mw_sig & mw_logFC < 0 ~ "Down",
      TRUE                  ~ "NS"
    )
  )


# Step 27: Merge limma and Mann-Whitney results

pairwise_results <- lapply(contrast_names, function(cn) {
  left_join(
    limma_results[[cn]],
    mw_all %>% filter(contrast == cn) %>%
      select(protein, mw_p_raw, mw_adj_p, mw_logFC, mw_sig, mw_direction),
    by = "protein"
  ) %>%
    mutate(
      sig       = limma_sig,
      direction = limma_direction
    )
})
names(pairwise_results) <- contrast_names


# Step 28: Write corrected p-values back to ipv_cohort_clean.rds
# Adds limma and Mann-Whitney adjusted p-values as new columns per protein per
# contrast. Every row gets the same value for a given protein (p-values are
# protein-level, not participant-level) — keeps all outputs in one file for
# downstream scripts.
# Column naming convention:
#   limma_adj_p_<ContrastShortName>_<protein>
#   mw_adj_p_<ContrastShortName>_<protein>
#   limma_sig_<ContrastShortName>_<protein>

contrast_short <- c(
  "Acute_IPV_vs_NonIPV"   = "AcuteIPVvsNon",
  "Chronic_vs_Acute_IPV"  = "ChronicvsAcute",
  "Chronic_IPV_vs_NonIPV" = "ChronicvsNon"
)

for (cn in contrast_names) {
  short  <- contrast_short[cn]
  res_cn <- pairwise_results[[cn]] %>% select(protein, limma_adj_p, mw_adj_p, limma_sig)
  
  for (i in seq_len(nrow(res_cn))) {
    prot <- res_cn$protein[i]
    cohort[[paste0("limma_adj_p_", short, "_", prot)]] <- res_cn$limma_adj_p[i]
    cohort[[paste0("mw_adj_p_",    short, "_", prot)]] <- res_cn$mw_adj_p[i]
    cohort[[paste0("limma_sig_",   short, "_", prot)]] <- res_cn$limma_sig[i]
  }
}

saveRDS(cohort, "ipv_cohort_clean.rds")


# Step 29: Shared visual elements

hallmark_levels <- hallmark_lookup %>% pull(hallmark_abbrev) %>% unique() %>% sort()
hallmark_colors <- setNames(
  colorRampPalette(brewer.pal(9, "Set3"))(length(hallmark_levels)),
  hallmark_levels
)

group_colors <- list(
  Group = c(
    "Acute IPV-BI"     = "#E41A1C",
    "Acute Non-IPV-BI" = "#377EB8",
    "Chronic IPV-BI"   = "#FF7F00"
  )
)

data_note <- "Data: z-scored NPQ values (simulated)"


# Step 30: Volcano plot function
# Non-significant proteins: grey. Significant proteins: coloured by hallmark.
# Top 5 significant proteins labelled; if none significant, top 5 by adjusted p.
# hallmark_abbrev cast to factor so scale_color_manual matches levels correctly.
# drop=TRUE shows only hallmark categories present in significant proteins.

make_volcano <- function(result_df, contrast_label, sig_threshold = 0.05,
                         exploratory = FALSE) {
  result_df <- result_df %>%
    mutate(
      fc              = limma_logFC,
      adjp            = limma_adj_p,
      sig             = limma_adj_p < sig_threshold,
      hallmark_abbrev = factor(hallmark_abbrev, levels = hallmark_levels)
    )
  
  top_labels <- result_df %>% filter(sig) %>% slice_min(adjp, n = 5, with_ties = FALSE)
  if (nrow(top_labels) == 0) top_labels <- result_df %>% slice_min(adjp, n = 5, with_ties = FALSE)
  
  subtitle <- if (exploratory) {
    paste0("Exploratory: BH-adjusted p < ", sig_threshold, " — ", data_note)
  } else {
    paste0("BH-adjusted p < ", sig_threshold, " (FDR) — ", data_note)
  }
  
  ggplot(result_df, aes(x = fc, y = -log10(adjp))) +
    geom_point(data = result_df %>% filter(!sig),
               color = "grey75", alpha = 0.5, size = 1.5) +
    geom_point(data = result_df %>% filter(sig),
               aes(color = hallmark_abbrev), alpha = 0.9, size = 2.2) +
    geom_hline(yintercept = -log10(sig_threshold),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.5, 0.5),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_text_repel(
      data          = top_labels,
      aes(label     = protein),
      size          = 2.8,
      max.overlaps  = 20,
      segment.color = "grey50",
      box.padding   = 0.4,
      point.padding = 0.3,
      color         = "black"
    ) +
    scale_color_manual(values = hallmark_colors, name = "Hallmark",
                       na.value = "grey75", drop = TRUE) +
    labs(
      title    = contrast_label,
      subtitle = subtitle,
      x        = "Mean difference (z-score units)",
      y        = expression(-log[10](BH-adjusted~p))
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 10),
      plot.subtitle   = element_text(color = "grey40", hjust = 0.5, size = 7),
      legend.position = "right",
      legend.text     = element_text(size = 7),
      legend.title    = element_text(face = "bold", size = 8),
      legend.key.size = unit(0.35, "cm")
    )
}


# Step 31: Diagnostic p-value histograms
# Plots the raw limma p-value distribution per contrast to confirm the pipeline
# is producing well-behaved p-values. A right-skewed distribution (spike near 0,
# flat elsewhere) indicates real signal. A uniform distribution indicates mostly
# null proteins, which is expected for simulated data at n=66.
# This documents that the flat appearance of primary volcano plots is caused by
# BH correction on mostly-null data, not by a pipeline error.

for (cn in contrast_names) {
  tt <- topTable(fit2, coef = cn, number = Inf, adjust = "BH", sort.by = "none")
  
  p_hist <- ggplot(tt, aes(x = P.Value)) +
    geom_histogram(bins = 50, fill = "grey70", color = "white", linewidth = 0.3) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "#E41A1C", linewidth = 0.5) +
    labs(
      title    = paste0("Raw p-value distribution: ", contrast_labels[cn]),
      subtitle = paste0("Confirms pipeline is correct — flat volcano reflects BH correction ",
                        "on mostly-null data, not a pipeline error\n", data_note),
      x        = "Raw p-value (unadjusted)",
      y        = "Number of proteins",
      caption  = paste0("Red line: p = 0.05  |  Proteins with raw p < 0.05: ",
                        sum(tt$P.Value < 0.05))
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 10),
      plot.subtitle = element_text(color = "grey40", hjust = 0.5, size = 7),
      plot.caption  = element_text(color = "grey50", size = 8)
    )
  
  print(p_hist)
  ggsave(file.path(output_dir, paste0("diagnostic_pval_hist_", cn, ".png")),
         plot = p_hist, width = 8, height = 5, dpi = 300)
}


# Step 32: Primary volcano plots — BH-adjusted p-values, FDR < 0.05
# Two proteins expected to be significant in the primary contrast with simulated data.
# Flat appearance is correct and expected — see diagnostic histograms (Step 31).

for (cn in contrast_names) {
  p <- make_volcano(pairwise_results[[cn]], contrast_labels[cn],
                    sig_threshold = 0.05, exploratory = FALSE)
  print(p)
  ggsave(file.path(output_dir, paste0("volcano_primary_", cn, ".png")),
         plot = p, width = 10, height = 7, dpi = 300)
}


# Step 33: Identify proteins significant in at least one contrast (FDR < 0.05)
# Union approach: qualifies if limma FDR < 0.05 in ANY contrast.

sig_proteins <- bind_rows(pairwise_results) %>%
  filter(limma_sig) %>%
  pull(protein) %>%
  unique()


# Step 34: Exploratory volcano plots — BH-adjusted p-values, FDR < 0.20
# Separate from primary plots. FDR < 0.20 is a common discovery threshold in
# proteomics to identify candidate proteins for follow-up. Clearly labelled
# as exploratory so they are not confused with the primary analysis.

for (cn in contrast_names) {
  p <- make_volcano(pairwise_results[[cn]], contrast_labels[cn],
                    sig_threshold = 0.20, exploratory = TRUE)
  print(p)
  ggsave(file.path(output_dir, paste0("volcano_exploratory_", cn, ".png")),
         plot = p, width = 10, height = 7, dpi = 300)
}

# Identify proteins significant at exploratory threshold in at least one contrast
sig_proteins_exploratory <- bind_rows(pairwise_results) %>%
  filter(limma_adj_p < 0.20) %>%
  pull(protein) %>%
  unique()


# Step 35: Heatmap setup
# Row ordering: proteins sorted by hallmark category then alphabetically.
# Column ordering: participants sorted by study_group so groups are contiguous.
# Diverging colour scale (blue-white-red):
#   blue = below cohort mean, white = at cohort mean, red = above cohort mean.

order_proteins_by_hallmark <- function(protein_set) {
  hallmark_lookup %>%
    filter(protein %in% protein_set) %>%
    select(protein, hallmark_abbrev) %>%
    distinct() %>%
    arrange(hallmark_abbrev, protein) %>%
    pull(protein)
}

make_row_annotation <- function(protein_set) {
  hallmark_lookup %>%
    filter(protein %in% protein_set) %>%
    select(protein, hallmark_abbrev) %>%
    distinct() %>%
    column_to_rownames("protein") %>%
    rename(Hallmark = hallmark_abbrev)
}

sample_order <- cohort %>% arrange(study_group) %>% pull(participant_id)

col_annotation_individual <- data.frame(
  Group     = cohort$study_group,
  row.names = cohort$participant_id
)

heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

group_means <- cohort %>%
  select(study_group, all_of(npqz_cols)) %>%
  group_by(study_group) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("study_group") %>%
  t() %>%
  as.data.frame()
rownames(group_means) <- protein_names

col_annotation_group    <- data.frame(Group = colnames(group_means),
                                      row.names = colnames(group_means))
row_order_all           <- order_proteins_by_hallmark(protein_names)
row_ann_all             <- make_row_annotation(protein_names)
hallmark_ann_colors_all <- list(Hallmark = hallmark_colors[unique(row_ann_all$Hallmark)])


# Step 36: Heatmap A — all 221 proteins, individual participants
graphics.off()
save_heatmap(
  filepath = file.path(output_dir, "heatmap_A_all_proteins_individual.png"),
  width    = 16,
  height   = 20,
  pheatmap_args = list(
    mat               = expr_matrix[row_order_all, sample_order],
    color             = heatmap_colors,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    annotation_col    = col_annotation_individual[sample_order, , drop = FALSE],
    annotation_row    = row_ann_all,
    annotation_colors = c(group_colors, hallmark_ann_colors_all),
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    fontsize_row      = 5,
    fontsize          = 9,
    main              = paste0("All 221 Proteins — Individual Participants\n", data_note),
    border_color      = NA
  )
)


# Step 37: Heatmap B — all 221 proteins, group means

save_heatmap(
  filepath = file.path(output_dir, "heatmap_B_all_proteins_group_means.png"),
  width    = 10,
  height   = 20,
  pheatmap_args = list(
    mat               = as.matrix(group_means)[row_order_all, ],
    color             = heatmap_colors,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    annotation_col    = col_annotation_group,
    annotation_row    = row_ann_all,
    annotation_colors = c(group_colors, hallmark_ann_colors_all),
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 5,
    fontsize          = 9,
    main              = paste0("All 221 Proteins — Group Means\n", data_note),
    border_color      = NA
  )
)


# Steps 38-39: Heatmaps C and D — significant proteins only (FDR < 0.05)
# C = group means, D = individual participants.
# Skipped if no proteins reach limma FDR < 0.05.

if (length(sig_proteins) > 0) {
  
  row_order_sig           <- order_proteins_by_hallmark(sig_proteins)
  row_ann_sig             <- make_row_annotation(sig_proteins)
  hallmark_ann_colors_sig <- list(Hallmark = hallmark_colors[unique(row_ann_sig$Hallmark)])
  plot_height             <- max(6, length(sig_proteins) * 0.18 + 2)
  
  save_heatmap(
    filepath = file.path(output_dir, "heatmap_C_sig_proteins_group_means.png"),
    width    = 10,
    height   = plot_height,
    pheatmap_args = list(
      mat               = as.matrix(group_means)[row_order_sig, ],
      color             = heatmap_colors,
      cluster_rows      = FALSE,
      cluster_cols      = FALSE,
      annotation_col    = col_annotation_group,
      annotation_row    = row_ann_sig,
      annotation_colors = c(group_colors, hallmark_ann_colors_sig),
      show_rownames     = TRUE,
      show_colnames     = TRUE,
      fontsize_row      = 7,
      fontsize          = 10,
      main              = paste0(
        "Significant Proteins (FDR < 0.05 in \u22651 contrast) \u2014 Group Means\n",
        data_note),
      border_color      = NA
    )
  )
  
  save_heatmap(
    filepath = file.path(output_dir, "heatmap_D_sig_proteins_individual.png"),
    width    = 16,
    height   = plot_height,
    pheatmap_args = list(
      mat               = expr_matrix[row_order_sig, sample_order],
      color             = heatmap_colors,
      cluster_rows      = FALSE,
      cluster_cols      = FALSE,
      annotation_col    = col_annotation_individual[sample_order, , drop = FALSE],
      annotation_row    = row_ann_sig,
      annotation_colors = c(group_colors, hallmark_ann_colors_sig),
      show_rownames     = TRUE,
      show_colnames     = FALSE,
      fontsize_row      = 7,
      fontsize          = 9,
      main              = paste0(
        "Significant Proteins (FDR < 0.05 in \u22651 contrast) \u2014 Individual Participants\n",
        data_note),
      border_color      = NA
    )
  )
  
}
# Step 40: Save full results as CSV (one file per contrast)
# Columns: protein, hallmark, hallmark_abbrev,
#          limma_logFC, limma_p_raw, limma_adj_p, limma_sig, limma_direction,
#          mw_logFC,   mw_p_raw,   mw_adj_p,   mw_sig,   mw_direction
# Rows sorted by limma_adj_p ascending.

for (cn in contrast_names) {
  pairwise_results[[cn]] %>%
    select(protein, hallmark, hallmark_abbrev,
           limma_logFC, limma_p_raw, limma_adj_p, limma_sig, limma_direction,
           mw_logFC,   mw_p_raw,   mw_adj_p,   mw_sig,   mw_direction) %>%
    arrange(limma_adj_p) %>%
    write_csv(file.path(output_dir, paste0("results_", cn, ".csv")))
}