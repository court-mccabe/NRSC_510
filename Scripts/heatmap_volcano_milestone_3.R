# Weeks 5-6: NULISA Proteomics — Limma Models, Volcano Plots, and Heatmaps

# Goal: Run differential expression models across the three study groups using
# limma, visualize results with volcano plots per contrast, and produce heatmaps
# organized by hallmark category of neurodegeneration.

# Pairwise contrasts by scientific priority:
#   Primary:   Acute IPV-BI vs Acute Non-IPV-BI
#              (effect of IPV-BI within the acute phase)
#   Secondary: Chronic IPV-BI vs Acute IPV-BI
#              (effect of time/chronicity within IPV-BI)
#   Tertiary:  Chronic IPV-BI vs Acute Non-IPV-BI
#              (cumulative difference: chronicity + IPV-BI vs acute non-injury)

# DATA NOTE — SIMULATED DATASET:
# The NPQz_ columns contain z-scored NPQ (NULISA Protein Quantification) values.
# In real Alamar NULISA data, the pipeline begins with raw NPQ values requiring
# log2 transformation before z-scoring. In this simulated dataset, values were
# generated as already log2-transformed and z-scored. All plots therefore carry
# the subtitle "Data: log2-transformed, z-scored NPQ values (simulated)" to
# accurately describe the input scale. When applying this pipeline to real
# NULISA data, insert the log2 transformation step noted in Step 17 before
# running the limma model.

# EXCLUSION NOTE:
# NPQz_MSLN was present in the simulated dataset but was not assigned a hallmark
# category in the Alamar Neuro 220 Panel documentation. This protein was
# incorrectly included during data simulation and is excluded from all downstream
# analyses. All results therefore reflect 221 proteins.

# DETECTION LIMIT NOTE:
# In real Alamar NULISA data, proteins must pass a 50% detection threshold per
# group before inclusion in differential expression analysis. Detection is
# determined by comparing each protein's NPQ signal to inter-plate controls
# (IPC) and intra-plate controls (IC) provided by Alamar's output pipeline.
# A protein is considered reliably detected if its signal exceeds background
# in >=50% of samples within each group. This simulated dataset does not include
# IPC or IC columns, so this QC step is scaffolded in Step 18 but cannot be
# executed on simulated data. Apply Step 18 fully before running real NULISA
# data through this pipeline.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("limma", update = FALSE, ask = FALSE)

# ggrepel and pheatmap are on CRAN but included here for completeness
# in case they are also missing
install.packages(
  c("ggrepel", "pheatmap", "RColorBrewer"),
  dependencies = TRUE
)

library(tidyverse)
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)

# Step 16: Load cleaned dataset and set working directory

setwd("C:/NRSC_510")
cohort <- readRDS("ipv_cohort_clean.rds")

output_dir <- "C:/NRSC_510/Outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# REAL DATA EXTENSION NOTE:
# If applying this pipeline to raw NPQ values (not yet log2-transformed),
# insert the following block here before proceeding to Step 18:
#
# npqz_cols_raw <- cohort %>% select(starts_with("NPQz_")) %>% names()
# cohort <- cohort %>%
#   mutate(across(
#     all_of(npqz_cols_raw),
#     ~ log2(.x + 1)
#   ))
#
# The +1 offset handles any zero values before log2 transformation.
# Note that for limma, log2-transformed values without z-scoring are preferred
# because z-scoring removes between-sample variance that limma relies on for
# modelling. Check Alamar's pipeline documentation for the recommended input
# format for your specific data version.

# Step 18: Detection limit QC scaffold (simulated data — cannot execute)
#
# In real Alamar NULISA data, replace this section with the following:
#
# detection_flags <- read_csv("alamar_detection_flags.csv")
#   This file contains one row per sample, one column per protein,
#   with binary values: 1 = detected above background, 0 = not detected.
#   It is provided by Alamar's standard output pipeline alongside NPQ values.
#
# detection_rate <- detection_flags %>%
#   left_join(cohort %>% select(participant_id, study_group),
#             by = "participant_id") %>%
#   group_by(study_group) %>%
#   summarise(across(starts_with("NPQz_"), ~ mean(.x, na.rm = TRUE)),
#             .groups = "drop")
#
# proteins_pass <- detection_rate %>%
#   select(-study_group) %>%
#   summarise(across(everything(), ~ all(.x >= 0.50))) %>%
#   pivot_longer(everything(),
#                names_to  = "protein",
#                values_to = "pass") %>%
#   filter(pass) %>%
#   pull(protein)
#
# cat("Proteins passing 50% detection in all groups:", length(proteins_pass), "\n")
# cat("Proteins excluded by detection filter:", 222 - length(proteins_pass), "\n")
# cohort <- cohort %>% select(-any_of(setdiff(
#   names(cohort)[startsWith(names(cohort), "NPQz_")], proteins_pass
# )))

cat("DETECTION NOTE: 50% detection QC skipped — simulated data lacks IPC/IC controls.\n")
cat("All 221 proteins retained (after MSLN exclusion). See header note.\n\n")

# Step 19: Build hallmark lookup table from Alamar Neuro 220 Panel documentation

# The panel CSV uses parenthetical aliases in protein names (e.g., "CCL2 (MCP-1)")
# that do not match the NPQz_ column names in the dataset. Names are cleaned by:
#   1. Stripping everything in parentheses
#   2. Replacing semicolons (for dual-gene entries like HBA1; HBA2) with underscores
#   3. Replacing remaining spaces with hyphens
# Two additional manual fixes:
#   - "IL-12A|IL-12B" cleaned target maps to dataset column NPQz_IL-12p70
#   - "Aberrant proteostasis" (lowercase in some CSV rows) is harmonised to
#     title case "Aberrant Proteostasis" for consistent grouping
# MSLN is excluded — no hallmark assignment, incorrectly included in simulation.

panel_raw <- read_csv("Neuro_220_Panel.csv", show_col_types = FALSE)

hallmark_lookup <- panel_raw %>%
  rename(
    target_full = `Target List`,
    hallmark    = `Hallmarks of Neurodegenerative Disease`,
    uniprot_id  = `Uniprot ID`
  ) %>%
  mutate(
    protein = target_full %>%
      str_replace_all("\\s*\\(.*?\\)", "") %>% # remove everything in parentheses and space before
      str_replace_all(";\\s*", "_") %>% # replace semicolon with underscore
      str_replace_all("\\|.*", "") %>% # remove everything after pipe character
      str_replace_all("\\s+", "-") %>% # replace spaces with - 
      str_trim(),
    protein = case_when(
      str_detect(target_full, "IL-12A") ~ "IL-12p70",
      TRUE                              ~ protein
    ),
    hallmark = case_when(
      hallmark == "Aberrant proteostasis" ~ "Aberrant Proteostasis",
      TRUE                               ~ hallmark
    ),
    hallmark_abbrev = case_when( # make names more code friendly 
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
  filter(!is.na(hallmark)) %>% # keep rows that is not NA
  select(protein, uniprot_id, hallmark, hallmark_abbrev, target_full)

print(table(hallmark_lookup$hallmark_abbrev))

# Step 20: Identify and validate NPQz protein columns

#subset NPQz columns and extract as character vector
npqz_cols_all  <- cohort %>% select(starts_with("NPQz_")) %>% names()
#transform NPQz columns by removing "NPQz" from start of column names
protein_names_all <- str_remove(npqz_cols_all, "^NPQz_")

excluded_proteins <- c("MSLN")

npqz_cols     <- npqz_cols_all[!protein_names_all %in% excluded_proteins]
protein_names <- str_remove(npqz_cols, "^NPQz_")

cat("Total NPQz columns in dataset:  ", length(npqz_cols_all), "\n")
cat("Excluded (no hallmark):          ", paste(excluded_proteins, collapse = ", "), "\n")
cat("Proteins retained for analysis:  ", length(npqz_cols), "\n\n")

#logical index: create vector of any protein names that exist in dataset
#but have no entry in hallmark_loop
unmatched <- protein_names[!protein_names %in% hallmark_lookup$protein]
if (length(unmatched) > 0) {
  warning("Proteins in dataset with no hallmark match: ",
          paste(unmatched, collapse = ", "))
} else {
  cat("All", length(protein_names),
      "retained proteins matched to hallmark lookup.\n\n")
}

# Step 21: Build expression matrix and design matrix for limma

# limma expects:
#   - expression matrix: proteins as rows, samples as columns (p x n)
#   - design matrix:     samples as rows, group/covariate columns (n x k)

# A no-intercept design (~0 + study_group) gives each group its own
# coefficient, making pairwise contrasts straightforward to specify.
# Age and sex are included as covariates, consistent with milestone 2 models.
# sex is converted to 0/1 numeric for limma compatibility.

# R's model.matrix converts factor level names to column-safe strings by
# replacing spaces and hyphens with dots, producing:
#   study_groupAcute.IPV.BI
#   study_groupAcute.Non.IPV.BI
#   study_groupChronic.IPV.BI
# These exact strings are used in the contrast specification in Step 22.

expr_matrix <- cohort %>%
  select(all_of(npqz_cols)) %>%
  as.matrix() %>% # create transposed matrix (limma expects proteins as rows)
  t()

rownames(expr_matrix) <- protein_names
colnames(expr_matrix) <- cohort$participant_id

cat("Expression matrix: ", nrow(expr_matrix), "proteins x",
    ncol(expr_matrix), "samples\n\n") # confirm 221 proteins and 66 samples

design <- model.matrix( # function telling limma how samples relate to each other
  ~ 0 + study_group + age_years + sex, # each row is one participant, column is group indicator or covariate
  data = cohort %>%
    mutate(sex = ifelse(sex == "Male", 1L, 0L))
)

colnames(design) <- make.names(colnames(design))

cat("Design matrix columns:\n")
print(colnames(design))
cat("\n")

# Step 22: Specify pairwise contrasts

# Column names from model.matrix use dot notation for special characters.
# The three contrasts below are listed in order of scientific priority.
# Positive logFC means the first-named group has higher expression.
print(colnames(design))

contrast_matrix <- makeContrasts(
  Acute_IPV_vs_NonIPV   = study_groupAcute.IPV.BI - study_groupAcute.Non.IPV.BI,
  Chronic_vs_Acute_IPV  = study_groupChronic.IPV.BI - study_groupAcute.IPV.BI,
  Chronic_IPV_vs_NonIPV = study_groupChronic.IPV.BI - study_groupAcute.Non.IPV.BI,
  levels = design # read column names from matrix built
)

contrast_labels <- c(
  "Acute_IPV_vs_NonIPV"   = "Acute IPV-BI vs Acute Non-IPV-BI",
  "Chronic_vs_Acute_IPV"  = "Chronic IPV-BI vs Acute IPV-BI",
  "Chronic_IPV_vs_NonIPV" = "Chronic IPV-BI vs Acute Non-IPV-BI"
)

contrast_names <- names(contrast_labels)

# Step 23: Fit limma model

# lmFit fits a linear model per protein using the design matrix.
# contrasts.fit recalculates coefficients for the specified contrasts.
# eBayes applies empirical Bayes moderation of variance estimates,
# borrowing strength across all 221 proteins to stabilise estimates
# at small n. trend = TRUE allows the prior variance to trend with
# mean expression; robust = TRUE down-weights outlier samples.

fit  <- lmFit(expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE) # must run after contrasts

cat("Limma model fitted successfully.\n\n") # confirm

# Step 25: Extract pairwise contrast results

# topTable is called once per contrast. Key output columns:
#   logFC     — mean difference in z-score units between the two groups
#               (positive = first-named group higher)
#   AveExpr   — average expression across all 66 samples
#   t         — moderated t-statistic
#   P.Value   — raw p-value
#   adj.P.Val — BH-FDR adjusted p-value across the 221 proteins
#   B         — log-odds that the protein is differentially expressed

# NOTE on logFC units: because input values are z-scored, logFC represents
# a difference in z-score units rather than a true fold-change in raw
# expression. A logFC of 1.0 means the first group is on average 1 SD
# higher than the second on that protein's normalized scale.

pairwise_results <- lapply(contrast_names, function(cn) { # return a list
  topTable(
    fit2,
    coef    = cn,
    number  = Inf,
    adjust  = "BH",
    sort.by = "P"
  ) %>%
    as.data.frame() %>%
    rownames_to_column("protein") %>%
    mutate(
      contrast       = cn,
      contrast_label = contrast_labels[cn],
      sig            = adj.P.Val < 0.05,
      direction      = case_when(
        sig & logFC > 0 ~ "Up",
        sig & logFC < 0 ~ "Down",
        TRUE            ~ "NS"
      )
    ) %>%
    left_join(
      hallmark_lookup %>% select(protein, hallmark, hallmark_abbrev),
      by = "protein"
    )
})

names(pairwise_results) <- contrast_names

cat("Pairwise contrast results summary:\n")
for (cn in contrast_names) {
  res <- pairwise_results[[cn]]
  cat(sprintf("  %-45s  sig: %d  (up: %d  down: %d)\n",
              contrast_labels[cn],
              sum(res$sig, na.rm = TRUE),
              sum(res$direction == "Up", na.rm = TRUE),
              sum(res$direction == "Down", na.rm = TRUE)))
}
cat("\n")

# Step 26: Hallmark category enrichment
# DATA LIMITATION: 
# Since only two proteins show as signficant, this block of code is kind of 
# useless. It is used to determine which hallmark would be more prominent 
# in a data set. For example, if 30 markers were significant and 20 were 
# inflammatory, we could test that here. Not using for this simualted set.

# For each contrast, a hypergeometric test asks whether each hallmark
# category contains more significant proteins than expected by chance,
# given the total number of proteins tested and the total number significant.
# This is the standard over-representation analysis (ORA) approach.
# BH-FDR applied across hallmark categories within each contrast separately.

#run_hallmark_enrichment <- function(result_df, contrast_label) {
#  total_n   <- nrow(result_df)
#  total_sig <- sum(result_df$sig, na.rm = TRUE)
  
#  result_df %>%
#    group_by(hallmark_abbrev) %>%
#    summarise(
#      n_category = n(),
#      n_sig      = sum(sig, na.rm = TRUE),
#      pct_sig    = round(n_sig / n_category * 100, 1),
#      .groups    = "drop"
#    ) %>%
#    mutate(
#      p_hypergeom = phyper(
#        q          = n_sig - 1,
#        m          = total_sig,
#        n          = total_n - total_sig,
#        k          = n_category,
#        lower.tail = FALSE
#      ),
#      p_adj_BH = p.adjust(p_hypergeom, method = "BH"),
#      contrast = contrast_label
#    ) %>%
#    arrange(p_hypergeom)
#}

#hallmark_enrichment <- bind_rows(
#  lapply(contrast_names, function(cn) {
#    run_hallmark_enrichment(pairwise_results[[cn]], contrast_labels[cn])
#  })
#)

#sig_enrichment <- hallmark_enrichment %>% filter(p_adj_BH < 0.05)
#if (nrow(sig_enrichment) > 0) {
#  print(sig_enrichment %>%
#          select(contrast, hallmark_abbrev, n_category, n_sig, pct_sig, p_adj_BH))
#} else {
#  cat("  None at FDR < 0.05 (expected with simulated data at n=66)\n")
#}
#cat("\n")

# Step 27: Define shared visual elements
#
# Hallmark color palette: colorblind-friendly via RColorBrewer Set3,
# extended to cover all 9 hallmark categories present in the panel.
# The same palette is used across all volcano plots and heatmap annotations
# so category colors are consistent throughout the figure set.

hallmark_levels <- hallmark_lookup %>%
  pull(hallmark_abbrev) %>%
  unique() %>%
  sort()

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

data_note <- "Data: log2-transformed, z-scored NPQ values (simulated)"

# Step 28: Volcano plots — one per contrast

# X-axis: logFC (z-score units — difference in group mean z-scores)
# Y-axis: -log10(BH-FDR adjusted p-value)
# Non-significant proteins: grey
# Significant proteins: colored by hallmark category
# Top 10 most significant proteins labelled with ggrepel
# Reference lines: FDR = 0.05 (horizontal), |logFC| = 0.5 (vertical)
# A logFC threshold of 0.5 SD is used as a practical effect size guide;
# it is not a hard filter — all proteins at FDR < 0.05 are considered
# significant regardless of logFC magnitude.

make_volcano <- function(result_df, contrast_label) { #result_df is function argument
  
  top_labels <- result_df %>% 
    filter(sig) %>%
    slice_min(adj.P.Val, n = 10, with_ties = FALSE) # sort on p value, keep 10 sig values
  
  ggplot(result_df,
         aes(x = logFC, y = -log10(adj.P.Val), color = hallmark_abbrev)) +
    geom_point(
      data  = result_df %>% filter(!sig), # not significant values become grey
      color = "grey75", alpha = 0.5, size = 1.8
    ) +
    geom_point(
      data  = result_df %>% filter(sig), # sig values follow hallmark colors
      aes(color = hallmark_abbrev), alpha = 0.85, size = 2.2 # alpha fxn is for less transparency
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype   = "dashed", color = "grey40", linewidth = 0.5
    ) +
    geom_vline(
      xintercept = c(-0.5, 0.5),
      linetype   = "dashed", color = "grey40", linewidth = 0.5
    ) +
    geom_text_repel(
      data          = top_labels,
      aes(label     = protein),
      size          = 3,
      max.overlaps  = 15,
      segment.color = "grey50",
      box.padding   = 0.4,
      point.padding = 0.3,
      color         = "black"
    ) +
    scale_color_manual(values = hallmark_colors,
                       name   = "Hallmark Category") +
    labs(
      title    = paste0("Volcano Plot: ", contrast_label),
      subtitle = data_note,
      x        = "log2 Fold Change (z-score units)",
      y        = expression(-log[10](FDR~adjusted~p-value)),
      caption  = paste0(
        "Dashed lines: FDR = 0.05 (horizontal), |logFC| = 0.5 (vertical)\n",
        "Top 10 most significant proteins labelled"
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(color = "grey40", hjust = 0.5, size = 9),
      plot.caption  = element_text(color = "grey50", size = 8),
      legend.position  = "right",
      legend.text      = element_text(size = 8),
      legend.title     = element_text(face = "bold", size = 9),
      legend.key.size  = unit(0.4, "cm")
    )
}

for (cn in contrast_names) {
  p        <- make_volcano(pairwise_results[[cn]], contrast_labels[cn])
  filename <- file.path(output_dir, paste0("volcano_", cn, ".png"))
  ggsave(filename, plot = p, width = 10, height = 7, dpi = 300)
  cat("Saved:", filename, "\n")
}

# Step 29: Identify proteins significant in at least one contrast

# A protein is included in the significant-protein heatmaps (C and D) if
# it reaches FDR < 0.05 in at least one of the three pairwise contrasts.
# This is a union across contrasts — significance in any one comparison
# qualifies a protein for inclusion.

sig_proteins <- bind_rows(pairwise_results) %>%
  filter(sig) %>%
  pull(protein) %>% # extract protein name as character vector
  unique()

cat("\nProteins significant in at least one contrast: ", length(sig_proteins), "\n\n")

# Step 30: Shared heatmap setup

# Row ordering: proteins sorted by hallmark category then alphabetically
# within category, so related biology is visually grouped.
# Column ordering for individual heatmaps: participants sorted by study_group
# so groups are contiguous without hierarchical clustering of columns.
# Diverging color scale (blue-white-red) reflects z-score direction:
#   blue  = below cohort mean
#   white = at cohort mean
#   red   = above cohort mean

order_proteins_by_hallmark <- function(protein_set) {
  hallmark_lookup %>%
    filter(protein %in% protein_set) %>% # subset based on condition
    select(protein, hallmark_abbrev) %>% # subset based on name
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

sample_order <- cohort %>%
  arrange(study_group) %>%
  pull(participant_id)

col_annotation_individual <- data.frame(
  Group     = cohort$study_group,
  row.names = cohort$participant_id
)

heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# Group mean expression matrix — used in heatmaps B and C
group_means <- cohort %>%
  select(study_group, all_of(npqz_cols)) %>%
  group_by(study_group) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("study_group") %>%
  t() %>%
  as.data.frame()

rownames(group_means) <- protein_names

col_annotation_group <- data.frame(
  Group     = colnames(group_means),
  row.names = colnames(group_means)
)

# Step 31: Heatmap A — all 221 proteins, individual participants

# Provides the most granular view of the full protein set.
# Rows ordered by hallmark category, columns ordered by group.
# No hierarchical clustering of rows or columns.

row_order_all <- order_proteins_by_hallmark(protein_names)
row_ann_all   <- make_row_annotation(protein_names)

hallmark_ann_colors_all <- list(
  Hallmark = hallmark_colors[unique(row_ann_all$Hallmark)]
)

png(
  file.path(output_dir, "heatmap_A_all_proteins_individual.png"),
  width  = 16, height = 20, units = "in", res = 300
)

pheatmap(
  expr_matrix[row_order_all, sample_order],
  color             = heatmap_colors,
  cluster_rows      = FALSE, # no heirarchical clustering
  cluster_cols      = FALSE,
  annotation_col    = col_annotation_individual[sample_order, , drop = FALSE],
  annotation_row    = row_ann_all, # add color coded sidebar
  annotation_colors = c(group_colors, hallmark_ann_colors_all),
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = 5,
  fontsize          = 9,
  main              = paste0("All 221 Proteins — Individual Participants\n", data_note),
  border_color      = NA
)

dev.off() # closes the png file to signal it is time to save
cat("Saved: heatmap_A_all_proteins_individual.png\n")

# Ste`q 1p 32: Heatmap B — all 221 proteins, group means
#
# Reduces within-group noise by averaging across participants.
# Useful for comparing the overall protein expression pattern
# across groups without individual-level variability.

png(
  file.path(output_dir, "heatmap_B_all_proteins_group_means.png"),
  width  = 10, height = 20, units = "in", res = 300
)

pheatmap(
  as.matrix(group_means)[row_order_all, ],
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

dev.off()
cat("Saved: heatmap_B_all_proteins_group_means.png\n")

# Step 33: Heatmap C — significant proteins only, group means
#
# Restricts rows to proteins reaching FDR < 0.05 in at least one contrast.
# Height scales dynamically with the number of significant proteins.
# If no proteins reach significance (possible at n=66 with simulated data),
# the step is skipped with an informative message.

if (length(sig_proteins) > 0) {
  
  row_order_sig <- order_proteins_by_hallmark(sig_proteins)
  row_ann_sig   <- make_row_annotation(sig_proteins)
  
  hallmark_ann_colors_sig <- list(
    Hallmark = hallmark_colors[unique(row_ann_sig$Hallmark)]
  )
  
  plot_height_C <- max(6, length(sig_proteins) * 0.18 + 2)
  
  png(
    file.path(output_dir, "heatmap_C_sig_proteins_group_means.png"),
    width  = 10, height = plot_height_C, units = "in", res = 300
  )
  
  pheatmap(
    as.matrix(group_means)[row_order_sig, ],
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
      data_note
    ),
    border_color = NA
  )
  
  dev.off()
  cat("Saved: heatmap_C_sig_proteins_group_means.png\n")
  
} else {
  cat("No significant proteins at FDR < 0.05 — heatmap C skipped.\n")
  cat("This is expected with simulated data at n=66.\n")
}

# Step 34: Heatmap D — significant proteins only, individual participants
#
# Same protein set as heatmap C but at individual resolution.
# Allows inspection of within-group variability and potential
# participant subgroups that group means would obscure.

if (length(sig_proteins) > 0) {
  
  plot_height_D <- max(6, length(sig_proteins) * 0.18 + 2)
  
  png(
    file.path(output_dir, "heatmap_D_sig_proteins_individual.png"),
    width  = 16, height = plot_height_D, units = "in", res = 300
  )
  
  pheatmap(
    expr_matrix[row_order_sig, sample_order],
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
      data_note
    ),
    border_color = NA
  )
  
  dev.off()
  cat("Saved: heatmap_D_sig_proteins_individual.png\n")
  
} else {
  cat("No significant proteins at FDR < 0.05 — heatmap D skipped.\n")
}

# Step 35: Export results to Excel
#
# One sheet per contrast for pairwise results, plus the global F-test
# and hallmark enrichment results. All sheets include hallmark annotations
# so results are interpretable without cross-referencing other files.
# Rows within each contrast sheet are sorted by adjusted p-value.

for (cn in contrast_names) {
  sheet_name <- substr(cn, 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(
    wb, sheet_name,
    pairwise_results[[cn]] %>%
      select(protein, hallmark, hallmark_abbrev, logFC, AveExpr,
             t, P.Value, adj.P.Val, sig, direction) %>%
      arrange(adj.P.Val)
  )
}

addWorksheet(wb, "Hallmark_Enrichment")
writeData(
  wb, "Hallmark_Enrichment",
  hallmark_enrichment %>%
    arrange(contrast, p_hypergeom) %>%
    select(contrast, hallmark_abbrev, n_category, n_sig,
           pct_sig, p_hypergeom, p_adj_BH)
)

saveWorkbook(
  wb,
  file.path(output_dir, "limma_proteomics_results.xlsx"),
  overwrite = TRUE
)

cat("\nAll outputs saved to:", output_dir, "\n")
cat("\nFiles produced:\n")
cat("  Volcano plots:\n")
cat("    volcano_Acute_IPV_vs_NonIPV.png\n")
cat("    volcano_Chronic_vs_Acute_IPV.png\n")
cat("    volcano_Chronic_IPV_vs_NonIPV.png\n")
cat("  Heatmaps:\n")
cat("    heatmap_A_all_proteins_individual.png\n")
cat("    heatmap_B_all_proteins_group_means.png\n")
cat("    heatmap_C_sig_proteins_group_means.png  (skipped if no sig proteins)\n")
cat("    heatmap_D_sig_proteins_individual.png   (skipped if no sig proteins)\n")
cat("  Tables:\n")
cat("    limma_proteomics_results.xlsx\n")