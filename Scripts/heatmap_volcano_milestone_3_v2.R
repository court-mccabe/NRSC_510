# Milestone 3: NULISA Proteomics — Limma, Mann-Whitney U, Volcano Plots, Heatmaps
# NPQ_ columns are log2-scale (Alamar ARGO pipeline output). No further
# transformation needed. logFC = log2 fold-change between group means.
# Contrasts: (1) Acute IPV-BI vs Acute Non-IPV-BI [primary]
#            (2) Chronic IPV-BI vs Acute IPV-BI    [secondary]
#            (3) Chronic IPV-BI vs Acute Non-IPV-BI [tertiary]


# Step 16: Packages
library(tidyverse)
library(limma)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)


# Step 17: Setup
setwd("C:/NRSC_510")
output_dir <- "C:/NRSC_510/Outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

run_label  <- "_v9"          # increment each run to avoid overwriting outputs
cohort     <- readRDS("ipv_cohort_clean.rds")
data_note  <- "Simulated data v7 — raw NPQ log2-scale"

qc_meta_cols <- c("nc_mean_npq", "nc_sd_npq", "nc_pass",
                  "ic_reads", "ic_pct_of_median", "ic_pass",
                  "ipc_cv_pct", "ipc_pass")


# Step 18: Data format check
# Halts if NPQ_ values look like z-scores instead of log2-scale NPQ.
npq_cols_all <- cohort %>%
  select(starts_with("NPQ_"), -any_of(qc_meta_cols)) %>% names()

npq_vals <- cohort %>% select(all_of(npq_cols_all)) %>%
  as.matrix() %>% as.vector() %>% na.omit()

if (abs(mean(npq_vals)) < 2 || sd(npq_vals) < 0.5 || min(npq_vals) < 0)
  stop("NPQ values do not look like log2-scale data (expected mean ~13, range ~6-20).\n",
       "Rebuild ipv_cohort_clean.rds from IPV_cohort_simulated_v7.xlsx.")


# Step 19: QC checks — NC, IC, IPC
# Thresholds from Alamar NULISAseqR User Guide, Chapter 2.
if (all(qc_meta_cols %in% names(cohort))) {
  if (any(!cohort$nc_pass,  na.rm = TRUE))
    warning("NC check: nc_mean_npq >= 3 on one or more plates.")
  if (any(!cohort$ic_pass,  na.rm = TRUE))
    warning(sum(!cohort$ic_pass, na.rm = TRUE), " sample(s) failed IC check: ",
            paste(cohort$participant_id[!cohort$ic_pass], collapse = ", "))
  if (any(!cohort$ipc_pass, na.rm = TRUE))
    warning("IPC CV check: ipc_cv_pct > 10% on one or more plates.")
  cohort <- cohort %>% mutate(qc_pass_all = nc_pass & ic_pass & ipc_pass)
} else {
  warning("QC columns missing — rebuild from IPV_cohort_simulated_v7.xlsx.")
}


# Step 20: Hallmark lookup
hallmark_lookup <- read_csv("C:/NRSC_510/Neuro_220_Panel.csv", show_col_types = FALSE) %>%
  rename(target_full = `Target List`,
         uniprot_id  = `Uniprot ID`,
         hallmark    = `Hallmarks of Neurodegenerative Disease`) %>%
  mutate(
    protein = target_full %>%
      str_replace_all("\\s*\\(.*?\\)", "") %>%
      str_replace_all(";\\s*", "_") %>%
      str_replace_all("\\|.*", "") %>%
      str_replace_all("\\s+", "-") %>%
      str_trim(),
    protein  = if_else(str_detect(target_full, "IL-12A"), "IL-12p70", protein),
    hallmark = if_else(hallmark == "Aberrant proteostasis", "Aberrant Proteostasis", hallmark),
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


# Step 21: Identify and validate protein columns
npq_cols     <- npq_cols_all[!str_remove(npq_cols_all, "^NPQ_") %in% c("MSLN")]
protein_names <- str_remove(npq_cols, "^NPQ_")

unmatched <- protein_names[!protein_names %in% hallmark_lookup$protein]
if (length(unmatched) > 0)
  stop("Proteins with no hallmark match: ", paste(unmatched, collapse = ", "))


# Step 22: Expression and design matrices
expr_matrix <- cohort %>% select(all_of(npq_cols)) %>% as.matrix() %>% t()
rownames(expr_matrix) <- protein_names
colnames(expr_matrix) <- cohort$participant_id

design <- model.matrix(
  ~ 0 + study_group + age_years + sex,
  data = cohort %>% mutate(sex = ifelse(sex == "Male", 1L, 0L))
)
colnames(design) <- make.names(colnames(design))


# Step 23: Contrasts
contrast_matrix <- makeContrasts(
  Acute_IPV_vs_NonIPV   = study_groupAcute.IPV.BI   - study_groupAcute.Non.IPV.BI,
  Chronic_vs_Acute_IPV  = study_groupChronic.IPV.BI - study_groupAcute.IPV.BI,
  Chronic_IPV_vs_NonIPV = study_groupChronic.IPV.BI - study_groupAcute.Non.IPV.BI,
  levels = design
)
contrast_labels <- c(
  Acute_IPV_vs_NonIPV   = "Acute IPV-BI vs Acute Non-IPV-BI",
  Chronic_vs_Acute_IPV  = "Chronic IPV-BI vs Acute IPV-BI",
  Chronic_IPV_vs_NonIPV = "Chronic IPV-BI vs Acute Non-IPV-BI"
)
contrast_names  <- names(contrast_labels)
contrast_groups <- list(
  Acute_IPV_vs_NonIPV   = c("Acute IPV-BI",   "Acute Non-IPV-BI"),
  Chronic_vs_Acute_IPV  = c("Chronic IPV-BI", "Acute IPV-BI"),
  Chronic_IPV_vs_NonIPV = c("Chronic IPV-BI", "Acute Non-IPV-BI")
)


# Step 24: Limma fit and BH pre-check
fit  <- lmFit(expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

# BH pre-check: halt if correction collapses on primary contrast
tt_check  <- topTable(fit2, coef = "Acute_IPV_vs_NonIPV", number = Inf,
                      adjust = "BH", sort.by = "none")
if (sum(tt_check$adj.P.Val >= 0.999, na.rm = TRUE) > 0.95 * nrow(tt_check))
  stop("BH-adjusted p-values are almost all 1.0 — check NPQ data format (Step 18).")


# Step 25: Extract limma results
limma_results <- lapply(contrast_names, function(cn) {
  topTable(fit2, coef = cn, number = Inf, adjust = "BH", sort.by = "P") %>%
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
    left_join(hallmark_lookup %>% select(protein, hallmark, hallmark_abbrev), by = "protein") %>%
    rename(limma_logFC = logFC, limma_p_raw = P.Value, limma_adj_p = adj.P.Val)
}) %>% setNames(contrast_names)


# Step 26: Mann-Whitney U
# Nonparametric sanity check; cannot adjust for covariates.
# mw_adj_p = BH across all 663 tests (221 proteins x 3 contrasts).
run_mw <- function(cn) {
  g   <- contrast_groups[[cn]]
  map_dfr(seq_along(npq_cols), function(i) {
    x1 <- cohort[[npq_cols[i]]][cohort$study_group == g[1]] %>% na.omit()
    x2 <- cohort[[npq_cols[i]]][cohort$study_group == g[2]] %>% na.omit()
    if (length(x1) < 3 || length(x2) < 3)
      return(tibble(protein = protein_names[i], mw_p_raw = NA_real_, mw_logFC = NA_real_))
    tibble(protein  = protein_names[i],
           mw_p_raw = wilcox.test(x1, x2, exact = FALSE)$p.value,
           mw_logFC = mean(x1) - mean(x2))
  }) %>% mutate(contrast = cn)
}

mw_all <- map_dfr(contrast_names, run_mw) %>%
  mutate(mw_adj_p     = p.adjust(mw_p_raw, method = "BH"),
         mw_sig       = mw_adj_p < 0.05,
         mw_direction = case_when(
           mw_sig & mw_logFC > 0 ~ "Up",
           mw_sig & mw_logFC < 0 ~ "Down",
           TRUE                  ~ "NS"
         ))


# Step 27: Merge limma and Mann-Whitney results
pairwise_results <- lapply(contrast_names, function(cn) {
  left_join(limma_results[[cn]],
            mw_all %>% filter(contrast == cn) %>%
              select(protein, mw_p_raw, mw_adj_p, mw_logFC, mw_sig, mw_direction),
            by = "protein")
}) %>% setNames(contrast_names)


# Step 28: Write p-values back to RDS
contrast_short <- c(Acute_IPV_vs_NonIPV   = "AcuteIPVvsNon",
                    Chronic_vs_Acute_IPV  = "ChronicvsAcute",
                    Chronic_IPV_vs_NonIPV = "ChronicvsNon")

for (cn in contrast_names) {
  res <- pairwise_results[[cn]] %>% select(protein, limma_adj_p, mw_adj_p, limma_sig)
  for (i in seq_len(nrow(res))) {
    p <- res$protein[i]; s <- contrast_short[cn]
    cohort[[paste0("limma_adj_p_", s, "_", p)]] <- res$limma_adj_p[i]
    cohort[[paste0("mw_adj_p_",    s, "_", p)]] <- res$mw_adj_p[i]
    cohort[[paste0("limma_sig_",   s, "_", p)]] <- res$limma_sig[i]
  }
}
saveRDS(cohort, "ipv_cohort_clean.rds")


# Step 29: Shared visual elements
hallmark_levels <- sort(unique(hallmark_lookup$hallmark_abbrev))
hallmark_colors <- setNames(
  colorRampPalette(brewer.pal(9, "Set3"))(length(hallmark_levels)),
  hallmark_levels
)
group_colors <- list(Group = c(
  "Acute IPV-BI"     = "#E41A1C",
  "Acute Non-IPV-BI" = "#377EB8",
  "Chronic IPV-BI"   = "#984EA3"
))


# Step 30: Volcano plot function
make_volcano <- function(result_df, contrast_label, sig_threshold = 0.05,
                         exploratory = FALSE) {
  df <- result_df %>%
    mutate(sig             = limma_adj_p < sig_threshold,
           hallmark_abbrev = factor(hallmark_abbrev, levels = hallmark_levels))
  
  labels <- if (any(df$sig)) slice_min(filter(df, sig), limma_adj_p, n = 5, with_ties = FALSE) else
    slice_min(df, limma_adj_p, n = 5, with_ties = FALSE)
  
  y_max <- max(max(-log10(df$limma_adj_p), na.rm = TRUE) * 1.10,
               -log10(sig_threshold) + 1.0)
  
  ggplot(df, aes(x = limma_logFC, y = -log10(limma_adj_p))) +
    geom_point(data = filter(df, !sig), color = "grey75", alpha = 0.5, size = 1.5) +
    geom_point(data = filter(df,  sig), aes(color = hallmark_abbrev), alpha = 0.9, size = 2.2) +
    geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    geom_text_repel(data = labels, aes(label = protein), size = 2.8,
                    max.overlaps = 20, segment.color = "grey50",
                    box.padding = 0.4, point.padding = 0.3, color = "black") +
    scale_color_manual(values = hallmark_colors, name = "Hallmark",
                       na.value = "grey75", drop = TRUE) +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(
      title    = contrast_label,
      subtitle = paste0(if (exploratory) "Exploratory: " else "",
                        "BH-adjusted p < ", sig_threshold, " — ", data_note),
      x = "log2 fold-change (NPQ units)",
      y = expression(-log[10](BH-adjusted~p))
    ) +
    theme_classic(base_size = 10) +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 10),
          plot.subtitle   = element_text(color = "grey40", hjust = 0.5, size = 7),
          legend.text     = element_text(size = 7),
          legend.title    = element_text(face = "bold", size = 8),
          legend.key.size = unit(0.35, "cm"))
}


# Step 31: P-value histograms
# One function produces both raw and BH-adjusted histogram for a given contrast.
save_pval_hist <- function(cn) {
  tt <- topTable(fit2, coef = cn, number = Inf, adjust = "BH", sort.by = "none")
  
  make_hist <- function(vals, x_lab, fill_col, threshold, caption_text) {
    ggplot(data.frame(p = vals), aes(x = p)) +
      geom_histogram(bins = 50, fill = fill_col, color = "white", alpha = 0.8,
                     linewidth = 0.3) +
      geom_vline(xintercept = threshold, linetype = "dashed",
                 color = "#E41A1C", linewidth = 0.5) +
      labs(title    = contrast_labels[cn],
           subtitle = data_note,
           x = x_lab, y = "Number of proteins",
           caption = caption_text) +
      theme_classic(base_size = 10) +
      theme(plot.title   = element_text(face = "bold", hjust = 0.5, size = 10),
            plot.subtitle = element_text(color = "grey40", hjust = 0.5, size = 7),
            plot.caption  = element_text(color = "grey50", size = 8))
  }
  
  p_raw <- make_hist(tt$P.Value,   "Raw p-value",        "grey70", 0.05,
                     paste0("p < 0.05: ", sum(tt$P.Value < 0.05), " proteins"))
  p_adj <- make_hist(tt$adj.P.Val, "BH-adjusted p-value","#377EB8", 0.05,
                     paste0("FDR < 0.05: ", sum(tt$adj.P.Val < 0.05),
                            "  |  FDR < 0.20: ", sum(tt$adj.P.Val < 0.20), " proteins"))
  
  ggsave(file.path(output_dir, paste0("pval_raw_", cn, run_label, ".png")),
         p_raw, width = 8, height = 5, dpi = 300)
  ggsave(file.path(output_dir, paste0("pval_adj_", cn, run_label, ".png")),
         p_adj, width = 8, height = 5, dpi = 300)
}

walk(contrast_names, save_pval_hist)


# Step 32: Volcano plots
save_volcano <- function(cn, threshold, exploratory) {
  p <- make_volcano(pairwise_results[[cn]], contrast_labels[cn],
                    sig_threshold = threshold, exploratory = exploratory)
  tag  <- if (exploratory) "volcano_exploratory_" else "volcano_primary_"
  ggsave(file.path(output_dir, paste0(tag, cn, run_label, ".png")),
         p, width = 10, height = 7, dpi = 300)
}

walk(contrast_names, save_volcano, threshold = 0.05, exploratory = FALSE)
walk(contrast_names, save_volcano, threshold = 0.20, exploratory = TRUE)

sig_proteins <- bind_rows(pairwise_results) %>%
  filter(limma_sig) %>% pull(protein) %>% unique()


# Step 32b: Volcano plots coloured by study group
# Significant points coloured by which group is upregulated; NS remains grey.
# Colours match group_colors defined in Step 29.

#make_volcano_group <- function(result_df, contrast_label, cn,
 #                              sig_threshold = 0.05, exploratory = FALSE) {
  #g      <- contrast_groups[[cn]]
  #up_col <- group_colors$Group[g[1]]
  #dn_col <- group_colors$Group[g[2]]
  
  # Use completely fresh vectors to avoid any column name conflicts
  #is_sig   <- result_df$limma_sig
  #is_up    <- is_sig & result_df$limma_logFC > 0
  #is_down  <- is_sig & result_df$limma_logFC < 0
  #pt_color <- ifelse(is_up, "Up", ifelse(is_down, "Down", "NS"))
  
  #plot_df <- data.frame(
  #  protein    = result_df$protein,
  #  fc         = result_df$limma_logFC,
   # adjp       = result_df$limma_adj_p,
   # pt_color   = pt_color,
  #  stringsAsFactors = FALSE
#  )
  
 # labels  <- if (any(is_sig)) plot_df[is_sig, ][order(plot_df[is_sig, "adjp"])[1:min(5, sum(is_sig))], ] else
  #  plot_df[order(plot_df$adjp)[1:5], ]
  
  #y_max   <- max(max(-log10(plot_df$adjp), na.rm = TRUE) * 1.10, -log10(0.05) + 1.0)
  
  #ggplot(plot_df, aes(x = fc, y = -log10(adjp), color = pt_color)) +
   # geom_point(alpha = 0.7, size = 1.8) +
  #  geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed",
  #             color = "grey40", linewidth = 0.5) +
  #  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
  #             color = "grey40", linewidth = 0.5) +
  #  geom_text_repel(data = labels, aes(label = protein), size = 2.8,
  #                  max.overlaps = 20, segment.color = "grey50",
  #                  box.padding = 0.4, point.padding = 0.3, color = "black") +
  #  scale_color_manual(
#      values = c("Up" = up_col, "Down" = dn_col, "NS" = "grey75"),
 #     limits = c("Up", "Down", "NS"),
#      labels = c("Up" = g[1], "Down" = g[2], "NS" = "NS"),
 #     name   = "Group",
  #    drop   = FALSE
   # ) +
  #  coord_cartesian(ylim = c(0, y_max)) +
  #  labs(
  #    title    = contrast_label,
   #   subtitle = paste0(if (exploratory) "Exploratory: " else "",
    #                    "BH-adjusted p < ", sig_threshold, " — ", data_note),
     # x = "log2 fold-change (NPQ units)",
      #y = expression(-log[10](BH-adjusted~p))
  #  ) +
   # theme_classic(base_size = 10) +
  #  theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 10),
   #       plot.subtitle   = element_text(color = "grey40", hjust = 0.5, size = 7),
  #        legend.text     = element_text(size = 7),
  #        legend.title    = element_text(face = "bold", size = 8),
  #        legend.key.size = unit(0.35, "cm"))
#}

#for (cn in contrast_names) {
 # p <- make_volcano_group(pairwise_results[[cn]], contrast_labels[cn], cn,
  #                        sig_threshold = 0.05, exploratory = FALSE)
  #ggsave(file.path(output_dir, paste0("volcano_group_", cn, run_label, ".png")),
   #      p, width = 10, height = 7, dpi = 300)
#}

# Step 33: Heatmap setup
save_heatmap <- function(filepath, width, height, pheatmap_args) {
  png(filepath, width = width, height = height, units = "in", res = 300)
  do.call(pheatmap, pheatmap_args)
  dev.off()
}

order_by_hallmark <- function(proteins) {
  hallmark_lookup %>% filter(protein %in% proteins) %>%
    distinct(protein, hallmark_abbrev) %>%
    arrange(hallmark_abbrev, protein) %>% pull(protein)
}

row_annotation <- function(proteins) {
  hallmark_lookup %>% filter(protein %in% proteins) %>%
    distinct(protein, hallmark_abbrev) %>%
    column_to_rownames("protein") %>% rename(Hallmark = hallmark_abbrev)
}

sample_order <- cohort %>% arrange(study_group) %>% pull(participant_id)

col_ann_individual <- data.frame(Group = cohort$study_group,
                                 row.names = cohort$participant_id)

group_means <- cohort %>% select(study_group, all_of(npq_cols)) %>%
  group_by(study_group) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("study_group") %>% t() %>% as.data.frame()
rownames(group_means) <- protein_names

col_ann_group  <- data.frame(Group = colnames(group_means), row.names = colnames(group_means))
heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)


# Step 34: Heatmap builder
# Single function handles all four heatmap variants.
make_heatmap <- function(proteins, mat, ann_col, filename, width, fontsize_row) {
  row_ord  <- order_by_hallmark(proteins)
  row_ann  <- row_annotation(proteins)
  ann_cols <- c(group_colors, list(Hallmark = hallmark_colors[unique(row_ann$Hallmark)]))
  height   <- max(6, length(proteins) * 0.18 + 2)
  
  save_heatmap(
    filepath = file.path(output_dir, filename),
    width = width, height = height,
    pheatmap_args = list(
      mat               = mat[row_ord, ],
      color             = heatmap_colors,
      cluster_rows      = FALSE,
      cluster_cols      = FALSE,
      annotation_col    = ann_col,
      annotation_row    = row_ann,
      annotation_colors = ann_cols,
      show_rownames     = TRUE,
      show_colnames     = width < 12,
      fontsize_row      = fontsize_row,
      fontsize          = 9,
      main              = paste0(filename %>% str_remove(".png") %>%
                                   str_replace_all("_", " "), "\n", data_note),
      border_color      = NA
    )
  )
}

graphics.off()

# A: all proteins — individual participants
make_heatmap(protein_names, expr_matrix[, sample_order],
             col_ann_individual[sample_order, , drop = FALSE],
             paste0("heatmap_A_all_proteins_individual", run_label, ".png"),
             width = 16, fontsize_row = 5)

# B: all proteins — group means
make_heatmap(protein_names, as.matrix(group_means),
             col_ann_group,
             paste0("heatmap_B_all_proteins_group_means", run_label, ".png"),
             width = 10, fontsize_row = 5)

# C & D: significant proteins only (skipped if none reach FDR < 0.05)
if (length(sig_proteins) > 0) {
  make_heatmap(sig_proteins, as.matrix(group_means),
               col_ann_group,
               paste0("heatmap_C_sig_proteins_group_means", run_label, ".png"),
               width = 10, fontsize_row = 7)
  make_heatmap(sig_proteins, expr_matrix[, sample_order],
               col_ann_individual[sample_order, , drop = FALSE],
               paste0("heatmap_D_sig_proteins_individual", run_label, ".png"),
               width = 16, fontsize_row = 7)
}


# Step 35: Export results CSVs
walk(contrast_names, function(cn) {
  pairwise_results[[cn]] %>%
    select(protein, hallmark, hallmark_abbrev,
           limma_logFC, limma_p_raw, limma_adj_p, limma_sig, limma_direction,
           mw_logFC,   mw_p_raw,   mw_adj_p,   mw_sig,   mw_direction) %>%
    arrange(limma_adj_p) %>%
    write_csv(file.path(output_dir, paste0("results_", cn, run_label, ".csv")))
})
saveRDS(pairwise_results, "pairwise_results.rds")
