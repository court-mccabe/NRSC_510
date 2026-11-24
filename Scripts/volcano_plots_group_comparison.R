# Volcano plots coloured by study group
# Standalone: loads pairwise_results.rds saved at the end of milestone 3.

library(ggplot2)
library(ggrepel)

setwd("C:/NRSC_510")
output_dir       <- "C:/NRSC_510/Outputs"
data_note        <- "Simulated data v7 — raw NPQ log2-scale"
pairwise_results <- readRDS("pairwise_results.rds")

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

run_label <- "_v15"   # increment each run

group_colors <- c(
  "Acute IPV-BI"     = "#E41A1C",
  "Acute Non-IPV-BI" = "#377EB8",
  "Chronic IPV-BI"   = "#984EA3"
)

for (cn in contrast_names) {
  g      <- contrast_groups[[cn]]
  up_col <- unname(group_colors[g[1]])
  dn_col <- unname(group_colors[g[2]])
  
  res      <- pairwise_results[[cn]]
  is_up    <- res$limma_sig & res$limma_logFC > 0
  is_down  <- res$limma_sig & res$limma_logFC < 0
  pt_color <- factor(
    ifelse(is_up, "Up", ifelse(is_down, "Down", "NS")),
    levels = c("Up", "Down", "NS")
  )
  
  plot_df <- data.frame(
    protein  = res$protein,
    fc       = res$limma_logFC,
    adjp     = res$limma_adj_p,
    pt_color = pt_color
  )
  
  sig_df    <- subset(plot_df, pt_color != "NS")
  labels_df <- sig_df[order(sig_df$adjp)[1:min(5, nrow(sig_df))], ]
  y_max     <- max(max(-log10(plot_df$adjp), na.rm = TRUE) * 1.10,
                   -log10(0.05) + 1.0)
  
  p <- ggplot(plot_df, aes(x = fc, y = -log10(adjp), color = pt_color)) +
    geom_point(data = subset(plot_df, pt_color == "NS"), alpha = 0.4, size = 1.5) +
    geom_point(data = subset(plot_df, pt_color != "NS"), alpha = 1.0, size = 3.0) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    geom_text_repel(data = labels_df, aes(label = protein), size = 2.8,
                    max.overlaps = 20, segment.color = "grey50",
                    box.padding = 0.4, point.padding = 0.3, color = "black") +
    scale_color_manual(
      values = c("Up" = up_col, "Down" = dn_col, "NS" = "grey75"),
      labels = c("Up" = g[1], "Down" = g[2], "NS" = "NS"),
      name   = "Group"
    ) +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(
      title    = contrast_labels[cn],
      subtitle = paste0("BH-adjusted p < 0.05 — ", data_note),
      x        = "log2 fold-change (NPQ units)",
      y        = expression(-log[10](BH-adjusted~p))
    ) +
    theme_classic(base_size = 10) +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 10),
          plot.subtitle   = element_text(color = "grey40", hjust = 0.5, size = 7),
          legend.text     = element_text(size = 7),
          legend.title    = element_text(face = "bold", size = 8),
          legend.key.size = unit(0.35, "cm"))
  
  ggsave(file.path(output_dir, paste0("volcano_group_", cn, run_label, ".png")),
         p, width = 10, height = 7, dpi = 300)
}