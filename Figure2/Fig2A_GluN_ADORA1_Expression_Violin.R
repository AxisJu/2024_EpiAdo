# ==============================================================================
# Script: Fig2A_GluN_ADORA1_Expression_Violin.R
# Purpose: ADORA1 expression and distribution in Glu.N across Control and Epilepsy groups
# Article: EpiAdo (Science Translational Medicine) - Figure 2A
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(gghalves)
  library(rstatix)
  library(ggpubr)
})

out_dir_sd <- "../sourcedata/Figure2"
out_dir_plot <- "../sourcedata/check_plots"
dir.create(out_dir_sd, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_plot, showWarnings = FALSE, recursive = TRUE)

cat("Loading glu_n_seurat...\n")
# [DATA REQUIRED] The line below originally loaded '20260309_glu_n_seurat.qs' from a private path.
# Please provide 'glu_n_seurat' by loading the appropriate data before running this script.
# Example: glu_n_seurat <- qread('path/to/20260309_glu_n_seurat.qs')
exp_adora1 <- GetAssayData(glu_n_seurat, layer = "data")["ADORA1", ]
meta <- glu_n_seurat@meta.data

# Single-cell data frame
sc_df <- data.frame(
  CellID = rownames(meta),
  Sample = meta$Sample_snSeq.processed,
  Group = meta$Group,
  BrainRegion = meta$BrainRegion,
  Cohort = meta$Cohort,
  ADORA1_LogCPM = as.numeric(exp_adora1)
)

# Pseudobulk per sample (for non-zero expression or all cells)
pb_df <- sc_df %>%
  dplyr::group_by(Sample, Group, BrainRegion, Cohort) %>%
  dplyr::summarise(
    CellCount = dplyr::n(),
    MeanExp_All = mean(ADORA1_LogCPM, na.rm = TRUE),
    MeanExp_NonZero = mean(ADORA1_LogCPM[ADORA1_LogCPM > 0], na.rm = TRUE),
    PercentExpressed = mean(ADORA1_LogCPM > 0) * 100,
    .groups = "drop"
  )

# Statistics
wt <- wilcox.test(ADORA1_LogCPM ~ Group, data = sc_df)
p_val_sc <- wt$p.value

cat(sprintf("Single-cell Wilcoxon test p-value: %e\n", p_val_sc))

# 1. Export SourceData
cat("Exporting SourceData for Figure 2A...\n")
# Export subsample of single cells (1 in 5) for manageable file size while preserving full distribution
write.csv(sc_df[seq(1, nrow(sc_df), by = 5), ], 
          file.path(out_dir_sd, "Fig2A_GluN_ADORA1_SingleCell_Expression_Subsample.csv"), row.names = FALSE)
write.csv(pb_df, file.path(out_dir_sd, "Fig2A_GluN_ADORA1_Pseudobulk_Expression.csv"), row.names = FALSE)

# 2. Generate Check Plot
cat("Generating Check Plot for Figure 2A...\n")
color_group <- c("Control" = "#989798", "Epilepsy" = "#205A9E")

# Filter for plotting distribution
plot_data <- sc_df %>%
  dplyr::filter(ADORA1_LogCPM > 0)

p_violin <- ggplot(plot_data, aes(x = Group, y = ADORA1_LogCPM, fill = Group, color = Group)) +
  geom_half_violin(
    side = "l", scale = "width", trim = FALSE, alpha = 0.4, linewidth = 0.6
  ) +
  geom_half_boxplot(
    side = "l", width = 0.25, outlier.shape = NA, alpha = 0.8, color = "black", linewidth = 0.6
  ) +
  geom_point(
    aes(x = as.numeric(factor(Group)) + 0.15),
    position = position_jitter(width = 0.08, height = 0),
    shape = 16, size = 0.5, alpha = 0.15
  ) +
  scale_fill_manual(values = color_group) +
  scale_color_manual(values = color_group) +
  scale_y_continuous(limits = c(0, 4), breaks = 0:4) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  annotate("text", x = 1.5, y = 3.8, label = "P < 1 x 10^-100", size = 4.5, fontface = "italic") +
  labs(
    title = "ADORA1 expression and distribution in Glu.N",
    x = NULL,
    y = "ADORA1 Expression (LogCPM)"
  )

ggsave(file.path(out_dir_plot, "check_Fig2A_GluN_ADORA1_Expression_Violin.png"), plot = p_violin, width = 4.5, height = 6, dpi = 300)
cat("Figure 2A finished successfully.\n")
