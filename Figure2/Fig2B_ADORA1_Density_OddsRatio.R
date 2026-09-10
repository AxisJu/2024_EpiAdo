# ==============================================================================
# Script: Fig2B_ADORA1_Density_OddsRatio.R
# Purpose: ADORA1 single-cell density distribution, Limma log2FC CI, and Logistic Regression OR
# Article: EpiAdo (Science Translational Medicine) - Figure 2B
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(grid)
  library(gridExtra)
})

out_dir_sd <- "../sourcedata/Figure2"
out_dir_plot <- "../sourcedata/check_plots"
dir.create(out_dir_sd, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_plot, showWarnings = FALSE, recursive = TRUE)

cat("Loading Fig 2B input files...\n")
# [DATA REQUIRED] The line below originally loaded 'Glu.N ADORA1 Threshold_OR.csv' from a private path.
# Please provide 'threshold_or' by loading the appropriate data before running this script.
# Example: threshold_or <- read.csv('path/to/Glu.N ADORA1 Threshold_OR.csv')
if ("X" %in% colnames(threshold_or)) threshold_or <- threshold_or[, !colnames(threshold_or) %in% c("X", "Unnamed..0")]

# [DATA REQUIRED] The line below originally loaded 'Glu.N Adjusted CI of Gene&AdoScore FC.csv' from a private path.
# Please provide 'de_ci' by loading the appropriate data before running this script.
# Example: de_ci <- read.csv('path/to/Glu.N Adjusted CI of Gene&AdoScore FC.csv')
adora1_ci <- de_ci[de_ci[, 1] == "ADORA1", ]

cat("Loading single-cell ADORA1 expression...\n")
# [DATA REQUIRED] The line below originally loaded '20260309_glu_n_seurat.qs' from a private path.
# Please provide 'glu_n_seurat' by loading the appropriate data before running this script.
# Example: glu_n_seurat <- qread('path/to/20260309_glu_n_seurat.qs')
exp_adora1 <- GetAssayData(glu_n_seurat, layer = "data")["ADORA1", ]
meta <- glu_n_seurat@meta.data

sc_df <- data.frame(
  CellID = rownames(meta),
  Group = meta$Group,
  Cohort = meta$Cohort,
  ADORA1_LogCPM = as.numeric(exp_adora1)
) %>%
  # NOTE: Cohorts "Tran", "YC", and "LWS" are excluded from the density plot because
  # they use a different sequencing platform or library preparation protocol, making
  # direct comparison of CPM-scale ADORA1 expression distributions non-comparable.
  # These cohorts are included in pseudobulk statistical analyses (Fig2A, Fig2C) with
  # appropriate batch correction. See Methods: "Single-nucleus RNA-seq data integration".
  dplyr::filter(!Cohort %in% c("Tran", "YC", "LWS")) %>%
  dplyr::filter(ADORA1_LogCPM > 0)

# 1. Export SourceData
cat("Exporting SourceData for Figure 2B...\n")
write.csv(threshold_or, file.path(out_dir_sd, "Fig2B_ADORA1_Threshold_OddsRatio.csv"), row.names = FALSE)

sd_ci <- data.frame(
  Gene = "ADORA1",
  Beta_log2FC = adora1_ci$beta,
  FoldChange = 2^(adora1_ci$beta),
  CI_Lower = adora1_ci$lower.un,
  CI_Upper = adora1_ci$upper.un,
  P_val_adj = adora1_ci$p.bh
)
write.csv(sd_ci, file.path(out_dir_sd, "Fig2B_ADORA1_Limma_FoldChange_CI.csv"), row.names = FALSE)

# Density points for reproducibility
dens_ctrl <- density(sc_df$ADORA1_LogCPM[sc_df$Group == "Control"])
dens_epi <- density(sc_df$ADORA1_LogCPM[sc_df$Group == "Epilepsy"])
sd_dens <- rbind(
  data.frame(Group = "Control", ADORA1_LogCPM = dens_ctrl$x, Density = dens_ctrl$y),
  data.frame(Group = "Epilepsy", ADORA1_LogCPM = dens_epi$x, Density = dens_epi$y)
)
write.csv(sd_dens, file.path(out_dir_sd, "Fig2B_ADORA1_Density_Distribution.csv"), row.names = FALSE)

# 2. Generate Check Plots
cat("Generating Check Plot for Figure 2B...\n")
color_group <- c("Control" = "#E5D2C2", "Epilepsy" = "#205A9E")

# Top panel: Limma FC CI
p_top <- ggplot(adora1_ci, aes(x = beta, y = 1)) +
  geom_point(size = 2.5, color = "black") +
  geom_errorbar(aes(xmin = lower.un, xmax = upper.un, y = 1), orientation = "y", width = 0.4, linewidth = 0.8, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_x_continuous(limits = c(-0.1, 0.1), breaks = c(-0.1, 0, 0.1)) +
  scale_y_continuous(limits = c(0, 2), breaks = NULL) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(5, 10, 0, 10)
  ) +
  labs(x = expression(log[2]~FC), title = expression("FC = 0.95, FDR < 5% (" * italic("limma") * ")"))

# Main panel: Density & OR
p_main <- ggplot(sc_df, aes(x = ADORA1_LogCPM, color = Group, fill = Group)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_color_manual(values = color_group) +
  scale_fill_manual(values = color_group) +
  scale_x_log10(
    limits = c(0.05, 10),
    breaks = c(0.1, 0.5, 1, 5, 10),
    labels = c(expression(10^-1), "0.5", expression(10^0), "5", expression(10^1))
  ) +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "red", linewidth = 0.8) +
  geom_point(data = filter(threshold_or, p < 0.05), aes(x = threshold, y = or), inherit.aes = FALSE, color = "black", size = 2) +
  geom_errorbar(data = filter(threshold_or, p < 0.05), aes(x = threshold, y = or, ymin = ci.L, ymax = ci.H), inherit.aes = FALSE, width = 0.05, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  scale_y_continuous(
    name = "Density",
    limits = c(0, 2.8),
    sec.axis = sec_axis(~ ., name = "Odds Ratio (OR)", breaks = c(0, 0.5, 1, 1.5, 2))
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.85, 0.85),
    legend.title = element_blank(),
    plot.margin = margin(0, 10, 10, 10)
  ) +
  annotate("text", x = 0.15, y = 1.2, label = "italic(ADORA1)^Low", parse = TRUE, color = "#7BC5B4", fontface = "bold", size = 4) +
  annotate("text", x = 2.5, y = 1.2, label = "italic(ADORA1)^High", parse = TRUE, color = "#205A9E", fontface = "bold", size = 4) +
  labs(x = "ADORA1 Expression (LogCPM)")

png(file.path(out_dir_plot, "check_Fig2B_ADORA1_Density_OddsRatio.png"), width = 2100, height = 2400, res = 300)
grid.arrange(p_top, p_main, ncol = 1, heights = c(1, 4))
dev.off()

cat("Figure 2B finished successfully.\n")
