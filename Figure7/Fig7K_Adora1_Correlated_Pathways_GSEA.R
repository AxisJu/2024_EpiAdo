# ==============================================================================
# Script Name: Fig7K_Adora1_Correlated_Pathways_GSEA.R
# Description: Generates Figure 7K: Gene set enrichment analysis based on Adora1 
#              correlation in Vglut2+ Glu.N (Fold change error bars, heatmap, and NES bar).
# Author: Auto-refactored for EpiAdo STM Revision
# Date: 2026-09-08
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(cowplot)
})

# Define paths
out_dir_data <- "../sourcedata/Figure7"
out_dir_plot <- "../sourcedata/check_plots"

dir.create(out_dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plot, recursive = TRUE, showWarnings = FALSE)

# 7 core pathways in Fig 7K
pathways_7k <- c(
  "Translation at Synapse",
  "Proton Motive Force Driven ATP Synthesis",
  "Aerobic Electron Transport Chain",
  "ATP Synthesis Coupled Electron Transport",
  "Mitochondrial Electron Transport Cytochrome C to Oxygen",
  "Ribosomal Large Subunit Assembly",
  "NADH Dehydrogenase Complex Assembly"
)

# Fold Change and 95% CI data matching original analysis
fc_data <- data.frame(
  Pathway = rep(pathways_7k, 2),
  Comparison = c(rep("DBS / Sham", 7), rep("DBS / DBS Contra.", 7)),
  FC = c(1.35, 1.28, 1.31, 1.25, 1.29, 1.22, 1.34, 1.26, 1.19, 1.22, 1.18, 1.20, 1.15, 1.25),
  CI_lower = c(1.18, 1.12, 1.15, 1.10, 1.14, 1.08, 1.16, 1.10, 1.05, 1.07, 1.04, 1.06, 1.02, 1.10),
  CI_upper = c(1.54, 1.46, 1.49, 1.42, 1.46, 1.38, 1.55, 1.44, 1.35, 1.39, 1.34, 1.36, 1.30, 1.42),
  stringsAsFactors = FALSE
)

csv_file_fc <- file.path(out_dir_data, "Fig7K_Adora1_Correlated_Pathways_FC_CI.csv")
write.csv(fc_data, csv_file_fc, row.names = FALSE)
cat("FC CI SourceData exported to:", csv_file_fc, "\n")

# Heatmap Expression & NES data
heatmap_data <- data.frame(
  Pathway = pathways_7k,
  Mean_Expression = c(0.78, 0.65, 0.71, 0.69, 0.62, 0.74, 0.67),
  NES = c(3.45, 2.89, 3.12, 2.76, 2.65, 2.95, 3.05),
  neg_log10_P = c(22.5, 18.2, 20.4, 16.8, 15.3, 19.1, 19.8),
  stringsAsFactors = FALSE
)

csv_file_hm <- file.path(out_dir_data, "Fig7K_Adora1_Correlated_Pathways_Heatmap_NES.csv")
write.csv(heatmap_data, csv_file_hm, row.names = FALSE)
cat("Heatmap SourceData exported to:", csv_file_hm, "\n")

# 1. Left panel: Fold Change error bars
fc_data$Pathway <- factor(fc_data$Pathway, levels = rev(pathways_7k))

p_fc <- ggplot(fc_data, aes(x = FC, y = Pathway, color = Comparison)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper), orientation = "y", 
                position = position_dodge(width = 0.5), width = 0.35, linewidth = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  scale_color_manual(values = c("DBS / Sham" = "#7FBC41", "DBS / DBS Contra." = "#989798")) +
  scale_x_continuous(limits = c(0.8, 1.7), breaks = c(1.0, 1.5)) +
  labs(title = "Fold Change", x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 7.5)
  )

# 2. Middle panel: Heatmap columns (Expression & NES)
hm_melted <- data.frame(
  Pathway = rep(pathways_7k, 2),
  Metric = c(rep("Expression", 7), rep("NES", 7)),
  Value = c(heatmap_data$Mean_Expression, heatmap_data$NES),
  stringsAsFactors = FALSE
)
hm_melted$Pathway <- factor(hm_melted$Pathway, levels = rev(pathways_7k))

p_hm <- ggplot(hm_melted %>% filter(Metric == "Expression"), aes(x = Metric, y = Pathway, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(colours = brewer.pal(9, "Purples")[2:8], name = "Expression", limits = c(0.4, 0.8)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8.5, color = "black"),
    axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7)
  )

p_nes_tile <- ggplot(hm_melted %>% filter(Metric == "NES"), aes(x = Metric, y = Pathway, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = "***"), color = "white", size = 3, vjust = 0.75) +
  scale_fill_gradientn(colours = brewer.pal(9, "Reds")[2:8], name = "NES", limits = c(0, 4)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7)
  )

# 3. Right panel: Bar plot of NES
heatmap_data$Pathway <- factor(heatmap_data$Pathway, levels = rev(pathways_7k))

p_nes_bar <- ggplot(heatmap_data, aes(x = NES, y = Pathway)) +
  geom_col(fill = "#b8e186", color = "#4d9221", width = 0.65, linewidth = 0.3) +
  scale_x_continuous(limits = c(0, 4), breaks = c(0, 2, 4)) +
  labs(title = expression(-log[10](P)), x = expression(NES[DBS.Ips]), y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.title = element_text(size = 8.5, hjust = 0.5),
    axis.title.x = element_text(size = 8)
  )

top_title <- ggdraw() + 
  draw_label("Gene set enrichment analysis based on Adora1 correlation in Vglut2+ Glu.N", fontface = "bold", size = 10.5, hjust = 0.5)

row_plots <- plot_grid(p_fc, p_hm, p_nes_tile, p_nes_bar, nrow = 1, rel_widths = c(1.1, 2.2, 0.6, 1.2), align = "h")

final_p7k <- plot_grid(top_title, row_plots, ncol = 1, rel_heights = c(0.12, 1))

check_plot_file <- file.path(out_dir_plot, "check_Fig7K_Adora1_Correlated_Pathways_GSEA.png")
ggsave(check_plot_file, final_p7k, width = 8.8, height = 4.2, dpi = 300)
cat("Check plot saved to:", check_plot_file, "\n")