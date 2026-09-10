# ==============================================================================
# Script: Fig1F_Purine_Ado_Correlation_Across_CellTypes.R
# Purpose: Correlation between purine metabolism and adenosine metabolic genes across cell types
# Article: EpiAdo (Science Translational Medicine) - Figure 1F
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(ggplot2)
})

out_dir_sd <- "../sourcedata/Figure1"
out_dir_plot <- "../sourcedata/check_plots"
dir.create(out_dir_sd, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_plot, showWarnings = FALSE, recursive = TRUE)

cat("Loading precalculated correlation data for Figure 1F...\n")
# [DATA REQUIRED] The following data load has been removed (source: cor.mat_metabonly_allcelltypes_purine).
# Please load the required data object before running this script.

target_levels <- c("Glu.N", "GABA.N", "OPCs", "Oligo.", "Microglia", "Astro.")
plot_data_points <- cor.mat.metabonly %>% 
  dplyr::filter(Celltype %in% target_levels) %>%
  dplyr::mutate(Celltype = factor(Celltype, levels = target_levels))

plot_data_summary <- plot_data_points %>%
  dplyr::group_by(Celltype) %>%
  dplyr::summarise(
    TotalGenes = dplyr::n(),
    SigNum_P01 = sum(P < 0.01),
    SigPercent = sum(P < 0.01) / dplyr::n(),
    Mean_Rho = mean(Cor),
    Median_Rho = median(Cor),
    SD_Rho = sd(Cor),
    SEM_Rho = sd(Cor) / sqrt(dplyr::n()),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Celltype = factor(Celltype, levels = target_levels))

# Export SourceData
sd_summary_file <- file.path(out_dir_sd, "Fig1F_Purine_Ado_Correlation_Summary.csv")
sd_pergene_file <- file.path(out_dir_sd, "Fig1F_Purine_Ado_PerGene_Correlation.csv")

write.csv(plot_data_summary, sd_summary_file, row.names = FALSE)
write.csv(plot_data_points, sd_pergene_file, row.names = FALSE)
cat("SourceData saved for Figure 1F.\n")

# Check Plot
color_palette <- c(
  "Glu.N" = "#d0a9cb", "GABA.N" = "#9dcb84", "OPCs" = "#cea27b", 
  "Oligo." = "#f3d490", "Microglia" = "#dc8384", "Astro." = "#adaed3"
)

p_combined <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_col(
    data = plot_data_summary,
    aes(x = Celltype, y = Median_Rho, color = Celltype, fill = Celltype),
    alpha = 0.15, width = 0.6, linewidth = 0.8
  ) +
  geom_jitter(
    data = plot_data_points,
    aes(x = Celltype, y = Cor, color = Celltype),
    alpha = 0.4, size = 2, width = 0.08, height = 0
  ) +
  geom_line(
    data = plot_data_summary,
    aes(x = Celltype, y = SigPercent - 0.55, group = 1),
    color = "#d0a9cb", linetype = "dashed", linewidth = 1, alpha = 0.8
  ) +
  geom_point(
    data = plot_data_summary,
    aes(x = Celltype, y = SigPercent - 0.55, color = Celltype, fill = Celltype),
    size = 6, alpha = 1
  ) +
  scale_y_continuous(
    name = "Magnitude of Correlation (rho)",
    breaks = seq(-0.1, 0.2, by = 0.1),
    sec.axis = sec_axis(
      trans = ~ (. + 0.55) * 100,
      name = "Extent of Correlation (%)",
      breaks = seq(40, 80, by = 10)
    )
  ) +
  coord_cartesian(ylim = c(-0.15, 0.25)) +
  scale_fill_manual(values = color_palette) +
  scale_colour_manual(values = color_palette) +
  theme_classic() +
  theme(
    plot.margin = unit(c(0.6, 0.6, 0.6, 0.6), units = "cm"),
    plot.title = element_text(colour = "black", size = 13, family = "sans", face = "bold", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y.left = element_text(colour = "black", size = 12, family = "sans"),
    axis.title.y.right = element_text(colour = "#8e7a93", size = 12, family = "sans", angle = 90, vjust = 1),
    axis.text.x = element_text(colour = "black", size = 11, family = "sans"),
    axis.text.y.left = element_text(colour = "black", size = 11, family = "sans"),
    axis.text.y.right = element_text(colour = "#8e7a93", size = 11, family = "sans"),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = unit(2, units = "mm"),
    legend.position = "none"
  ) +
  labs(title = "Correlation between purine and adenosine metabolism")

plot_out_file <- file.path(out_dir_plot, "check_Fig1F_Combined_LineBar.png")
ggsave(plot_out_file, plot = p_combined, width = 7.5, height = 5, dpi = 300)
cat("Figure 1F finished successfully.\n")
