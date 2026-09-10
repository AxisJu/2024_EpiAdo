# ==============================================================================
# Script: FigS7B_PatchSeq_EFeatures_UMAP.R
# Description: UMAP projection of 12 electrophysiological features (E-Features)
#              for patched human neurons.
# Panel: Supplementary Figure 7B
# Source Data:
#   - FigS7B_PatchSeq_EFeatures_UMAP.csv
# Output:
#   - ../../sourcedata-260907/check_plots/check_FigS7B_PatchSeq_EFeatures_UMAP.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(grid)
})

base_dir <- here::here()  # Project root; use here::here() or setwd() as appropriate
code_dir <- "."  # SuppFig7
data_dir <- file.path("..", "sourcedata", "SuppFig7")
plot_dir <- file.path("..", "sourcedata", "check_plots")

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

cache_path <- file.path(code_dir, "ps_s7_cache.rds")
s7_df <- readRDS(cache_path)

# 12 electrophysiological features (matching published figure order and labels)
feature_list <- list(
  list(raw = "MembraneCapacitance", scaled = "MembraneCapacitance_Scaled", title = "Membrane Capacitance"),
  list(raw = "RestingMembranePotential", scaled = "RestingMembranePotential_Scaled", title = "Resting Membrane Potential"),
  list(raw = "FirstSpikeLatency", scaled = "FirstSpikeLatency_Scaled", title = "First Spike Latency"),
  list(raw = "AP_UpStroke", scaled = "AP_UpStroke_Scaled", title = "AP UpStroke"),
  list(raw = "AP_DownStroke", scaled = "AP_DownStroke_Scaled", title = "AP DownStroke"),
  list(raw = "AP_HalfWidth", scaled = "AP_HalfWidth_Scaled", title = "AP Half-Width"),
  list(raw = "AP_Threshold", scaled = "AP_Threshold_Scaled", title = "AP Threshold"),
  list(raw = "AP_Duration", scaled = "AP_Duration_Scaled", title = "AP Duration"),
  list(raw = "AP_Amplitude", scaled = "AP_Amplitude_Scaled", title = "AP Amplitude"),
  list(raw = "AfterHyperpolarizationAmplitude", scaled = "AfterHyperpolarizationAmplitude_Scaled", title = "After Hyperpolarization Amplitude"),
  list(raw = "SagDepth", scaled = "SagDepth_Scaled", title = "Sag Depth"),
  list(raw = "Rheobase", scaled = "Rheobase_Scaled", title = "Rheobase")
)

# ==============================================================================
# Export Source Data
# ==============================================================================
raw_cols <- sapply(feature_list, function(x) x$raw)
scaled_cols <- sapply(feature_list, function(x) x$scaled)

s7b_export <- s7_df %>%
  dplyr::select(Cell_ID, Sample, Patient, UMAP_1, UMAP_2, all_of(raw_cols), all_of(scaled_cols))

write.csv(s7b_export, file.path(data_dir, "FigS7B_PatchSeq_EFeatures_UMAP.csv"), row.names = FALSE)
cat("Exported SourceData to:", file.path(data_dir, "FigS7B_PatchSeq_EFeatures_UMAP.csv"), "\n")

# ==============================================================================
# Plot generation
# ==============================================================================

# Color gradient: light green -> white -> magenta (PiYG 3:8 reversed)
cols_palette <- colorRampPalette(rev(brewer.pal(11, "PiYG")[3:8]))(100)

plot_list <- list()

for (i in seq_along(feature_list)) {
  item <- feature_list[[i]]
  f_raw <- item$raw
  f_scaled <- item$scaled
  f_title <- item$title
  
  df_sub <- s7_df %>%
    dplyr::select(UMAP_1, UMAP_2, Scaled = all_of(f_scaled)) %>%
    arrange(!is.na(Scaled), Scaled)  # NA values plotted on bottom; valid values on top
  
  df_na <- df_sub %>% filter(is.na(Scaled))
  df_val <- df_sub %>% filter(!is.na(Scaled))
  
  p <- ggplot() +
    # NA missing values: grey
    geom_point(data = df_na, aes(x = UMAP_1, y = UMAP_2), color = "grey60", size = 1.4, alpha = 0.8) +
    # Valid data points: continuous color gradient
    geom_point(data = df_val, aes(x = UMAP_1, y = UMAP_2, color = Scaled), size = 1.8, alpha = 0.95) +
    scale_color_gradientn(
      colours = cols_palette,
      limits = c(1.0, 2.0),
      breaks = c(1.0, 2.0),
      labels = c("1.00", "2.00")
    ) +
    labs(title = f_title) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9.5, margin = margin(b = 4)),
      legend.position = "none"
    )
  
  # Add right-angle axes to the lower-left of panel at row 2, col 1 (AP Threshold, i = 7)
  if (i == 7) {
    arrow_style <- arrow(length = unit(0.18, "cm"), type = "closed")
    axis_df <- data.frame(
      x = c(-12.5, -12.5),
      y = c(-8.5, -8.5),
      xend = c(-6.5, -12.5),
      yend = c(-8.5, -2.5)
    )
    p <- p +
      geom_segment(data = axis_df[1,], aes(x = x, y = y, xend = xend, yend = yend),
                   arrow = arrow_style, color = "black", linewidth = 0.6, inherit.aes = FALSE) +
      geom_segment(data = axis_df[2,], aes(x = x, y = y, xend = xend, yend = yend),
                   arrow = arrow_style, color = "black", linewidth = 0.6, inherit.aes = FALSE) +
      annotate("text", x = -9.5, y = -10.2, label = "UMAP1", size = 2.8, fontface = "bold") +
      annotate("text", x = -14.2, y = -5.5, label = "UMAP2", size = 2.8, fontface = "bold", angle = 90) +
      coord_cartesian(clip = "off")
  }
  
  plot_list[[i]] <- p
}

# Construct standalone right-side legend
dummy_df <- data.frame(x = 1, y = 1, Scaled = c(1.0, 2.0))
p_legend_dummy <- ggplot(dummy_df, aes(x = x, y = y, color = Scaled)) +
  geom_point(alpha = 0) +
  scale_color_gradientn(
    colours = cols_palette,
    limits = c(1.0, 2.0),
    breaks = c(1.0, 2.0),
    labels = c("2.00", "1.00"),
    name = "Scaled\nE-Features"
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9.5, face = "bold", hjust = 0),
    legend.text = element_text(size = 9, face = "bold"),
    legend.key.height = unit(1.8, "cm"),
    legend.key.width = unit(0.35, "cm")
  )

legend_grob <- cowplot::get_legend(p_legend_dummy)

# Assemble 2-row x 6-column grid layout
grid_12 <- wrap_plots(plot_list, nrow = 2, ncol = 6)

p_final <- cowplot::plot_grid(
  grid_12, p_legend_dummy,
  nrow = 1, rel_widths = c(20, 1.2)
)

check_png <- file.path(plot_dir, "check_FigS7B_PatchSeq_EFeatures_UMAP.png")
ggsave(check_png, plot = p_final, width = 19.5, height = 6.8, dpi = 300)
cat("Check plot saved to:", check_png, "\n")
