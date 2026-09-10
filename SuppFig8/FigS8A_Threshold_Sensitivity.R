# ==============================================================================
# Script: FigS8A_Threshold_Sensitivity.R
# Description: Sensitivity analysis of Cohen's d for 4 key electrophysiological
#              properties across varying ADORA1 expression thresholds.
# Panel: Supplementary Figure 8A
# Source Data:
#   - FigS8A_Threshold_Sensitivity_Statistics.csv
# Output:
#   - ../../sourcedata-260907/check_plots/check_FigS8A_Threshold_Sensitivity.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- here::here()  # Project root; use here::here() or setwd() as appropriate
code_dir <- "."  # SuppFig8
data_dir <- file.path("..", "sourcedata", "SuppFig8")
plot_dir <- file.path("..", "sourcedata", "check_plots")

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

cache_path <- file.path(code_dir, "ps_s8_cache.rds")
s8_cache <- readRDS(cache_path)
df_thresh <- s8_cache$threshold_res

# ==============================================================================
# Export Source Data
# ==============================================================================
s8a_export <- df_thresh %>%
  dplyr::select(Property, Step, StepNum, Threshold, CohenD, Pval) %>%
  dplyr::mutate(
    MinusLog10P = -log10(Pval),
    Significance = ifelse(Pval < 0.05, "p < 0.05", "n.s.")
  )

write.csv(s8a_export, file.path(data_dir, "FigS8A_Threshold_Sensitivity_Statistics.csv"), row.names = FALSE)
cat("SourceData exported to:", file.path(data_dir, "FigS8A_Threshold_Sensitivity_Statistics.csv"), "\n")

# ==============================================================================
# Plot generation
# ==============================================================================
features <- c("AP_UpStroke", "AP_UpDownRatio", "AP_HalfWidth", "AP_Duration")
plot_list <- list()

for (i in seq_along(features)) {
  prop <- features[i]
  sub_df <- s8a_export %>% filter(Property == prop)
  
  p <- ggplot(sub_df, aes(x = StepNum, y = CohenD)) +
    # Background horizontal reference line y=0
    geom_hline(yintercept = 0, linetype = "solid", color = "grey80", linewidth = 0.5) +
    # Vertical reference line x=0 (baseline Mean)
    geom_vline(xintercept = 0, linetype = "dashed", color = "#d73027", alpha = 0.5) +
    # Trend line
    geom_line(color = "#2c7fb8", linewidth = 1, alpha = 0.4) +
    # Scatter points
    geom_point(aes(size = MinusLog10P, fill = Significance),
               shape = 21, color = "white", stroke = 0.5) +
    scale_fill_manual(values = c("p < 0.05" = "#2c7fb8", "n.s." = "grey70")) +
    scale_size_continuous(range = c(2, 5), limits = c(0.1, 1.7), breaks = c(0.4, 0.8, 1.2, 1.6)) +
    scale_x_continuous(
      breaks = seq(-0.1, 0.1, 0.1),
      labels = c("-0.1", "Mean", "+0.1"),
      limits = c(-0.16, 0.16)
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(color = "black", size = 9),
      panel.grid.major.y = element_line(color = "grey95"),
      legend.position = if (i == 4) "right" else "none"
    ) +
    labs(
      title = prop,
      x = if (i == 2) "ADORA1 Expression Threshold Shift" else NULL,
      y = if (i == 1 || i == 3) "Effect Size (Cohen's d)" else "Effect Size (Cohen's d)",
      fill = "Significance",
      size = "-log10(p)"
    )
  
  plot_list[[i]] <- p
}

p_s8a <- wrap_plots(plot_list, nrow = 1) +
  plot_layout(guides = "collect")

check_png <- file.path(plot_dir, "check_FigS8A_Threshold_Sensitivity.png")
ggsave(check_png, plot = p_s8a, width = 16, height = 4.2, dpi = 300)
cat("Check plot saved to:", check_png, "\n")
