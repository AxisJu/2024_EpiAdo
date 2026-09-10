# ==============================================================================
# Script Name: FigS5A_GluN_ADORA1_Groups_UMAP.R
# Description: UMAP of all glutamatergic neurons divided into Control,
#              ADORA1^High Epilepsy, and ADORA1^Low Epilepsy.
# Author: Antigravity Agent
# Date: 2026-09-09
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# Path definitions
script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) "./SuppFig5"
)
cache_path <- file.path(script_dir, "glun_s5_cache.rds")
out_csv <- "../sourcedata/SuppFig5/FigS5A_GluN_ADORA1_Groups_UMAP.csv"
check_plot_path <- "../sourcedata/check_plots/check_FigS5A_GluN_ADORA1_Groups_UMAP.png"

# 1. Load cached UMAP data
if (!file.exists(cache_path)) {
  stop("Cache file not found: ", cache_path)
}
cache <- readRDS(cache_path)
umap_df <- cache$umap_df

# Order dataframe so Control is at the bottom, then High, then Low on top
umap_df <- umap_df %>%
  arrange(factor(A1Group, levels = c("Control", "ADORA1^High Epilepsy", "ADORA1^Low Epilepsy")))

# 2. Export SourceData CSV
sourcedata_df <- umap_df %>%
  transmute(
    Cell_ID = Cell,
    UMAP_1 = round(UMAP_1, 4),
    UMAP_2 = round(UMAP_2, 4),
    Clinical_Group = Group,
    ADORA1_Expression_LogCPM = round(ADORA1_Exp, 4),
    ADORA1_Subgroup = as.character(A1Group)
  )
write.csv(sourcedata_df, out_csv, row.names = FALSE)
cat("SourceData saved to:", out_csv, "with", nrow(sourcedata_df), "cells.\n")

# 3. Plotting
color_palette <- c(
  "Control" = "#989798",
  "ADORA1^High Epilepsy" = "#205A9E",
  "ADORA1^Low Epilepsy" = "#7BC5B4"
)

# Custom small axes at bottom-left corner
x_min <- min(umap_df$UMAP_1)
y_min <- min(umap_df$UMAP_2)
axis_len <- 3

p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = A1Group)) +
  geom_point(size = 0.4, alpha = 0.85, stroke = 0) +
  scale_color_manual(values = color_palette, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
  # Draw custom L-shaped corner axis lines
  annotate("segment", x = x_min, xend = x_min + axis_len, y = y_min, yend = y_min,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 0.5) +
  annotate("segment", x = x_min, y = y_min, xend = x_min, yend = y_min + axis_len,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 0.5) +
  annotate("text", x = x_min + axis_len / 2, y = y_min - 0.5, label = "UMAP1", size = 3.2, family = "sans") +
  annotate("text", x = x_min - 0.5, y = y_min + axis_len / 2, label = "UMAP2", size = 3.2, angle = 90, family = "sans") +
  theme_void(base_size = 12) +
  theme(
    legend.position = c(0.25, 0.92),
    legend.text = element_text(size = 9, family = "sans"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_fixed(ratio = 1)

# Save check plot
ggsave(check_plot_path, plot = p, width = 5.5, height = 5.5, dpi = 300)
cat("Check plot saved to:", check_plot_path, "\n")
