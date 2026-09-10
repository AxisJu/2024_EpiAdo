# ==============================================================================
# Script: FigS7A_PatchSeq_Metadata_UMAP.R
# Description: UMAP visualization of clinical and anatomical metadata 
#              (Age, Gender, Brain Region, Layer) for patched cells.
# Panel: Supplementary Figure 7A
# Source Data:
#   - FigS7A_PatchSeq_Metadata_UMAP.csv
# Output:
#   - ../../sourcedata-260907/check_plots/check_FigS7A_PatchSeq_Metadata_UMAP.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
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

# ==============================================================================
# Export Source Data
# ==============================================================================
s7a_export <- s7_df %>%
  dplyr::select(Cell_ID, Sample, Patient, UMAP_1, UMAP_2, Age, Gender, BrainRegion, Layer)

write.csv(s7a_export, file.path(data_dir, "FigS7A_PatchSeq_Metadata_UMAP.csv"), row.names = FALSE)
cat("Exported SourceData to:", file.path(data_dir, "FigS7A_PatchSeq_Metadata_UMAP.csv"), "\n")

# ==============================================================================
# Plot generation
# ==============================================================================

# Base UMAP theme settings
theme_umap_base <- theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13, margin = margin(b = 6)),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.45, "cm"),
    legend.margin = margin(l = 2, r = 8)
  )

# 1. Age
p_age <- ggplot(s7_df, aes(x = UMAP_1, y = UMAP_2, color = Age)) +
  geom_point(size = 2.0, alpha = 0.95) +
  scale_color_gradientn(
    colours = colorRampPalette(c("#fce2d8", "#c89d8d"))(31),
    breaks = c(20, 40, 60),
    limits = c(min(s7_df$Age, na.rm = TRUE), max(s7_df$Age, na.rm = TRUE))
  ) +
  labs(title = "Age", color = NULL) +
  theme_umap_base

# 2. Gender
cols_gender <- c("Female" = "#fac5b1", "Male" = "#9fa1cb")
p_gender <- ggplot(s7_df, aes(x = UMAP_1, y = UMAP_2, color = Gender)) +
  geom_point(size = 2.0, alpha = 0.95) +
  scale_color_manual(values = cols_gender) +
  labs(title = "Gender", color = NULL) +
  theme_umap_base +
  guides(color = guide_legend(override.aes = list(size = 3.0)))

# 3. Brain Region
cols_region <- c(
  "Frontal" = "#D6A897",
  "Parietal" = "#FAC5B1",
  "Temporal" = "#FACAB8",
  "Ocipital" = "#FBD4C6",
  "Hippocampus" = "#9fa1cb"
)
p_region <- ggplot(s7_df, aes(x = UMAP_1, y = UMAP_2, color = BrainRegion)) +
  geom_point(size = 2.0, alpha = 0.95) +
  scale_color_manual(values = cols_region) +
  labs(title = "Brain Region", color = NULL) +
  theme_umap_base +
  guides(color = guide_legend(override.aes = list(size = 3.0)))

# 4. Layers
cols_layer <- c(
  "L1" = "#a3dac8",
  "L2" = "#81bc61",
  "L2/3" = "#9dcb84",
  "L3" = "#cea27b",
  "L4" = "#f3d490",
  "L5" = "#bd86b6",
  "L5/6" = "#d0a9cb",
  "DG" = "#9b9cbd",
  "CA" = "#c5d3e0",
  "NOS" = "#808080"
)
p_layer <- ggplot(s7_df, aes(x = UMAP_1, y = UMAP_2, color = Layer)) +
  geom_point(size = 2.0, alpha = 0.95) +
  scale_color_manual(values = cols_layer) +
  labs(title = "Layers", color = NULL) +
  theme_umap_base +
  guides(color = guide_legend(override.aes = list(size = 3.0), ncol = 1))

# Add right-angle axes (UMAP1, UMAP2) to the first panel
arrow_style <- arrow(length = unit(0.20, "cm"), type = "closed")
axis_df <- data.frame(
  x = c(-12.5, -12.5),
  y = c(-8.5, -8.5),
  xend = c(-6.5, -12.5),
  yend = c(-8.5, -2.5)
)

p_age_axis <- p_age +
  geom_segment(data = axis_df[1,], aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow_style, color = "black", linewidth = 0.65, inherit.aes = FALSE) +
  geom_segment(data = axis_df[2,], aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow_style, color = "black", linewidth = 0.65, inherit.aes = FALSE) +
  annotate("text", x = -9.5, y = -9.8, label = "UMAP1", size = 3.2, fontface = "bold") +
  annotate("text", x = -13.8, y = -5.5, label = "UMAP2", size = 3.2, fontface = "bold", angle = 90) +
  coord_cartesian(clip = "off")

# Assemble Panel A
p_s7a <- p_age_axis + p_gender + p_region + p_layer +
  plot_layout(nrow = 1, widths = c(1, 1, 1, 1))

check_png <- file.path(plot_dir, "check_FigS7A_PatchSeq_Metadata_UMAP.png")
ggsave(check_png, plot = p_s7a, width = 15, height = 4.2, dpi = 300)
cat("Check plot saved to:", check_png, "\n")
