# ==============================================================================
# Script Name: FigS4A_GluN_ADORA1_Groups_UMAP.R
# Description: UMAP visualization of glutamatergic neurons color-coded by defined
#              ADORA1 expression groups (Control, Epi. ADORA1-High, Epi. ADORA1-Low)
#              (Supplementary Figure 4A).
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggrastr)

out_data_dir <- "../sourcedata/SuppFig4"
out_check_dir <- "../sourcedata/check_plots"

cache_path <- file.path(out_data_dir, "glun_s4_cache.rds")
if (!file.exists(cache_path)) {
  stop("glun_s4_cache.rds not found. Please ensure cache is built.")
}

cat("Loading SuppFig4 cache...\n")
s4_cache <- readRDS(cache_path)
cell_meta <- s4_cache$cell_meta

# Color palette from original publication
cols_grp <- c(
  "Control" = "#989798",
  "Epi_ADORA1_High" = "#205A9E",
  "Epi_ADORA1_Low" = "#7BC5B4"
)

# Export SourceData
cat("Exporting SourceData for Fig S4A...\n")
sourcedata_a <- cell_meta %>%
  dplyr::select(
    Cell_Barcode = Cell,
    UMAP_1,
    UMAP_2,
    ADORA1_Group = Compare_Group,
    ADORA1_CPM = ADORA1_Val
  )

write.csv(
  sourcedata_a,
  file = file.path(out_data_dir, "FigS4A_GluN_ADORA1_Groups_UMAP.csv"),
  row.names = FALSE
)

cat("Generating Check Plot for Fig S4A...\n")
# Shuffle points to avoid layer bias
set.seed(42)
cell_meta_shuffled <- cell_meta[sample(nrow(cell_meta)), ]

p_s4a <- ggplot(cell_meta_shuffled, aes(x = UMAP_1, y = UMAP_2, color = Compare_Group)) +
  rasterise(geom_point(size = 0.5, alpha = 0.8), dpi = 300) +
  scale_color_manual(
    values = cols_grp,
    labels = c("Control" = "Control", "Epi_ADORA1_High" = "Epi. ADORA1^High", "Epi_ADORA1_Low" = "Epi. ADORA1^Low"),
    name = "ADORA1 Groups"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  ) +
  labs(title = "ADORA1 Groups")

ggsave(
  filename = file.path(out_check_dir, "check_FigS4A_GluN_ADORA1_Groups_UMAP.png"),
  plot = p_s4a,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Fig S4A completed successfully.\n")
