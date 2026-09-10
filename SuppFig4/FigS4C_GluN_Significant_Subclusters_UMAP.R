# ==============================================================================
# Script Name: FigS4C_GluN_Significant_Subclusters_UMAP.R
# Description: UMAP highlighting the specific glutamatergic subclusters that are
#              significantly enriched for the ADORA1-Low state (Supplementary Figure 4C).
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggthemes)

out_data_dir <- "../sourcedata/SuppFig4"
out_check_dir <- "../sourcedata/check_plots"

cache_path <- file.path(out_data_dir, "glun_s4_cache.rds")
if (!file.exists(cache_path)) {
  stop("glun_s4_cache.rds not found. Please ensure cache is built.")
}

cat("Loading SuppFig4 cache...\n")
s4_cache <- readRDS(cache_path)
cell_meta <- s4_cache$cell_meta

sig_clusters <- c("9", "20", "17", "18", "6")
all_clusters <- as.character(sort(as.numeric(unique(cell_meta$seurat_clusters))))

# Palette: Tableau 20 colors for the 5 significant clusters, grey90 for the rest
sig_palette <- as.character(ggthemes::tableau_color_pal("Tableau 20")(length(sig_clusters)))
cluster_colors <- setNames(rep("grey90", length(all_clusters)), all_clusters)
cluster_colors[sig_clusters] <- sig_palette

# Export SourceData
cat("Exporting SourceData for Fig S4C...\n")
write.csv(
  cell_meta %>%
    dplyr::mutate(Is_Significant_Enriched = seurat_clusters %in% sig_clusters) %>%
    dplyr::select(
      Cell_Barcode = Cell,
      UMAP_1,
      UMAP_2,
      Subcluster = seurat_clusters,
      Is_Significant_Enriched
    ),
  file = file.path(out_data_dir, "FigS4C_GluN_Significant_Subclusters_UMAP.csv"),
  row.names = FALSE
)

# Calculate label centers for all clusters
centers <- cell_meta %>%
  group_by(seurat_clusters) %>%
  summarize(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

# Separate non-sig and sig points for layering
df_nonsig <- cell_meta %>% filter(!seurat_clusters %in% sig_clusters)
df_sig    <- cell_meta %>% filter(seurat_clusters %in% sig_clusters)

cat("Generating Check Plot for Fig S4C...\n")
p_s4c <- ggplot() +
  theme_void() +
  rasterise(
    geom_point(data = df_nonsig, aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5),
    dpi = 300
  ) +
  rasterise(
    geom_point(data = df_sig, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters), size = 0.6),
    dpi = 300
  ) +
  geom_text(
    data = centers,
    aes(x = UMAP_1, y = UMAP_2, label = seurat_clusters),
    size = 4,
    color = "black"
  ) +
  scale_color_manual(
    values = cluster_colors,
    breaks = all_clusters,
    name = "Significant ADORA1Low\nEnriched Clusters"
  ) +
  guides(
    color = guide_legend(
      ncol = 2,
      override.aes = list(size = 3, shape = 16)
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  ) +
  labs(title = "Significant ADORA1-Low Enriched Clusters")

ggsave(
  filename = file.path(out_check_dir, "check_FigS4C_GluN_Significant_Subclusters_UMAP.png"),
  plot = p_s4c,
  width = 9,
  height = 7,
  dpi = 300
)

cat("Fig S4C completed successfully.\n")
