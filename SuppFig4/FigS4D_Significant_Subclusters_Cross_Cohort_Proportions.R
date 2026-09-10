# ==============================================================================
# Script Name: FigS4D_Significant_Subclusters_Cross_Cohort_Proportions.R
# Description: Stacked bar plots illustrating the proportional composition of
#              significant ADORA1-Low enriched subclusters across multiple
#              patients, brain regions, and epilepsy pathologies (Supplementary Figure 4D).
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)
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

sig_clusters <- c("6", "9", "17", "18", "20")
sig_df <- cell_meta %>%
  dplyr::filter(seurat_clusters %in% sig_clusters) %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = sig_clusters))

# Function to compute proportions and plot
calc_props <- function(data, group_var) {
  data %>%
    dplyr::group_by(seurat_clusters, !!sym(group_var)) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::mutate(proportion = count / sum(count)) %>%
    dplyr::ungroup()
}

dist_patient <- calc_props(sig_df, "Patient")
dist_region  <- calc_props(sig_df, "BrainRegion")
dist_disease <- calc_props(sig_df, "Disease")

# Export SourceData
cat("Exporting SourceData for Fig S4D...\n")
write.csv(
  dist_patient %>% dplyr::select(Subcluster = seurat_clusters, Patient, Cell_Count = count, Proportion = proportion),
  file = file.path(out_data_dir, "FigS4D_Significant_Subclusters_Patient_Proportions.csv"),
  row.names = FALSE
)

write.csv(
  dist_region %>% dplyr::select(Subcluster = seurat_clusters, BrainRegion, Cell_Count = count, Proportion = proportion),
  file = file.path(out_data_dir, "FigS4D_Significant_Subclusters_BrainRegion_Proportions.csv"),
  row.names = FALSE
)

write.csv(
  dist_disease %>% dplyr::select(Subcluster = seurat_clusters, Disease, Cell_Count = count, Proportion = proportion),
  file = file.path(out_data_dir, "FigS4D_Significant_Subclusters_Disease_Proportions.csv"),
  row.names = FALSE
)

# Plotting
cat("Generating Check Plot for Fig S4D...\n")
# Palettes: exactly matching 260205_revSTM.R line 1973
pal_tableau29 <- colorRampPalette(ggthemes::tableau_color_pal("Tableau 20")(20))(29)

p_patient <- ggplot(dist_patient, aes(x = seurat_clusters, y = proportion, fill = Patient)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = pal_tableau29) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_blank(),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm")
  ) +
  guides(fill = guide_legend(ncol = 2)) +
  labs(y = "proportion")

p_region <- ggplot(dist_region, aes(x = seurat_clusters, y = proportion, fill = BrainRegion)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = pal_tableau29) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_blank(),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm")
  ) +
  labs(y = "proportion", fill = "BrainRegion")

p_disease <- ggplot(dist_disease, aes(x = seurat_clusters, y = proportion, fill = Disease)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = pal_tableau29) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_blank(),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm")
  ) +
  labs(y = "proportion", fill = "Disease")

p_s4d <- p_patient + p_region + p_disease + plot_layout(nrow = 1)

ggsave(
  filename = file.path(out_check_dir, "check_FigS4D_Significant_Subclusters_Proportions.png"),
  plot = p_s4d,
  width = 11,
  height = 4,
  dpi = 300
)

cat("Fig S4D completed successfully.\n")
