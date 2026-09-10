# ==============================================================================
# Script Name: FigS3B_StereoSeq_Celltype_QC_Violin.R
# Description: Quality control metrics (nFeature_RNA and log10 nCount_RNA) for
#              human hippocampal Stereo-seq datasets at the annotated cell-type
#              level (Supplementary Figure 3B).
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

out_data_dir <- "../sourcedata/SuppFig3"
out_check_dir <- "../sourcedata/check_plots"

cache_data <- file.path(out_data_dir, "stereoseq_s3_cache.rds")

if (!file.exists(cache_data)) {
  stop("Please run FigS3A first to generate stereoseq_s3_cache.rds cache.")
}

cat("Loading Stereo-seq cache...\n")
combined_data <- readRDS(cache_data)

celltype_levels <- c(
  "Glu.N", "GABA.N", "OPCs", "Oligo.", 
  "Mo/Mφ", "Microglia", "Astro.", "Endo.", "Epi."
)

combined_data <- combined_data %>%
  dplyr::filter(Celltype %in% celltype_levels)

combined_data$Celltype <- factor(combined_data$Celltype, levels = celltype_levels)

# Author palette from 260205_revSTM.R
color_res.coarse <- c(
  "Glu.N" = "#d0a9cb",
  "GABA.N" = "#9dcb84",
  "OPCs" = "#cea27b",
  "Oligo." = "#f3d490",
  "Mo/Mφ" = "#adaed3",
  "Microglia" = "#adc1d3",
  "Astro." = "#d25f61",
  "Endo." = "#dc8384",
  "Epi." = "#c5d3e0"
)

# Subsample for SourceData CSV
set.seed(42)
qc_sub <- combined_data %>%
  group_by(Celltype) %>%
  sample_n(min(n(), 3000)) %>%
  ungroup() %>%
  dplyr::select(Celltype, Sample, SampleID, SpotID, nCount_RNA, nFeature_RNA)

qc_stats <- combined_data %>%
  group_by(Celltype) %>%
  summarise(
    Spot_Count = n(),
    nFeature_Median = median(nFeature_RNA),
    nFeature_IQR = IQR(nFeature_RNA),
    nCount_Median = median(nCount_RNA),
    nCount_IQR = IQR(nCount_RNA),
    log10_nCount_Median = median(log10(nCount_RNA)),
    .groups = "drop"
  )

# Export SourceData
cat("Exporting SourceData...\n")
write.csv(
  qc_sub,
  file = file.path(out_data_dir, "FigS3B_StereoSeq_Celltype_QC_Metrics.csv"),
  row.names = FALSE
)

write.csv(
  qc_stats,
  file = file.path(out_data_dir, "FigS3B_StereoSeq_Celltype_QC_Statistics.csv"),
  row.names = FALSE
)

cat("Plotting Check Figures...\n")
p_vln_feature_cell <- ggplot(combined_data, aes(x = Celltype, y = nFeature_RNA, fill = Celltype)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.4) +
  scale_fill_manual(values = color_res.coarse) +
  theme_classic() +
  labs(y = "nFeature_RNA", x = "") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 11, face = "bold"),
    plot.margin = unit(c(0.5, 0.3, 0.5, 0.5), "cm")
  ) +
  ylim(0, 5000)

p_vln_count_cell <- ggplot(combined_data, aes(x = Celltype, y = log10(nCount_RNA), fill = Celltype)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.4) +
  scale_fill_manual(values = color_res.coarse) +
  theme_classic() +
  labs(y = "log10 nCount_RNA", x = "") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 11, face = "bold"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.3), "cm")
  ) +
  ylim(0, 4.5)

p_s3b <- (p_vln_feature_cell | p_vln_count_cell) +
  plot_layout(widths = c(1, 1))

ggsave(
  filename = file.path(out_check_dir, "check_FigS3B_Celltype_QC_Violin.png"),
  plot = p_s3b,
  width = 8.5,
  height = 4.5,
  dpi = 300
)

cat("FigS3B completed successfully.\n")
