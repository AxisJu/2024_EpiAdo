# ==============================================================================
# Script Name: FigS3A_StereoSeq_Sample_QC_Violin.R
# Description: Quality control metrics (nFeature_RNA and log10 nCount_RNA) for
#              human hippocampal Stereo-seq datasets at the individual sample
#              level (Supplementary Figure 3A).
# ==============================================================================

library(qs)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(GSEABase)

out_data_dir <- "../sourcedata/SuppFig3"
out_check_dir <- "../sourcedata/check_plots"

dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_check_dir, recursive = TRUE, showWarnings = FALSE)

cache_data <- file.path(out_data_dir, "stereoseq_s3_cache.rds")

sample_map <- c(
  "T1036" = "Ctrl.1",
  "T1035" = "Ctrl.2",
  "T940"  = "Epi.1",
  "T941"  = "Epi.2"
)

slice_titles <- c(
  "T1036" = "Slice: Control-1",
  "T1035" = "Slice: Control-2",
  "T940"  = "Slice: Epilepsy-1",
  "T941"  = "Slice: Epilepsy-2"
)

samples_order <- c("T1036", "T1035", "T940", "T941")
sample_labels <- c("Ctrl.1", "Ctrl.2", "Epi.1", "Epi.2")

celltype_levels <- c(
  "Glu.N", "GABA.N", "OPCs", "Oligo.", 
  "Mo/Mφ", "Microglia", "Astro.", "Endo.", "Epi.", "LowQuality"
)

if (!file.exists(cache_data)) {
  library(Seurat)
  cat("Loading GMT for Purine Metabolism...\n")
  gmt_path <- "Z:/2023_EpiAdo/bin/KEGG_metabolism_nc.gmt"
  purine_genes <- getGmt(gmt_path)[["Purine metabolism"]]@geneIds
  
  cat("Extracting Stereo-seq sample data...\n")
  data_list <- list()
  
  for (sm in samples_order) {
    cat("Loading sample:", sm, "...\n")
    obj_path <- paste0("Z:/2023_EpiAdo/data/StereoSeq_MSH/processed/", sm)
    seuratobj <- qread(obj_path)
    
    # Gene expression
    avail_purine <- intersect(purine_genes, rownames(seuratobj))
    purine_exp <- colMeans(as.matrix(GetAssayData(seuratobj, assay = "RNA", layer = "data")[avail_purine, ]))
    
    adora1_exp <- if ("ADORA1" %in% rownames(seuratobj)) {
      as.numeric(GetAssayData(seuratobj, assay = "RNA", layer = "data")["ADORA1", ])
    } else {
      rep(0, ncol(seuratobj))
    }
    
    df <- data.frame(
      SampleID = sm,
      Sample = sample_map[sm],
      Slice = slice_titles[sm],
      SpotID = colnames(seuratobj),
      X = seuratobj$x_center,
      Y = -seuratobj$y_center,
      CellSize = seuratobj$cell_size,
      Celltype = seuratobj$predicted.id,
      nCount_RNA = seuratobj$nCount_RNA,
      nFeature_RNA = seuratobj$nFeature_RNA,
      ADORA1 = adora1_exp,
      PurineMetabolism = purine_exp
    )
    
    data_list[[sm]] <- df
    rm(seuratobj); gc()
  }
  
  combined_data <- bind_rows(data_list)
  saveRDS(combined_data, file = cache_data)
  cat("Saved Stereo-seq cache successfully.\n")
} else {
  cat("Loading Stereo-seq cache...\n")
  combined_data <- readRDS(cache_data)
}

combined_data$Sample <- factor(combined_data$Sample, levels = sample_labels)
combined_data$Celltype <- factor(combined_data$Celltype, levels = celltype_levels)

# Subsample for SourceData CSV (up to 5,000 spots per sample)
set.seed(42)
qc_sub <- combined_data %>%
  group_by(Sample) %>%
  sample_n(min(n(), 5000)) %>%
  ungroup() %>%
  dplyr::select(Sample, SampleID, SpotID, Celltype, nCount_RNA, nFeature_RNA)

qc_stats <- combined_data %>%
  group_by(Sample) %>%
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
  file = file.path(out_data_dir, "FigS3A_StereoSeq_Sample_QC_Metrics.csv"),
  row.names = FALSE
)

write.csv(
  qc_stats,
  file = file.path(out_data_dir, "FigS3A_StereoSeq_Sample_QC_Statistics.csv"),
  row.names = FALSE
)

cat("Plotting Check Figures...\n")
p_vln_feature_sample <- ggplot(combined_data, aes(x = Sample, y = nFeature_RNA, fill = Sample)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.4) +
  scale_fill_viridis_d(option = "turbo") +
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

p_vln_count_sample <- ggplot(combined_data, aes(x = Sample, y = log10(nCount_RNA), fill = Sample)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.4) +
  scale_fill_viridis_d(option = "turbo") +
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

p_s3a <- (p_vln_feature_sample | p_vln_count_sample) +
  plot_layout(widths = c(1, 1))

ggsave(
  filename = file.path(out_check_dir, "check_FigS3A_Sample_QC_Violin.png"),
  plot = p_s3a,
  width = 6.5,
  height = 4.5,
  dpi = 300
)

cat("FigS3A completed successfully.\n")
