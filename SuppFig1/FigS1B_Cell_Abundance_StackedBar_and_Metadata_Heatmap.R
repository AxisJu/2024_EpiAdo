# ==============================================================================
# Script Name: FigS1B_Cell_Abundance_StackedBar_and_Metadata_Heatmap.R
# Description: Stacked bar plots of main cell types across samples with average
#              abundance and clinical metadata heatmap (Supplementary Figure 1B).
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)

out_data_dir <- "../sourcedata/SuppFig1"
out_check_dir <- "../sourcedata/check_plots"

source_csv <- file.path(out_data_dir, "FigS1B_Sample_Cell_Fraction_and_Metadata.csv")
avg_csv <- file.path(out_data_dir, "FigS1B_Group_Average_Cell_Fraction.csv")

if (!file.exists(source_csv) || !file.exists(avg_csv)) {
  library(Seurat)
  library(qs)
  obj_path <- "Z:/2023_EpiAdo/data/scSeq_MSH/Merged/merged_ps_log_pca_harmony_umap_final_seuratobj"
  cat("Loading Seurat metadata...\n")
  merged_seuratobj <- qread(obj_path)
  meta <- merged_seuratobj@meta.data
  
  sample_meta <- meta %>%
    dplyr::select(
      Sample = Sample_snSeq.processed,
      Patient, Age, Gender, Group, Disease, BrainRegion,
      OTS_Time, FocalSeizure = FS_Freq, GeneralizedSeizure = GS_Freq,
      GS_Burden.calc, FS_Burden.calc
    ) %>%
    distinct(Sample, .keep_all = TRUE)
  
  sample_meta$TotalSeizureBurden <- log(sample_meta$GS_Burden.calc + sample_meta$FS_Burden.calc + 1)
  
  cell_counts <- meta %>%
    group_by(Sample = Sample_snSeq.processed, Celltype = celltype_coarse) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    ungroup()
  
  microglia_frac <- cell_counts %>%
    dplyr::filter(Celltype == "Microglia") %>%
    dplyr::select(Sample, MicrogliaFrac = Fraction)
  
  sample_meta <- left_join(sample_meta, microglia_frac, by = "Sample")
  sample_meta$MicrogliaFrac[is.na(sample_meta$MicrogliaFrac)] <- 0
  
  ctrl_samples <- sample_meta %>% dplyr::filter(Group == "Control") %>% arrange(MicrogliaFrac) %>% pull(Sample)
  epi_samples <- sample_meta %>% dplyr::filter(Group == "Epilepsy") %>% arrange(MicrogliaFrac) %>% pull(Sample)
  ordered_samples <- c(ctrl_samples, epi_samples)
  
  sample_meta$Sample <- factor(sample_meta$Sample, levels = ordered_samples)
  sample_meta <- sample_meta %>% arrange(Sample)
  
  cell_frac_wide <- cell_counts %>%
    dplyr::select(Sample, Celltype, Fraction) %>%
    pivot_wider(names_from = Celltype, values_from = Fraction, values_fill = 0)
  
  sample_full_df <- left_join(sample_meta, cell_frac_wide, by = "Sample")
  
  write.csv(sample_full_df, file = source_csv, row.names = FALSE)
  
  group_avg <- meta %>%
    group_by(Group, Celltype = celltype_coarse) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(Group) %>%
    mutate(AverageFraction = Count / sum(Count)) %>%
    ungroup()
  
  write.csv(group_avg, file = avg_csv, row.names = FALSE)
}

cat("Reading SourceData...\n")
sample_full_df <- read.csv(source_csv)
group_avg <- read.csv(avg_csv)

# Palette for Cell Types (matches paper)
color_res_coarse <- c(
  "#d88d9e", "#c3c88c", "#cea27b", "#f3d490",
  "#d25f61", "#dc8384", "#adaed3", "#adc1d3", "#c5d3e0",
  "#EEEEEE"
)
celltypes <- c(
  "Glu.N", "GABA.N", "OPCs", "Oligo.",
  "Mo/Mφ", "Microglia", "Astro.", "Endo.", "Epi.", "LowQuality"
)
names(color_res_coarse) <- celltypes

# Melt cell fraction
cell_cols <- intersect(celltypes, colnames(sample_full_df))
cell_frac_long <- sample_full_df %>%
  dplyr::select(Sample, Group, all_of(cell_cols)) %>%
  pivot_longer(cols = all_of(cell_cols), names_to = "Celltype", values_to = "Fraction")

cell_frac_long$Celltype <- factor(cell_frac_long$Celltype, levels = celltypes)
cell_frac_long$Sample <- factor(cell_frac_long$Sample, levels = sample_full_df$Sample)

# Plot Stacked Bar per Sample (Separated by Group)
p_stacked <- ggplot(cell_frac_long, aes(x = Sample, y = Fraction, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill", width = 0.9, alpha = 0.9) +
  facet_grid(~Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = color_res_coarse) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right"
  ) +
  labs(y = "Cell Fraction (%)")

# Plot Group Average Bar
group_avg$Celltype <- factor(group_avg$Celltype, levels = celltypes)
group_avg$Group <- factor(group_avg$Group, levels = c("Control", "Epilepsy"))

p_avg <- ggplot(group_avg, aes(x = Group, y = AverageFraction, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6, alpha = 0.9) +
  scale_fill_manual(values = color_res_coarse) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Average Abundance")

# Metadata heatmap
meta_long <- sample_full_df %>%
  dplyr::select(
    Sample, Group, Disease, BrainRegion, Gender,
    Age, OTS_Time, TotalSeizureBurden
  ) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(
    cols = c(-Sample, -Group),
    names_to = "Feature",
    values_to = "Value"
  )

meta_long$Sample <- factor(meta_long$Sample, levels = sample_full_df$Sample)
meta_long$Feature <- factor(
  meta_long$Feature,
  levels = rev(c("Disease", "Gender", "Age", "BrainRegion", "OTS_Time", "TotalSeizureBurden"))
)

p_meta <- ggplot(meta_long, aes(x = Sample, y = Feature, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(~Group, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )

p_combined <- (p_meta / p_stacked) + plot_layout(heights = c(1, 2))

ggsave(
  filename = file.path(out_check_dir, "check_FigS1B_Cell_Abundance_and_Metadata.png"),
  plot = p_combined,
  width = 12,
  height = 7,
  dpi = 300
)

cat("FigS1B completed successfully.\n")
