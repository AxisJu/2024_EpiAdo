# ==============================================================================
# Script Name: FigS2A_GluN_Adenosine_Metabolism_Genes_Violin.R
# Description: Violin plots of logCPM for key adenosine metabolism genes in
#              glutamatergic neurons (Supplementary Figure 2A).
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)

out_data_dir <- "../sourcedata/SuppFig2"
out_check_dir <- "../sourcedata/check_plots"

dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_check_dir, recursive = TRUE, showWarnings = FALSE)

cache_data <- file.path(out_data_dir, "ado_genes_glun_astro_meta_exp.rds")

ado_genes_all <- c(
  "ADORA1", "ADORA2A", "ENTPD1", "ENPP1", "NT5E", "NT5C2",
  "SLC29A1", "SLC29A2", "SLC28A2", "SLC28A3", "ADK", "ADA"
)

if (!file.exists(cache_data)) {
  library(Seurat)
  library(qs)
  obj_path <- "Z:/2023_EpiAdo/data/scSeq_MSH/Merged/merged_ps_log_pca_harmony_umap_final_seuratobj"
  cpm_path <- "Z:/2023_EpiAdo/data/scSeq_MSH/Merged/merged_ps_log_pca_harmony_umap_final_logcpm"
  
  cat("Loading Seurat metadata and CPM matrix for targeted genes...\n")
  merged_seuratobj <- qread(obj_path)
  cpm <- qread(cpm_path)
  
  avail_genes <- intersect(ado_genes_all, rownames(cpm))
  sub_cpm <- as.matrix(cpm[avail_genes, ])
  rm(cpm); gc()
  
  meta_df <- merged_seuratobj@meta.data
  rm(merged_seuratobj); gc()
  
  save_payload <- list(meta = meta_df, cpm = sub_cpm)
  saveRDS(save_payload, file = cache_data)
  cat("Saved shared cache successfully.\n")
} else {
  cat("Loading targeted data from cache...\n")
  save_payload <- readRDS(cache_data)
  meta_df <- save_payload$meta
  sub_cpm <- save_payload$cpm
}

target_genes <- c("ENTPD1", "ENPP1", "NT5E", "NT5C2", "SLC29A1", "SLC29A2", "SLC28A2", "SLC28A3", "ADK", "ADA")
target_genes <- intersect(target_genes, rownames(sub_cpm))

# Subset to Glu.N cells
glun_cells <- rownames(meta_df)[meta_df$celltype_coarse == "Glu.N"]
glun_cells <- intersect(glun_cells, colnames(sub_cpm))

cat("Extracting Glu.N single-cell data...\n")
expr_list <- list()
stat_list <- list()

for (gene in target_genes) {
  df_gene <- data.frame(
    CellID = glun_cells,
    Gene = gene,
    Group = meta_df[glun_cells, "Group"],
    Sample = meta_df[glun_cells, "Sample_snSeq.processed"],
    Exp = sub_cpm[gene, glun_cells]
  ) %>%
    dplyr::filter(Exp > 0)
  
  df_gene$Group <- factor(df_gene$Group, levels = c("Control", "Epilepsy"))
  expr_list[[gene]] <- df_gene
  
  # Wilcoxon rank-sum test
  w_res <- wilcox.test(Exp ~ Group, data = df_gene)
  p_val <- w_res$p.value
  
  # Label format matching paper
  if (p_val < 1e-100) {
    p_txt <- "P < 1 x 10^-100"
  } else if (p_val < 1e-20) {
    p_txt <- "P < 1 x 10^-20"
  } else if (p_val < 1e-7) {
    p_txt <- "P < 1 x 10^-07"
  } else {
    p_txt <- sprintf("P = %.3f", p_val)
  }
  
  stat_list[[gene]] <- data.frame(
    Gene = gene,
    Control_Cells = sum(df_gene$Group == "Control"),
    Epilepsy_Cells = sum(df_gene$Group == "Epilepsy"),
    Control_Mean = mean(df_gene$Exp[df_gene$Group == "Control"]),
    Epilepsy_Mean = mean(df_gene$Exp[df_gene$Group == "Epilepsy"]),
    PValue = p_val,
    Label = p_txt
  )
}

expr_df <- do.call(rbind, expr_list)
stat_df <- do.call(rbind, stat_list)

# Subsample for SourceData CSV export
set.seed(42)
expr_sub <- expr_df %>%
  group_by(Gene, Group) %>%
  sample_n(min(n(), 5000)) %>%
  ungroup()

# Export SourceData
cat("Exporting SourceData...\n")
write.csv(
  expr_sub,
  file = file.path(out_data_dir, "FigS2A_GluN_Ado_Genes_SingleCell_Expression.csv"),
  row.names = FALSE
)

write.csv(
  stat_df,
  file = file.path(out_data_dir, "FigS2A_GluN_Ado_Genes_Statistics.csv"),
  row.names = FALSE
)

cat("Plotting Check Figures...\n")
# Exact colors from 250211 ChangeColor.R: Control=#989798, Epilepsy=#fac5b1
color_group <- c("Control" = "#989798", "Epilepsy" = "#fac5b1")

plots <- list()
for (i in seq_along(target_genes)) {
  gene <- target_genes[i]
  df_g <- expr_df %>% dplyr::filter(Gene == gene)
  st_g <- stat_df %>% dplyr::filter(Gene == gene)
  p_label <- st_g$Label[1]
  
  max_y <- max(df_g$Exp, na.rm = TRUE)
  bracket_y <- max_y * 1.08
  text_y <- max_y * 1.16
  upper_lim <- max_y * 1.25
  
  p <- ggplot(df_g, aes(x = Group, y = Exp, fill = Group, color = Group)) +
    geom_violin(linewidth = 0.4, alpha = 0.4, width = 0.8, trim = TRUE) +
    scale_fill_manual(values = color_group) +
    scale_color_manual(values = color_group) +
    # Bracket and P-value label
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, color = "black", linewidth = 0.3) +
    annotate("segment", x = 1, xend = 1, y = bracket_y - max_y*0.03, yend = bracket_y, color = "black", linewidth = 0.3) +
    annotate("segment", x = 2, xend = 2, y = bracket_y - max_y*0.03, yend = bracket_y, color = "black", linewidth = 0.3) +
    annotate("text", x = 1.5, y = text_y, label = p_label, size = 2.4, color = "black") +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = if (i == 1) element_text(size = 9, color = "black") else element_blank(),
      axis.text.y = element_text(size = 8, color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      legend.position = if (i == length(target_genes)) "right" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      plot.margin = unit(c(0.2, 0.1, 0.2, 0.1), "cm")
    ) +
    labs(title = gene, y = "log2CPM") +
    ylim(0, upper_lim)
  
  plots[[gene]] <- p
}

p_combined <- wrap_plots(plots, nrow = 1)
ggsave(
  filename = file.path(out_check_dir, "check_FigS2A_Ado_Genes_Violin.png"),
  plot = p_combined,
  width = 17,
  height = 3.5,
  dpi = 300
)

cat("FigS2A completed successfully.\n")
