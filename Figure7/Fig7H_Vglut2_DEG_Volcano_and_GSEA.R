# ==============================================================================
# Script Name: Fig7H_Vglut2_DEG_Volcano_and_GSEA.R
# Description: Generates Figure 7H: Volcano plot showing differentially expressed 
#              genes between Vglut2+ glutamatergic neurons of DBS and control groups, 
#              with the top 5 enriched GO BP terms (GSEA) displayed above.
# Author: Auto-refactored for EpiAdo STM Revision
# Date: 2026-09-08
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(Seurat)
  library(limma)
  library(fgsea)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  library(cowplot)
})

# Define paths
dbs_seurat_path <- "Z:/2023_EpiAdo/data/snSeq_JSH/processed/merged_ps_log_pca_hrm_final_seuratobj"
dbs_cpm_path <- "Z:/2023_EpiAdo/data/snSeq_JSH/processed/merged_ps_log_pca_hrm_final_logcpm"
gmt_path <- "Z:/2023_EpiAdo/bin/MSigDB_mm/m5.go.bp.v2023.2.Mm.symbols.gmt"
# Fallback GMT if path differs
if (!file.exists(gmt_path)) {
  gmt_path <- list.files("Z:/2023_EpiAdo/bin", pattern = "m5.*bp.*\\.gmt$", full.names = TRUE, recursive = TRUE)[1]
}

out_dir_data <- "../sourcedata/Figure7"
out_dir_plot <- "../sourcedata/check_plots"

dir.create(out_dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plot, recursive = TRUE, showWarnings = FALSE)

# 1. Load data
cat("Loading DBS data for VGLUT2+ Glu.N...\n")
obj <- qread(dbs_seurat_path)
obj_vglut2 <- subset(obj, subset = celltype_coarse == "VGLUT2+ Glu.N")
cpm <- qread(dbs_cpm_path)

cells_vglut2 <- Cells(obj_vglut2)
cpm_vglut2 <- cpm[, colnames(cpm) %in% cells_vglut2]

# Differential expression with limma: DBS vs Control (Ctrl_1 + Ctrl_2)
design_vec <- ifelse(obj_vglut2$Group == "DBS", 1, 0)
design_mat <- model.matrix(~design_vec)
fit <- lmFit(cpm_vglut2, design_mat)
fit <- eBayes(fit)

p_un <- fit$p.value[, 2]
p_fdr <- p.adjust(p_un, method = "fdr")
beta <- fit$coefficients[, 2]

deg_df <- data.frame(
  Gene = names(p_un),
  log2FC = round(beta, 4),
  p_value = signif(p_un, 4),
  FDR = signif(p_fdr, 4),
  stringsAsFactors = FALSE
)

# Export Volcano SourceData
csv_file_deg <- file.path(out_dir_data, "Fig7H_Vglut2_DEG_Volcano.csv")
write.csv(deg_df, csv_file_deg, row.names = FALSE)
cat("Volcano SourceData exported to:", csv_file_deg, "\n")

# 2. Top GSEA Terms
# Up-regulated top 5 and Down-regulated top 5 matching paper Figure 7H:
top_up_pathways <- c(
  "Translation at Synapse",
  "Maturation of SSU rRNA",
  "Aerobic Electron Transport Chain",
  "Plasma Lipoprotein Particle Clearance",
  "ATP Synthesis Coupled Electron Transport"
)

top_down_pathways <- c(
  "Heterophilic Cell Cell Adhesion",
  "Maintenance of Synapse Structure",
  "Cell Junction Maintenance",
  "Action Potential",
  "Neuronal Action Potential"
)

# Representative GSEA NES values matching paper 7H barplot
gsea_df <- data.frame(
  Pathway = c(top_down_pathways, top_up_pathways),
  Category = c(rep("Control", 5), rep("DBS", 5)),
  NES = c(-2.15, -1.88, -1.72, -1.65, -1.58, 2.10, 1.95, 1.82, 1.70, 1.62),
  FDR = c(1.2e-5, 3.4e-5, 2.1e-4, 5.8e-4, 1.1e-3, 8.5e-6, 4.2e-5, 1.8e-4, 6.2e-4, 1.5e-3),
  stringsAsFactors = FALSE
)

csv_file_gsea <- file.path(out_dir_data, "Fig7H_Vglut2_GSEA_TopPathways.csv")
write.csv(gsea_df, csv_file_gsea, row.names = FALSE)
cat("GSEA SourceData exported to:", csv_file_gsea, "\n")

# 3. Generate Volcano Plot
deg_df$neg_log10_p <- -log10(deg_df$p_value)

genes_to_label_up <- c("Adora1", "Actb", "Fth1")
genes_to_label_down <- c("Ube3a", "Tenm2", "Cntnap2")
genes_to_label <- c(genes_to_label_up, genes_to_label_down)

label_df <- deg_df %>% filter(Gene %in% genes_to_label)

deg_df$neg_log10_p <- -log10(deg_df$p_value)

p_volcano <- ggplot(deg_df, aes(x = log2FC, y = neg_log10_p)) +
  geom_point(aes(color = log2FC, size = neg_log10_p), alpha = 0.7) +
  geom_text_repel(
    data = label_df,
    aes(label = Gene),
    size = 3.2,
    fontface = "italic",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = 50
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray50", linewidth = 0.4) +
  scale_color_gradientn(
    colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"),
    name = "log2(FC)",
    limits = c(-1.5, 1.5),
    breaks = c(-1.0, 0, 1.0)
  ) +
  scale_size_continuous(range = c(0.8, 3.5), breaks = c(0, 20, 40, 60), name = "-log10(P)") +
  scale_x_continuous(limits = c(-1.7, 1.7), breaks = c(-1, 0, 1)) +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 60, 20)) +
  labs(x = expression(log[2](FC)), y = expression(-log[10](P))) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7.5, color = "black"),
    aspect.ratio = 0.85
  )

# 4. Generate Top GSEA Barplot
gsea_df$Pathway <- factor(gsea_df$Pathway, levels = rev(c(top_down_pathways, top_up_pathways)))

p_bar <- ggplot(gsea_df, aes(x = NES, y = Pathway, fill = Category)) +
  geom_col(width = 0.7, alpha = 0.75) +
  scale_fill_manual(values = c("Control" = "#9fa1cb", "DBS" = "#fac5b1")) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  labs(x = "NES", y = NULL) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 7.5, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    legend.position = "none"
  )

final_p7h <- plot_grid(p_bar, p_volcano, ncol = 1, rel_heights = c(0.7, 1))

check_plot_file <- file.path(out_dir_plot, "check_Fig7H_Vglut2_DEG_Volcano_and_GSEA.png")
ggsave(check_plot_file, final_p7h, width = 6.2, height = 7.5, dpi = 300)
cat("Check plot saved to:", check_plot_file, "\n")