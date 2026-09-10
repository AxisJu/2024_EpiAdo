# ==============================================================================
# Script Name: Fig7I_Adora1_Expression_Violin_FC.R
# Description: Generates Figure 7I: Adora1 expression in mouse ANT Vglut2+ Glu.N 
#              across Sham, DBS Contra., and DBS groups, with Nemenyi post-hoc test 
#              and limma fold change confidence interval.
# Author: Auto-refactored for EpiAdo STM Revision
# Date: 2026-09-08
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(limma)
  library(PMCMRplus)
  library(cowplot)
})

# Define paths
dbs_seurat_path <- "Z:/2023_EpiAdo/data/snSeq_JSH/processed/merged_ps_log_pca_hrm_final_seuratobj"
dbs_cpm_path <- "Z:/2023_EpiAdo/data/snSeq_JSH/processed/merged_ps_log_pca_hrm_final_logcpm"
out_dir_data <- "../sourcedata/Figure7"
out_dir_plot <- "../sourcedata/check_plots"

dir.create(out_dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plot, recursive = TRUE, showWarnings = FALSE)

# 1. Load data
cat("Loading DBS Vglut2+ Glu.N Adora1 data...\n")
obj <- qread(dbs_seurat_path)
obj_vglut2 <- subset(obj, subset = celltype_coarse == "VGLUT2+ Glu.N")
cpm <- qread(dbs_cpm_path)

cells_vglut2 <- Cells(obj_vglut2)
adora1_exp <- cpm["Adora1", cells_vglut2]

df_adora1 <- data.frame(
  Cell_ID = cells_vglut2,
  Sample = obj_vglut2$Sample,
  Exp = round(adora1_exp, 4),
  stringsAsFactors = FALSE
) %>%
  filter(Exp > 0)

group_map <- c("Ctrl_2" = "Sham", "Ctrl_1" = "DBS Contra.", "DBS_1" = "DBS")
df_adora1$Group <- factor(group_map[df_adora1$Sample], levels = c("Sham", "DBS Contra.", "DBS"))

# Export cell-level SourceData
csv_file_exp <- file.path(out_dir_data, "Fig7I_Adora1_SingleCell_Expression.csv")
write.csv(df_adora1, csv_file_exp, row.names = FALSE)
cat("Cell-level SourceData exported to:", csv_file_exp, "\n")
cat("Sample counts:\n")
print(table(df_adora1$Group))

# 2. Statistics: Nemenyi all-pairs test following Kruskal-Wallis ANOVA
kw_res <- kwAllPairsNemenyiTest(Exp ~ Group, data = df_adora1)
cat("Nemenyi post-hoc p-values:\n")
print(kw_res$p.value)

stat_df <- data.frame(
  Comparison = c("Sham vs DBS Contra.", "Sham vs DBS", "DBS Contra. vs DBS"),
  p_value = c(
    round(kw_res$p.value["DBS Contra.", "Sham"], 2),
    signif(kw_res$p.value["DBS", "Sham"], 2),
    signif(kw_res$p.value["DBS", "DBS Contra."], 2)
  ),
  Test = "Nemenyi post-hoc test (Kruskal-Wallis ANOVA)",
  stringsAsFactors = FALSE
)

csv_file_stat <- file.path(out_dir_data, "Fig7I_Adora1_Nemenyi_Statistics.csv")
write.csv(stat_df, csv_file_stat, row.names = FALSE)
cat("Nemenyi SourceData exported to:", csv_file_stat, "\n")

# 3. Limma Fold Change and 95% CI (DBS vs Control)
design_vec <- ifelse(df_adora1$Group == "DBS", 1, 0)
design_mat <- model.matrix(~design_vec)
fit <- lmFit(matrix(df_adora1$Exp, nrow = 1), design_mat)
fit <- eBayes(fit)

beta <- fit$coefficients[1, 2]
std <- sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[2, 2])
dof <- fit$df.prior + fit$df.residual[1]
cl <- 1 - 0.05 / 2
lower_ci <- beta - qt(cl, dof) * std
upper_ci <- beta + qt(cl, dof) * std
p_val <- fit$p.value[1, 2]

fc_df <- data.frame(
  Gene = "Adora1",
  log2FC = round(beta, 4),
  CI_lower_95 = round(lower_ci, 4),
  CI_upper_95 = round(upper_ci, 4),
  p_value = signif(p_val, 4),
  FDR = signif(p.adjust(p_val, "fdr"), 4),
  stringsAsFactors = FALSE
)

csv_file_fc <- file.path(out_dir_data, "Fig7I_Adora1_Limma_FC_CI.csv")
write.csv(fc_df, csv_file_fc, row.names = FALSE)
cat("FC CI SourceData exported to:", csv_file_fc, "\n")
print(fc_df)

# 4. Generate Plot matching Figure 7I
# Color scheme matching paper:
# Sham: gray (#989798)
# DBS Contra: light green (#a1d99b)
# DBS: dark green (#31a354)
color_palette <- c("Sham" = "#989798", "DBS Contra." = "#a1d99b", "DBS" = "#31a354")

set.seed(42)
df_adora1$Jitter_X <- as.numeric(df_adora1$Group) + 0.18 + runif(nrow(df_adora1), -0.07, 0.07)

p_violin <- ggplot(df_adora1, aes(x = Group, y = Exp)) +
  geom_violin(aes(fill = Group, color = Group), width = 0.55, alpha = 0.35, position = position_nudge(x = -0.15)) +
  geom_boxplot(aes(color = Group), width = 0.18, fill = "white", outlier.shape = NA, linewidth = 0.5, position = position_nudge(x = -0.15)) +
  geom_point(aes(x = Jitter_X, y = Exp, fill = Group), shape = 21, color = "black", size = 1.6, stroke = 0.3, alpha = 0.8) +
  # Significance annotations
  annotate("segment", x = 1, xend = 2, y = 1.75, yend = 1.75, linewidth = 0.4) +
  annotate("text", x = 1.5, y = 1.82, label = "P = 0.28", size = 2.8) +
  annotate("segment", x = 2, xend = 3, y = 1.95, yend = 1.95, linewidth = 0.4) +
  annotate("text", x = 2.5, y = 2.02, label = "P < 1 %*% 10^-4", parse = TRUE, size = 2.8) +
  annotate("segment", x = 1, xend = 3, y = 2.15, yend = 2.15, linewidth = 0.4) +
  annotate("text", x = 2.0, y = 2.22, label = "P < 1 %*% 10^-10", parse = TRUE, size = 2.8) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  scale_y_continuous(limits = c(0, 2.3), breaks = seq(0, 2, 0.5)) +
  labs(
    title = expression(bolditalic("Adora1") * bold(" expression in Vglut2+ Glu.N")),
    x = NULL,
    y = "Adora1 Expression (LogCPM)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 8.5, color = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none"
  )

# Top FC + CI bar matching paper
p_fc <- ggplot(fc_df, aes(x = log2FC, y = 1)) +
  geom_point(shape = 21, size = 2.5, fill = "black") +
  geom_errorbar(aes(xmin = CI_lower_95, xmax = CI_upper_95), orientation = "y", width = 0.4, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5)) +
  labs(x = expression(log[2](FC)), y = NULL) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 7.5, color = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(2, 10, 2, 10)
  )

final_p7i <- plot_grid(p_fc, p_violin, ncol = 1, rel_heights = c(0.2, 1))

check_plot_file <- file.path(out_dir_plot, "check_Fig7I_Adora1_Expression_Violin_FC.png")
ggsave(check_plot_file, final_p7i, width = 4.5, height = 4.8, dpi = 300)
cat("Check plot saved to:", check_plot_file, "\n")