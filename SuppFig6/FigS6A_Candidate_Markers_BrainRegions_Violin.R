# ==============================================================================
# Script: FigS6A_Candidate_Markers_BrainRegions_Violin.R
# Description: Identification and verification of markers for ADORA1^Low Glu.N
#              across Amygdala, Cortex, and Hippocampus (Single-Cell & Pseudo-Bulk).
# Panel: Supplementary Figure 6A
# Source Data:
#   - FigS6A_Candidate_Markers_SC_Expression.csv
#   - FigS6A_Candidate_Markers_PB_Sample_Expression.csv
#   - FigS6A_Candidate_Markers_Effect_Sizes_Statistics.csv
# Output:
#   - ../../sourcedata-260907/check_plots/check_FigS6A_Candidate_Markers_BrainRegions_Violin.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(effsize)
  library(ggrastr)
})

# Set working paths
base_dir <- here::here()  # Project root; use here::here() or setwd() as appropriate
code_dir <- "."  # SuppFig6
data_dir <- file.path("..", "sourcedata", "SuppFig6")
plot_dir <- file.path("..", "sourcedata", "check_plots")

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Load lightweight pre-extracted single-cell cache
cache_path <- file.path(code_dir, "glun_s6_cache.rds")
if (!file.exists(cache_path)) {
  stop("Cache file not found: ", cache_path)
}
cache <- readRDS(cache_path)
meta <- cache$meta
expr_mat <- cache$expr

target_genes <- c("BMPER", "GCSH", "RTN4")
regions <- c("Amygdala", "Cortex", "Hippocampus")
grp_levels <- c("Control", "Epi_ADORA1_Low", "Epi_ADORA1_High")
cols_grp <- c("Control" = "#989798", "Epi_ADORA1_Low" = "#7BC5B4", "Epi_ADORA1_High" = "#205A9E")

set.seed(42)

df_sc_list <- list()
df_pb_list <- list()
df_pb_mean_list <- list()
df_pval_list <- list()

for (g in target_genes) {
  for (reg in regions) {
    cells_in_reg <- rownames(meta)[meta$BrainRegion == reg & meta$Compare_Group %in% grp_levels]
    if (length(cells_in_reg) == 0) next
    
    expr_all <- as.numeric(expr_mat[g, cells_in_reg])
    groups_all <- meta[cells_in_reg, "Compare_Group"]
    samples_all <- meta[cells_in_reg, "Sample_snSeq.processed"]
    
    # 1. 1. Single-cell (SC) data
    tmp_sc <- data.frame(
      Cell_ID = cells_in_reg,
      Sample = samples_all,
      Gene = g,
      Region = reg,
      Group = factor(groups_all, levels = grp_levels),
      Expression = expr_all
    )
    tmp_sc$X_Num <- as.numeric(tmp_sc$Group) - 0.15
    tmp_sc$X_Jitter <- tmp_sc$X_Num + runif(nrow(tmp_sc), -0.08, 0.08)
    df_sc_list[[length(df_sc_list) + 1]] <- tmp_sc
    
    # 2. 2. Pseudo-bulk (PB) sample-level mean data
    tmp_pb <- data.frame(
      Expression = expr_all,
      Sample = samples_all,
      Group = factor(groups_all, levels = grp_levels)
    ) %>%
      group_by(Sample, Group) %>%
      summarize(Expression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
      mutate(Gene = g, Region = reg)
    
    tmp_pb$X_Num <- as.numeric(tmp_pb$Group) + 0.15
    tmp_pb$X_Jitter <- tmp_pb$X_Num + runif(nrow(tmp_pb), -0.03, 0.03)
    df_pb_list[[length(df_pb_list) + 1]] <- tmp_pb
    
    # 3. Group-level PB grand mean
    tmp_pb_mean <- tmp_pb %>%
      group_by(Group) %>%
      summarize(Mean_Exp = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
      mutate(Gene = g, Region = reg, X_Num = as.numeric(Group) + 0.15)
    df_pb_mean_list[[length(df_pb_mean_list) + 1]] <- tmp_pb_mean
    
    # 4. Statistical tests (Cohen's d and Wilcoxon test)
    for (target_ref in c("Control", "Epi_ADORA1_High")) {
      sc_low <- tmp_sc$Expression[tmp_sc$Group == "Epi_ADORA1_Low"]
      sc_ref <- tmp_sc$Expression[tmp_sc$Group == target_ref]
      
      pb_low <- tmp_pb$Expression[tmp_pb$Group == "Epi_ADORA1_Low"]
      pb_ref <- tmp_pb$Expression[tmp_pb$Group == target_ref]
      
      p_sc <- wilcox.test(sc_low, sc_ref)$p.value
      d_sc <- ifelse(length(sc_low) >= 3 & length(sc_ref) >= 3, cohen.d(sc_low, sc_ref)$estimate, NA)
      
      p_pb <- NA; d_pb <- NA
      if (length(pb_low) >= 3 & length(pb_ref) >= 3) {
        p_pb <- wilcox.test(pb_low, pb_ref, exact = FALSE)$p.value
        d_pb <- cohen.d(pb_low, pb_ref)$estimate
      }
      
      fmt_star <- function(p) {
        ifelse(is.na(p), "NA", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
      }
      
      label_text <- paste0(
        "SC: d=", ifelse(is.na(d_sc), "NA", sprintf("%.2f", d_sc)), " (", fmt_star(p_sc), ")\n",
        "PB: d=", ifelse(is.na(d_pb), "NA", sprintf("%.2f", d_pb)), " (", fmt_star(p_pb), ")"
      )
      
      center_x <- ifelse(target_ref == "Control", 1.5, 2.5)
      
      df_pval_list[[length(df_pval_list) + 1]] <- data.frame(
        Gene = g,
        Region = reg,
        Comparison = ifelse(target_ref == "Control", "Low_vs_Ctrl", "Low_vs_High"),
        SC_Cohen_d = d_sc,
        SC_P_Value = p_sc,
        SC_Significance = fmt_star(p_sc),
        PB_Cohen_d = d_pb,
        PB_P_Value = p_pb,
        PB_Significance = fmt_star(p_pb),
        label_text = label_text,
        Center_X = center_x
      )
    }
  }
}

plot_data_sc <- bind_rows(df_sc_list) %>%
  mutate(
    Gene = factor(Gene, levels = target_genes),
    Region = factor(Region, levels = regions)
  )

plot_data_pb <- bind_rows(df_pb_list) %>%
  mutate(
    Gene = factor(Gene, levels = target_genes),
    Region = factor(Region, levels = regions)
  )

plot_data_pb_mean <- bind_rows(df_pb_mean_list) %>%
  mutate(
    Gene = factor(Gene, levels = target_genes),
    Region = factor(Region, levels = regions)
  )

plot_data_pval <- bind_rows(df_pval_list) %>%
  mutate(
    Gene = factor(Gene, levels = target_genes),
    Region = factor(Region, levels = regions)
  )

# ==============================================================================
# Export Source Data
# ==============================================================================
sc_export <- plot_data_sc %>%
  dplyr::select(Cell_ID, Sample, Region, Gene, Group, Expression)
write.csv(sc_export, file.path(data_dir, "FigS6A_Candidate_Markers_SC_Expression.csv"), row.names = FALSE)

pb_export <- plot_data_pb %>%
  dplyr::select(Sample, Region, Gene, Group, Expression)
write.csv(pb_export, file.path(data_dir, "FigS6A_Candidate_Markers_PB_Sample_Expression.csv"), row.names = FALSE)

stat_export <- plot_data_pval %>%
  dplyr::select(Gene, Region, Comparison, SC_Cohen_d, SC_P_Value, SC_Significance, PB_Cohen_d, PB_P_Value, PB_Significance)
write.csv(stat_export, file.path(data_dir, "FigS6A_Candidate_Markers_Effect_Sizes_Statistics.csv"), row.names = FALSE)

cat("Source data successfully exported to:", data_dir, "\n")

# ==============================================================================
# Plot generation (matching published figure)
# ==============================================================================
p_expr <- ggplot() +
  # 1. SC violin
  geom_violin(data = plot_data_sc, aes(x = X_Num, y = Expression, fill = Group, group = X_Num),
              alpha = 0.35, color = NA, width = 0.26) +
  # 2. SC points (rasterized for performance)
  geom_point_rast(data = plot_data_sc, aes(x = X_Jitter, y = Expression, color = Group),
                  size = 0.35, alpha = 0.25, raster.dpi = 300) +
  # 3. PB individual sample mean points
  geom_point(data = plot_data_pb, aes(x = X_Jitter, y = Expression, color = Group),
             size = 1.4, alpha = 0.45, shape = 16) +
  # 4. PB group mean points (solid diamond)
  geom_point(data = plot_data_pb_mean, aes(x = X_Num, y = Mean_Exp, fill = Group),
             shape = 23, size = 4.0, color = "black", stroke = 1.0) +
  # 5. Statistical annotation labels
  geom_text(data = plot_data_pval, aes(x = Center_X, y = Inf, label = label_text),
            vjust = 1.25, size = 3.1, fontface = "bold", lineheight = 1.05) +
  facet_grid(Gene ~ Region, scales = "free_y") +
  scale_fill_manual(values = cols_grp) +
  scale_color_manual(values = cols_grp) +
  scale_x_continuous(
    breaks = 1:3,
    labels = c("Control", expression(paste("Epi. ", ADORA1^Low)), expression(paste("Epi. ", ADORA1^High)))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Normalized Expression") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.8),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
    legend.position = "none",
    panel.spacing.x = unit(2.2, "lines"),
    panel.spacing.y = unit(0.6, "lines")
  )

check_png <- file.path(plot_dir, "check_FigS6A_Candidate_Markers_BrainRegions_Violin.png")
ggsave(check_png, plot = p_expr, width = 10.5, height = 7.5, dpi = 300)
cat("Check plot saved to:", check_png, "\n")
