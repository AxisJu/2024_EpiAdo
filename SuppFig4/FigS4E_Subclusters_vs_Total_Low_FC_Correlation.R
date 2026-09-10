# ==============================================================================
# Script Name: FigS4E_Subclusters_vs_Total_Low_FC_Correlation.R
# Description: Scatter plot showing strong positive correlation between the log2
#              fold-change of genes in the significantly enriched subclusters and
#              the total ADORA1-Low population (Supplementary Figure 4E).
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggpubr)

out_data_dir <- "../sourcedata/SuppFig4"
out_check_dir <- "../sourcedata/check_plots"

cache_path <- file.path(out_data_dir, "glun_s4_cache.rds")
if (!file.exists(cache_path)) {
  stop("glun_s4_cache.rds not found. Please ensure cache is built.")
}

cat("Loading SuppFig4 cache...\n")
s4_cache <- readRDS(cache_path)
corr_df <- s4_cache$corr_df

# Filter floor values (< -10)
corr_df_clean <- corr_df %>%
  filter(Total_Low_FC > -10 & Subcluster_FC > -10)

# Calculate Pearson correlation
cor_res <- cor.test(corr_df_clean$Total_Low_FC, corr_df_clean$Subcluster_FC, method = "pearson")
cat(sprintf("Pearson correlation: R = %.2f, p-value = %e\n", cor_res$estimate, cor_res$p.value))

# Export SourceData
cat("Exporting SourceData for Fig S4E...\n")
write.csv(
  corr_df_clean %>%
    dplyr::select(
      Gene,
      Total_Low_log2FC = Total_Low_FC,
      Significant_Subclusters_log2FC = Subcluster_FC
    ),
  file = file.path(out_data_dir, "FigS4E_Subclusters_vs_Total_Low_FC_Correlation.csv"),
  row.names = FALSE
)

# Generate Check Plot
cat("Generating Check Plot for Fig S4E...\n")
p_s4e <- ggplot(corr_df_clean, aes(x = Total_Low_FC, y = Subcluster_FC)) +
  geom_point(alpha = 0.3, size = 0.8, color = "grey50") +
  geom_point(
    data = subset(corr_df_clean, abs(Total_Low_FC) > 1),
    aes(color = Subcluster_FC),
    alpha = 0.8,
    size = 1
  ) +
  scale_color_gradient2(
    low = "#4575b4",
    mid = "grey80",
    high = "#d73027",
    name = "Subcluster_FC"
  ) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = FALSE) +
  stat_cor(
    method = "pearson",
    size = 5.5,
    label.x = -2.5,
    label.y = 4.8
  ) +
  theme_bw() +
  labs(
    x = expression(log[2]*FC ~ "(Total " * italic(ADORA1)^Low * ")"),
    y = expression(log[2]*FC ~ "(Significant Subclusters)")
  ) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  )

ggsave(
  filename = file.path(out_check_dir, "check_FigS4E_Subclusters_vs_Total_Low_Correlation.png"),
  plot = p_s4e,
  width = 6.5,
  height = 5.5,
  dpi = 300
)

cat("Fig S4E completed successfully.\n")
