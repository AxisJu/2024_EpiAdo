# ==============================================================================
# Script: FigS8C_Gene_Permutation_Null_Distributions.R
# Description: Distribution of Cohen's d obtained by applying binary grouping logic
#              to 1,000 randomly selected background genes.
# Panel: Supplementary Figure 8C
# Source Data:
#   - FigS8C_Gene_Permutation_Distributions.csv
#   - FigS8C_Gene_Permutation_Statistics.csv
# Output:
#   - ../../sourcedata-260907/check_plots/check_FigS8C_Gene_Permutation_Null_Distributions.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- here::here()  # Project root; use here::here() or setwd() as appropriate
code_dir <- "."  # SuppFig8
data_dir <- file.path("..", "sourcedata", "SuppFig8")
plot_dir <- file.path("..", "sourcedata", "check_plots")

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

cache_path <- file.path(code_dir, "ps_s8_cache.rds")
s8_cache <- readRDS(cache_path)
permutation_data <- s8_cache$gene_permutation

# ==============================================================================
# Export Source Data
# ==============================================================================
stat_list <- list()
dist_list <- list()

for (prop in names(permutation_data)) {
  item <- permutation_data[[prop]]
  stat_list[[prop]] <- data.frame(
    Property = prop,
    Observed_Cohen_d = item$Obs_d,
    Standardized_Effect_Size = item$SES
  )
  dist_list[[prop]] <- data.frame(
    Gene_Index = 1:length(item$Dist),
    Property = prop,
    Random_Gene_Cohen_d = item$Dist
  )
}

stat_df <- bind_rows(stat_list)
dist_df <- bind_rows(dist_list)

write.csv(stat_df, file.path(data_dir, "FigS8C_Gene_Permutation_Statistics.csv"), row.names = FALSE)
write.csv(dist_df, file.path(data_dir, "FigS8C_Gene_Permutation_Distributions.csv"), row.names = FALSE)
cat("SourceData exported to:", data_dir, "\n")

# ==============================================================================
# Plot generation
# ==============================================================================
features <- c("AP_UpStroke", "AP_UpDownRatio", "AP_HalfWidth", "AP_Duration")
plot_list <- list()

for (i in seq_along(features)) {
  prop <- features[i]
  item <- permutation_data[[prop]]
  df_sub <- data.frame(d = item$Dist)
  obs_d <- item$Obs_d
  ses <- item$SES
  
  p <- ggplot(df_sub, aes(x = d)) +
    geom_histogram(bins = 40, fill = "grey80", color = "white") +
    geom_vline(xintercept = obs_d, color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Observed Cohen's d: %.3f\nStandardized Effect Size: %.3f", obs_d, ses),
             color = "black", vjust = 1.8, hjust = 1.05, fontface = "bold", size = 3.2) +
    theme_classic(base_size = 11) +
    labs(
      title = paste("Gene Permutation:", prop),
      x = "Cohen's d (Random Genes)",
      y = "Frequency"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(color = "black", size = 9)
    )
  
  plot_list[[i]] <- p
}

p_s8c <- wrap_plots(plot_list, nrow = 1)

check_png <- file.path(plot_dir, "check_FigS8C_Gene_Permutation_Null_Distributions.png")
ggsave(check_png, plot = p_s8c, width = 17, height = 4.2, dpi = 300)
cat("Check plot saved to:", check_png, "\n")
