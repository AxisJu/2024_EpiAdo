# ==============================================================================
# Script Name: FigS1E_Metabolic_Reprogramming_Magnitude_Boxplot.R
# Description: Boxplots comparing the absolute Cohen's d of metabolic
#              reprogramming across cell types (Supplementary Figure 1E).
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

out_data_dir <- "../sourcedata/SuppFig1"
out_check_dir <- "../sourcedata/check_plots"

cache_stat <- file.path(out_data_dir, "stat_res_all_regions.rds")

if (!file.exists(cache_stat)) {
  stop("Please run FigS1C first to generate stat_res_all_regions.rds cache.")
}

cat("Loading cached statistics...\n")
stat_res <- readRDS(cache_stat)

stat_res$AbsCohenD <- abs(stat_res$CohenD)
stat_res$BrainRegion <- factor(
  stat_res$BrainRegion,
  levels = c("Total", "Cortex", "Hippocampus", "Amygdala")
)
stat_res$Celltype <- factor(stat_res$Celltype, levels = c("Glu.N", "GABA.N", "Astro."))

# Export SourceData
cat("Exporting SourceData...\n")
write.csv(
  stat_res,
  file = file.path(out_data_dir, "FigS1E_Metabolic_Reprogramming_AbsCohenD.csv"),
  row.names = FALSE
)

# Pairwise Wilcoxon tests per region
my_comparisons <- list(c("Glu.N", "Astro."), c("GABA.N", "Astro."), c("Glu.N", "GABA.N"))
stat_compare <- stat_res %>%
  group_by(BrainRegion) %>%
  rstatix::wilcox_test(AbsCohenD ~ Celltype, comparisons = my_comparisons) %>%
  rstatix::add_significance()

write.csv(
  stat_compare,
  file = file.path(out_data_dir, "FigS1E_Metabolic_Reprogramming_Statistics.csv"),
  row.names = FALSE
)

cat("Plotting Check Figures...\n")
color_celltype <- c("Glu.N" = "#d0a9cb", "GABA.N" = "#9dcb84", "Astro." = "#adaed3")

p_box <- ggplot(stat_res, aes(x = Celltype, y = AbsCohenD, fill = Celltype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6, color = "black", linewidth = 0.3) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.2, width = 0.2) +
  facet_wrap(~BrainRegion, ncol = 4, scales = "fixed") +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    step_increase = 0.1
  ) +
  scale_fill_manual(values = color_celltype) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey95"),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Degree of Perturbation (|Cohen's d|)",
    title = "Magnitude of Metabolic Reprogramming"
  )

ggsave(
  filename = file.path(out_check_dir, "check_FigS1E_Reprogramming_Magnitude_Boxplot.png"),
  plot = p_box,
  width = 10,
  height = 5,
  dpi = 300
)

cat("FigS1E completed successfully.\n")
