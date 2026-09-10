# ==============================================================================
# Script Name: FigS1D_Sig_Pathways_Count_Barplot.R
# Description: Bar charts showing total number of significantly upregulated (red)
#              and downregulated (blue) metabolic pathways (FDR < 0.05) across
#              regions (Total, Cortex, HPC, AMY) (Supplementary Figure 1D).
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

out_data_dir <- "../sourcedata/SuppFig1"
out_check_dir <- "../sourcedata/check_plots"

cache_stat <- file.path(out_data_dir, "stat_res_all_regions.rds")

if (!file.exists(cache_stat)) {
  stop("Please run FigS1C first to generate stat_res_all_regions.rds cache.")
}

cat("Loading cached statistics...\n")
stat_res <- readRDS(cache_stat)

# Count significant pathways
stat_res$BrainRegion <- factor(
  stat_res$BrainRegion,
  levels = c("Total", "Cortex", "Hippocampus", "Amygdala"),
  labels = c("Total", "Cortex", "HPC", "AMY")
)
stat_res$Celltype <- factor(stat_res$Celltype, levels = c("Glu.N", "GABA.N", "Astro."))

sig_counts <- stat_res %>%
  dplyr::filter(AdjP < 0.05) %>%
  dplyr::mutate(Direction = ifelse(CohenD > 0, "Upregulated", "Downregulated")) %>%
  dplyr::group_by(BrainRegion, Celltype, Direction) %>%
  dplyr::summarise(Count = n(), .groups = "drop")

# Export SourceData
cat("Exporting SourceData...\n")
write.csv(
  sig_counts,
  file = file.path(out_data_dir, "FigS1D_Significant_Pathways_Counts.csv"),
  row.names = FALSE
)

write.csv(
  stat_res %>% dplyr::filter(AdjP < 0.05),
  file = file.path(out_data_dir, "FigS1D_Significant_Pathways_Details.csv"),
  row.names = FALSE
)

cat("Plotting Check Figures...\n")
color_direction <- c("Upregulated" = "#D6604D", "Downregulated" = "#4393C3")

p_bar <- ggplot(sig_counts, aes(x = Celltype, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", linewidth = 0.2) +
  facet_wrap(~BrainRegion, ncol = 4, scales = "fixed") +
  scale_fill_manual(values = color_direction) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  ) +
  labs(
    x = NULL,
    y = "Number of Significant Pathways (FDR < 0.05)"
  )

ggsave(
  filename = file.path(out_check_dir, "check_FigS1D_Sig_Pathways_Barplot.png"),
  plot = p_bar,
  width = 9,
  height = 5,
  dpi = 300
)

cat("FigS1D completed successfully.\n")
