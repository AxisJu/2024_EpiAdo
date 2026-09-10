# ==============================================================================
# Script Name: FigS2D_ADORA1_Leave_One_Out_Sensitivity.R
# Description: Scatter plot showing effect size of ADORA1 downregulation when
#              each individual patient/sample is removed (Supplementary Figure 2D).
# ==============================================================================

library(dplyr)
library(ggplot2)

out_data_dir <- "../sourcedata/SuppFig2"
out_check_dir <- "../sourcedata/check_plots"

cache_data <- file.path(out_data_dir, "ado_genes_glun_astro_meta_exp.rds")

if (!file.exists(cache_data)) {
  stop("Please run FigS2A first to generate ado_genes_glun_astro_meta_exp.rds cache.")
}

cat("Loading targeted data from cache...\n")
save_payload <- readRDS(cache_data)
meta_df <- save_payload$meta
sub_cpm <- save_payload$cpm

glun_cells <- rownames(meta_df)[meta_df$celltype_coarse == "Glu.N"]
glun_cells <- intersect(glun_cells, colnames(sub_cpm))

test_data <- data.frame(
  IsEpilepsy = meta_df[glun_cells, "Group"],
  Exp = sub_cpm["ADORA1", glun_cells],
  Sample = meta_df[glun_cells, "Sample_snSeq.processed"]
) %>%
  dplyr::filter(Exp > 0)

test_data$IsEpilepsy <- factor(test_data$IsEpilepsy, levels = c("Control", "Epilepsy"))

all_samples <- sort(unique(test_data$Sample))

cat("Running leave-one-out sensitivity models...\n")
loo_results <- lapply(all_samples, function(s_to_remove) {
  sub_data <- test_data %>% filter(Sample != s_to_remove)
  fit <- lm(Exp ~ IsEpilepsy, data = sub_data)
  s <- summary(fit)$coefficients
  ci <- confint(fit)
  
  data.frame(
    Removed_Sample = s_to_remove,
    Estimate = s["IsEpilepsyEpilepsy", "Estimate"],
    PValue = s["IsEpilepsyEpilepsy", "Pr(>|t|)"],
    LCI = ci["IsEpilepsyEpilepsy", 1],
    UCI = ci["IsEpilepsyEpilepsy", 2]
  )
}) %>% bind_rows()

# Baseline Full Data
full_fit <- lm(Exp ~ IsEpilepsy, data = test_data)
full_s <- summary(full_fit)$coefficients
full_ci <- confint(full_fit)

full_res <- data.frame(
  Removed_Sample = "None (Full Data)",
  Estimate = full_s["IsEpilepsyEpilepsy", "Estimate"],
  PValue = full_s["IsEpilepsyEpilepsy", "Pr(>|t|)"],
  LCI = full_ci["IsEpilepsyEpilepsy", 1],
  UCI = full_ci["IsEpilepsyEpilepsy", 2]
)

loo_plot_df <- rbind(full_res, loo_results)

# Export SourceData
cat("Exporting SourceData...\n")
write.csv(
  loo_plot_df,
  file = file.path(out_data_dir, "FigS2D_ADORA1_Leave_One_Out_Sensitivity.csv"),
  row.names = FALSE
)

cat("Plotting Check Figures...\n")
loo_plot_df$Removed_Sample <- factor(
  loo_plot_df$Removed_Sample,
  levels = c("None (Full Data)", all_samples)
)

# Plotting with exact author style from 260205_revSTM.R:
# geom_rect for full data CI, error bars in #9fa1cb, points in #2c7fb8, black diamond for baseline
p_loo <- ggplot(loo_plot_df, aes(x = Removed_Sample, y = Estimate)) +
  geom_rect(
    aes(ymin = full_res$LCI, ymax = full_res$UCI, xmin = -Inf, xmax = Inf),
    fill = "grey90", alpha = 0.3, inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.25, color = "#9fa1cb", linewidth = 0.6) +
  geom_point(aes(color = (PValue < 0.05)), size = 2) +
  geom_point(
    data = filter(loo_plot_df, Removed_Sample == "None (Full Data)"),
    color = "black", size = 3, shape = 18
  ) +
  scale_color_manual(values = c("TRUE" = "#2c7fb8", "FALSE" = "grey70"), labels = c("TRUE" = "P < 0.05")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  labs(
    x = NULL,
    y = "Effect Size (Epilepsy - Control)"
  ) +
  ylim(-0.075, 0.005)

ggsave(
  filename = file.path(out_check_dir, "check_FigS2D_ADORA1_Leave_One_Out.png"),
  plot = p_loo,
  width = 13,
  height = 4.5,
  dpi = 300
)

cat("FigS2D completed successfully.\n")
