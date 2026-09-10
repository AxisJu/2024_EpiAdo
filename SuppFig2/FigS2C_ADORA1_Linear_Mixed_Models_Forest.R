# ==============================================================================
# Script Name: FigS2C_ADORA1_Linear_Mixed_Models_Forest.R
# Description: Forest plot displaying estimated effect size of ADORA1 expression
#              across four progressively rigorous linear models (Supplementary Figure 2C).
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
  Disease = meta_df[glun_cells, "Disease"],
  Exp = sub_cpm["ADORA1", glun_cells],
  Age = as.numeric(meta_df[glun_cells, "Age"]),
  Gender = meta_df[glun_cells, "Gender"],
  Platform = meta_df[glun_cells, "Platform"],
  Cohort = meta_df[glun_cells, "Cohort"],
  Sample = meta_df[glun_cells, "Sample_snSeq.processed"]
) %>%
  dplyr::filter(Exp > 0)

test_data$IsEpilepsy <- factor(test_data$IsEpilepsy, levels = c("Control", "Epilepsy"))

models <- list(
  "Model 1\nCrude" = Exp ~ IsEpilepsy,
  "Model 2\n+ Cohort" = Exp ~ IsEpilepsy + Cohort,
  "Model 3\nCohort + Platform" = Exp ~ IsEpilepsy + Cohort + Platform,
  "Model 4\nCohort + Platform\n+ Age + Gender" = Exp ~ IsEpilepsy + Age + Gender + Platform + Cohort
)

cat("Fitting progressively adjusted linear models...\n")
results <- lapply(names(models), function(m_name) {
  fit <- lm(models[[m_name]], data = test_data)
  s <- summary(fit)$coefficients
  ci <- confint(fit)
  target_row <- "IsEpilepsyEpilepsy"
  
  data.frame(
    Model = m_name,
    Estimate = s[target_row, "Estimate"],
    StdError = s[target_row, "Std. Error"],
    PValue = s[target_row, "Pr(>|t|)"],
    LCI = ci[target_row, 1],
    UCI = ci[target_row, 2]
  )
}) %>% bind_rows()

# Export SourceData
cat("Exporting SourceData...\n")
write.csv(
  results,
  file = file.path(out_data_dir, "FigS2C_ADORA1_Models_Effect_Size_CI.csv"),
  row.names = FALSE
)

cat("Plotting Check Figures...\n")
results$Model <- factor(results$Model, levels = rev(names(models)))

p_forest <- ggplot(results, aes(x = Estimate, y = Model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2, color = "#2c7fb8", linewidth = 0.8) +
  geom_point(size = 4, color = "#2c7fb8") +
  geom_text(
    aes(label = sprintf("Est: %.2f\nP=%.2f", Estimate, PValue)),
    vjust = -0.8, size = 3.3, lineheight = 0.9
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.title.x = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  labs(
    title = "ADORA1 Downregulation in Glu.N",
    x = "Effect Size (Epilepsy - Control)",
    y = NULL
  ) +
  expand_limits(x = c(-0.08, 0.01), y = length(models) + 0.6)

ggsave(
  filename = file.path(out_check_dir, "check_FigS2C_ADORA1_Forest.png"),
  plot = p_forest,
  width = 5.5,
  height = 5,
  dpi = 300
)

cat("FigS2C completed successfully.\n")
