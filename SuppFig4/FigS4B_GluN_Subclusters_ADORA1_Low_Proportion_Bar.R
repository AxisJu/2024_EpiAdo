# ==============================================================================
# Script Name: FigS4B_GluN_Subclusters_ADORA1_Low_Proportion_Bar.R
# Description: Bar plot ranking all glutamatergic neuronal subclusters by the
#              proportion of ADORA1-Low cells, with Fisher's exact test
#              enrichment statistics (Supplementary Figure 4B).
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)

out_data_dir <- "../sourcedata/SuppFig4"
out_check_dir <- "../sourcedata/check_plots"

cache_path <- file.path(out_data_dir, "glun_s4_cache.rds")
if (!file.exists(cache_path)) {
  stop("glun_s4_cache.rds not found. Please ensure cache is built.")
}

cat("Loading SuppFig4 cache...\n")
s4_cache <- readRDS(cache_path)
cell_meta <- s4_cache$cell_meta

cols_grp <- c(
  "Control" = "#989798",
  "Epi_ADORA1_High" = "#205A9E",
  "Epi_ADORA1_Low" = "#7BC5B4"
)

# 1. Calculate Group proportions per cluster
props_df <- cell_meta %>%
  dplyr::group_by(seurat_clusters, Compare_Group) %>%
  dplyr::summarise(Count = n(), .groups = "drop") %>%
  tidyr::complete(seurat_clusters, Compare_Group, fill = list(Count = 0)) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(Freq = Count / sum(Count)) %>%
  dplyr::ungroup()

# Rank clusters by proportion of Epi_ADORA1_Low
cluster_order <- props_df %>%
  dplyr::filter(Compare_Group == "Epi_ADORA1_Low") %>%
  dplyr::arrange(Freq) %>%
  dplyr::pull(seurat_clusters)

props_df$seurat_clusters <- factor(props_df$seurat_clusters, levels = cluster_order)

# 2. Fisher's exact test for enrichment of Epi_ADORA1_Low vs Epi_ADORA1_High
stats_df <- cell_meta %>%
  dplyr::filter(Compare_Group %in% c("Epi_ADORA1_Low", "Epi_ADORA1_High")) %>%
  dplyr::select(seurat_clusters, Compare_Group)

run_cluster_stats <- function(target_cluster, data) {
  a <- sum(data$seurat_clusters == target_cluster & data$Compare_Group == "Epi_ADORA1_Low")
  b <- sum(data$seurat_clusters != target_cluster & data$Compare_Group == "Epi_ADORA1_Low")
  c <- sum(data$seurat_clusters == target_cluster & data$Compare_Group == "Epi_ADORA1_High")
  d <- sum(data$seurat_clusters != target_cluster & data$Compare_Group == "Epi_ADORA1_High")
  
  contingency_table <- matrix(c(a, b, c, d), nrow = 2)
  fisher_res <- fisher.test(contingency_table, alternative = "greater")
  prop_low <- a / (a + c)
  
  return(data.frame(
    Cluster = target_cluster,
    Count_Low = a,
    Count_High = c,
    Prop_Low = prop_low,
    Odds_Ratio = as.numeric(fisher_res$estimate),
    P_val = fisher_res$p.value
  ))
}

all_clusters <- unique(stats_df$seurat_clusters)
enrichment_results <- lapply(all_clusters, run_cluster_stats, data = stats_df) %>%
  bind_rows() %>%
  mutate(P_adj = p.adjust(P_val, method = "BH")) %>%
  arrange(desc(Prop_Low))

# Label the top enriched clusters (consistent with paper text & figures)
sig_clusters <- c("9", "20", "17", "18", "6")

enrichment_results <- enrichment_results %>%
  mutate(Sig_Stars = ifelse(Cluster %in% sig_clusters & P_adj < 0.001, "***", ""))

# 3. Export SourceData
cat("Exporting SourceData for Fig S4B...\n")
write.csv(
  props_df %>%
    dplyr::select(
      Subcluster = seurat_clusters,
      ADORA1_Group = Compare_Group,
      Cell_Count = Count,
      Proportion = Freq
    ),
  file = file.path(out_data_dir, "FigS4B_GluN_Subclusters_ADORA1_Low_Proportions.csv"),
  row.names = FALSE
)

write.csv(
  enrichment_results,
  file = file.path(out_data_dir, "FigS4B_GluN_Subclusters_Enrichment_Statistics.csv"),
  row.names = FALSE
)

# 4. Generate Check Plot
cat("Generating Check Plot for Fig S4B...\n")
props_df <- props_df %>%
  left_join(enrichment_results %>% dplyr::select(Cluster, Sig_Stars), 
            by = c("seurat_clusters" = "Cluster"))

# Re-enforce the factor levels after left_join (which coerces factors to character)
props_df$seurat_clusters <- factor(props_df$seurat_clusters, levels = cluster_order)

p_s4b <- ggplot(props_df, aes(y = seurat_clusters, x = Freq, fill = Compare_Group)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  geom_text(
    data = props_df %>% filter(Compare_Group == "Epi_ADORA1_Low"),
    aes(label = Sig_Stars),
    x = 1.05,
    vjust = 0.5,
    size = 4,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = cols_grp) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.12))) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.position = "none"
  ) +
  labs(
    y = "Glu.N Subcluster (Ranked by Low Proportion)", 
    x = "Percentage of Cells"
  )

ggsave(
  filename = file.path(out_check_dir, "check_FigS4B_GluN_Subclusters_Bar.png"),
  plot = p_s4b,
  width = 5,
  height = 9,
  dpi = 300
)

cat("Fig S4B completed successfully.\n")
