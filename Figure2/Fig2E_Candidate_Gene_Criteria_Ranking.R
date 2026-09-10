# ==============================================================================
# Script: Fig2E_Candidate_Gene_Criteria_Ranking.R
# Purpose: Dot plot ranking candidate genes based on the number of satisfied criteria (0 to 12)
# Article: EpiAdo (Science Translational Medicine) - Figure 2E
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

out_dir_sd <- "../sourcedata/Figure2"
out_dir_plot <- "../sourcedata/check_plots"
dir.create(out_dir_sd, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_plot, showWarnings = FALSE, recursive = TRUE)

de_dir <- "Z:/2023_EpiAdo/results/20260105"
regions <- c("Cortex", "Hippocampus", "Amygdala")
comparisons <- c("Low_vs_Ctrl" = "Control", "Low_vs_High" = "Epi_ADORA1_High")
target_genes <- c("RTN4", "GCSH", "BMPER")

cat("Loading glu_n_seurat (3.3 GB)...\n")
# [DATA REQUIRED] The line below originally loaded '20260309_glu_n_seurat.qs' from a private path.
# Please provide 'glu_n_seurat' by loading the appropriate data before running this script.
# Example: glu_n_seurat <- qread('path/to/20260309_glu_n_seurat.qs')
all_genes <- rownames(glu_n_seurat)

# Fast vectorized Cohen's d functions
calc_cohen_d_sparse <- function(mat1, mat2) {
  n1 <- ncol(mat1); n2 <- ncol(mat2)
  if (n1 < 2 | n2 < 2) return(rep(NA, nrow(mat1)))
  m1 <- Matrix::rowMeans(mat1); m2 <- Matrix::rowMeans(mat2)
  v1 <- (Matrix::rowSums(mat1 * mat1) - n1 * m1^2) / (n1 - 1)
  v2 <- (Matrix::rowSums(mat2 * mat2) - n2 * m2^2) / (n2 - 1)
  v1[v1 < 0] <- 0; v2[v2 < 0] <- 0
  sd_pool <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
  d <- (m1 - m2) / sd_pool
  d[sd_pool == 0 | is.na(sd_pool)] <- 0
  return(d)
}

calc_cohen_d_dense <- function(mat1, mat2) {
  n1 <- ncol(mat1); n2 <- ncol(mat2)
  if (n1 < 2 | n2 < 2) return(rep(NA, nrow(mat1)))
  m1 <- rowMeans(mat1, na.rm = TRUE); m2 <- rowMeans(mat2, na.rm = TRUE)
  v1 <- apply(mat1, 1, var, na.rm = TRUE); v2 <- apply(mat2, 1, var, na.rm = TRUE)
  sd_pool <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
  d <- (m1 - m2) / sd_pool
  d[sd_pool == 0 | is.na(sd_pool)] <- 0
  return(d)
}

get_pb_mat <- function(mat, samples) {
  uniq_s <- unique(samples)
  res <- sapply(uniq_s, function(s) {
    idx <- which(samples == s)
    if (length(idx) == 1) return(mat[, idx])
    Matrix::rowMeans(mat[, idx, drop = FALSE])
  })
  rownames(res) <- rownames(mat)
  colnames(res) <- uniq_s
  return(res)
}

cat("Evaluating 12 criteria across brain regions and comparison groups...\n")
score_mat <- matrix(FALSE, nrow = length(all_genes), ncol = 12)
cohen_mat <- matrix(0, nrow = length(all_genes), ncol = 12)
rownames(score_mat) <- all_genes
rownames(cohen_mat) <- all_genes
col_names <- c()

pb_meta <- glu_n_seurat@meta.data %>%
  dplyr::filter(Compare_Group %in% c("Control", "Epi_ADORA1_Low", "Epi_ADORA1_High"))
expr_data <- GetAssayData(glu_n_seurat, slot = "data")

col_idx <- 1
for (reg in regions) {
  for (comp_name in names(comparisons)) {
    ref_group <- comparisons[[comp_name]]
    sc_col <- paste0("SC_", reg, "_", comp_name)
    pb_col <- paste0("PB_", reg, "_", comp_name)
    col_names <- c(col_names, sc_col, pb_col)
    
    cells_low <- rownames(pb_meta)[pb_meta$BrainRegion == reg & pb_meta$Compare_Group == "Epi_ADORA1_Low"]
    cells_ref <- rownames(pb_meta)[pb_meta$BrainRegion == reg & pb_meta$Compare_Group == ref_group]
    
    mat_low <- expr_data[, cells_low, drop = FALSE]
    mat_ref <- expr_data[, cells_ref, drop = FALSE]
    
    # 1. SC level
    fname <- file.path(de_dir, paste0("7.2.2_DE_", reg, "_", comp_name, ".csv"))
    if (file.exists(fname)) {
      sc_res <- read.csv(fname, row.names = 1)
      sig_genes <- rownames(sc_res)[sc_res$p_val_adj < 0.05 & abs(sc_res$avg_log2FC) > 0.1]
      score_mat[sig_genes, col_idx] <- TRUE
    }
    cohen_mat[, col_idx] <- calc_cohen_d_sparse(mat_low, mat_ref)
    col_idx <- col_idx + 1
    
    # 2. PB level
    samples_low <- pb_meta[cells_low, "Sample_snSeq.processed"]
    samples_ref <- pb_meta[cells_ref, "Sample_snSeq.processed"]
    pb_mat_low <- get_pb_mat(mat_low, samples_low)
    pb_mat_ref <- get_pb_mat(mat_ref, samples_ref)
    
    n_low <- ncol(pb_mat_low)
    n_ref <- ncol(pb_mat_ref)
    if (n_low >= 3 & n_ref >= 3) {
      pb_pvals <- apply(cbind(pb_mat_low, pb_mat_ref), 1, function(x) {
        val_low <- x[1:n_low]
        val_ref <- x[(n_low + 1):(n_low + n_ref)]
        if (sum(!is.na(val_low)) >= 3 & sum(!is.na(val_ref)) >= 3) {
          return(suppressWarnings(wilcox.test(val_low, val_ref, exact = FALSE)$p.value))
        } else {
          return(NA)
        }
      })
      pb_sig <- names(which(pb_pvals < 0.05))
      score_mat[pb_sig, col_idx] <- TRUE
    }
    cohen_mat[, col_idx] <- calc_cohen_d_dense(pb_mat_low, pb_mat_ref)
    col_idx <- col_idx + 1
  }
}
colnames(score_mat) <- col_names
colnames(cohen_mat) <- col_names

# Summarize metrics
gene_scores <- rowSums(score_mat, na.rm = TRUE)
sum_cohen <- rowSums(cohen_mat, na.rm = TRUE)
abs_sum_cohen <- rowSums(abs(cohen_mat), na.rm = TRUE)

master_df <- data.frame(
  Gene = all_genes,
  Score = gene_scores,
  Total_Cohen = sum_cohen,
  Abs_Total_Cohen = abs_sum_cohen
) %>%
  dplyr::mutate(
    Label = ifelse(Gene %in% target_genes, Gene, NA),
    Highlight = ifelse(Gene %in% target_genes, "Target", "Other")
  )

df_plot1 <- master_df %>%
  dplyr::arrange(Score, Abs_Total_Cohen) %>%
  dplyr::mutate(Rank = dplyr::row_number())

# 1. Export SourceData
cat("Exporting SourceData for Figure 2E...\n")
sd_ranking <- df_plot1 %>%
  dplyr::select(
    Rank,
    Gene,
    Number_of_Satisfied_Criteria = Score,
    Sum_Abs_CohenD = Abs_Total_Cohen,
    Total_CohenD = Total_Cohen,
    Is_Target = Highlight
  )
write.csv(sd_ranking, file.path(out_dir_sd, "Fig2E_Candidate_Gene_Criteria_Ranking.csv"), row.names = FALSE)

# 2. Generate Check Plot
cat("Generating Check Plot for Figure 2E...\n")
p1 <- ggplot(df_plot1, aes(x = Rank, y = Score)) +
  geom_point(
    data = filter(df_plot1, Highlight == "Other"),
    aes(size = Abs_Total_Cohen, color = Highlight, fill = Highlight, alpha = Highlight),
    shape = 21
  ) +
  geom_point(
    data = filter(df_plot1, Highlight == "Target"),
    aes(size = Abs_Total_Cohen, color = Highlight, fill = Highlight, alpha = Highlight),
    shape = 21
  ) +
  scale_fill_manual(values = c("Target" = "#D6604D", "Other" = "transparent")) +
  scale_color_manual(values = c("Target" = "black", "Other" = "grey85")) +
  scale_alpha_manual(values = c("Target" = 1, "Other" = 0.4)) +
  scale_size_continuous(range = c(0.5, 5), name = "Sum of |Cohen's D|") +
  geom_text_repel(
    data = filter(df_plot1, Highlight == "Target"),
    aes(label = Label), color = "black", size = 4.5, fontface = "bold",
    box.padding = 0.8, point.padding = 0.5, max.overlaps = Inf,
    segment.color = "black", segment.size = 0.5
  ) +
  scale_y_continuous(breaks = seq(0, 12, 2), limits = c(0, 12.5)) +
  theme_classic() +
  labs(x = "Gene Rank", y = "Number of Satisfied Criteria") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(out_dir_plot, "check_Fig2E_Candidate_Gene_Criteria_Ranking.png"), plot = p1, width = 5, height = 7, dpi = 300)
cat("Figure 2E finished successfully.\n")
