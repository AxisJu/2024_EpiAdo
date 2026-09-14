# 2024_EpiAdo

Bioinformatics analysis pipelines and statistical workflows for single-nucleus RNA-seq, Stereo-seq, and Patch-seq profiling of human drug-resistant epilepsy tissue and mouse DBS models.

## Repository Architecture

```text
2024_EpiAdo/
├── Figure1/                    # Human snRNA-seq atlas and adenosine metabolism landscape
│   ├── Fig1A_snRNA_UMAP_RadialHeatmap.R
│   ├── Fig1B_Marker_BubblePlot.R
│   ├── Fig1C_KEGG_Metabolic_Pathway_Heatmap.R
│   ├── Fig1D_Adenosine_Metabolism_Gene_Heatmap.R
│   ├── Fig1E_GluN_Ado_KEGG_Spearman_Barplot.R
│   └── Fig1F_Purine_Ado_Correlation_Across_CellTypes.R
├── Figure2/                    # ADORA1 expression analysis and candidate target screening
│   ├── Fig2A_GluN_ADORA1_Expression_Violin.R
│   ├── Fig2B_ADORA1_Density_OddsRatio.R
│   ├── Fig2C_ADORA1_Correlated_Pathways_GSEA.R
│   ├── Fig2D_Candidate_Gene_Heatmap_Correlation_CohenD.R
│   └── Fig2E_Candidate_Gene_Criteria_Ranking.R
├── Figure3/                    # Patch-seq electrophysiology-transcriptomics and clinical correlations
│   ├── Fig3B_PatchSeq_UMAP_Mapping.R
│   ├── Fig3C_PatchSeq_Marker_BubblePlot.R
│   ├── Fig3D_GluN_EFeatures_Epilepsy_vs_Control.R
│   ├── Fig3E_PatchSeq_ActionPotential_Traces.R
│   ├── Fig3F_GluN_EFeatures_ADORA1_Low_vs_High.R
│   └── Fig3G_Clinical_SeizureBurden_ADORA1_Correlation.R
├── Figure7/                    # Mouse DBS model: single-nucleus RNA-seq and transcriptional rescue
│   ├── Fig7G_snRNA_UMAP_Mouse_DBS.R
│   ├── Fig7H_Vglut2_DEG_Volcano_and_GSEA.R
│   ├── Fig7I_Adora1_Expression_Violin_FC.R
│   ├── Fig7J_Ado_Pathways_Expression_Violin.R
│   ├── Fig7K_Adora1_Correlated_Pathways_GSEA.R
│   ├── Fig7L_Dysfunctional_Pathways_Recovery.R
│   └── Fig7M_Ion_Channels_Expression_Violin.R
├── SuppFig1/                   # snRNA-seq QC, cell composition, and metabolic reprogramming overview
├── SuppFig2/                   # ADORA1 stability across neuropathological subtypes and linear mixed models
├── SuppFig3/                   # Stereo-seq spatial transcriptomics QC and cell-type mapping
├── SuppFig4/                   # Excitatory neuron subclusters and cross-cohort validation
├── SuppFig5/                   # Purine pathways and electrophysiological gene correlations
├── SuppFig6/                   # Candidate marker expression across anatomical brain regions
├── SuppFig7/                   # Patch-seq clinical metadata and electrophysiology feature projections
├── SuppFig8/                   # Threshold sensitivity, label shuffling, and GABAergic control analyses
├── install_packages.R          # Automated R package installer
├── .gitignore
└── LICENSE
```

## Prerequisites

- R >= 4.2.0
- Bioconductor >= 3.17

Install required CRAN and Bioconductor packages by running:

```R
source("install_packages.R")
```

Key dependencies:
- Single-cell analysis: `Seurat`, `qs`
- Differential expression: `DESeq2`, `limma`, `lme4`, `lmerTest`
- Pathway and functional enrichment: `clusterProfiler`, `fgsea`, `enrichplot`
- Visualization and statistics: `ggplot2`, `ggpubr`, `ComplexHeatmap`, `pheatmap`, `cowplot`, `patchwork`, `effsize`, `rstatix`

## Publication

Kun Song, Xiaoshuai Ji, Sihan Ju, Zehan Wu, Lingzhao Min, Fangzhou Li, Ming Chen, Shuhao Mei, Liyi Qian, Jieming Li, Yuhao Xu, Shasha Yang, Ziyang Lin, Yijun Huang, Siheng Feng, Xiang Zou, Liang Chen, and Ying Mao. "Glutamatergic neuron–associated adenosine deficits define a therapeutically addressable epileptic state." *Science Translational Medicine* (Under Review).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
