# 🧠 2024_EpiAdo: Glutamatergic Neuron–Associated Adenosine Deficits in Drug-Resistant Epilepsy

**Official Bioinformatics Codebase for the EpiAdo Study**

[![R 4.2+](https://img.shields.io/badge/R-4.2+-blue.svg)](https://www.r-project.org/)
[![Science Translational Medicine](https://img.shields.io/badge/Science%20Translational%20Medicine-AAAS-red.svg)](https://www.science.org/journal/stm)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status: Under Review](https://img.shields.io/badge/Status-Under%20Review-orange.svg)](#-citation)
[![Seurat](https://img.shields.io/badge/Seurat-5.x-lightgrey.svg)](https://satijalab.org/seurat/)

---

## 📖 Paper Information & Citation

This repository contains the complete analytical pipeline for:

> **Song, Kun†, Xiaoshuai Ji†, Sihan Ju†, Zehan Wu, Lingzhao Min, Fangzhou Li, Ming Chen, Shuhao Mei, Liyi Qian, Jieming Li, Yuhao Xu, Shasha Yang, Ziyang Lin, Yijun Huang, Siheng Feng, Xiang Zou\*, Liang Chen\*, and Ying Mao\*.**
> *"Glutamatergic neuron–associated adenosine deficits define a therapeutically addressable epileptic state."*
> ***Science Translational Medicine*** (Under Review).
>
> †Equal contributions. \*Corresponding authors.

**Affiliations:** Department of Neurosurgery, Huashan Hospital, Fudan University; National Center for Neurological Disorders; State Key Laboratory of Medical Neurobiology and MOE Frontiers Center for Brain Science, Fudan University, Shanghai, China.

```bibtex
@article{song2024epiadо,
  title={Glutamatergic neuron–associated adenosine deficits define a therapeutically addressable epileptic state},
  author={Song, Kun and Ji, Xiaoshuai and Ju, Sihan and Wu, Zehan and Min, Lingzhao and Li, Fangzhou and Chen, Ming and Mei, Shuhao and Qian, Liyi and Li, Jieming and others},
  journal={Science Translational Medicine},
  year={2026},
  note={Under review}
}
```

---

## 🌟 Study Overview

Drug-resistant epilepsy (DRE) affects ~30% of all epilepsy patients and remains a major therapeutic challenge. This study presents a comprehensive multi-omics framework to dissect the cellular and molecular basis of DRE, identifying a conserved adenosine-related deficit in glutamatergic neurons as a key pathological feature, and demonstrating its engagement by deep brain stimulation (DBS) of the anterior thalamic nucleus (ANT).

### Key Findings

1. **Multi-omic landscape of human epileptogenic tissue** (Figure 1):
   - Single-nucleus RNA-seq atlas of **421,714 cells** from multi-region human epileptogenic cortex and hippocampus.
   - Metabolic reprogramming detected across cell types; adenosine metabolism specifically dysregulated in glutamatergic neurons (Glu.N).

2. **Conserved *ADORA1* downregulation in excitatory neurons** (Figure 2):
   - *ADORA1* (adenosine A1 receptor) is consistently downregulated in epileptic Glu.N across major DRE subtypes (FCD, MTLE, LEAT) and multiple cohorts.
   - *ADORA1*^Low Glu.N are enriched for transcriptional signatures of metabolic stress and synaptic dysfunction.

3. **Patch-seq electrophysiology–transcriptomics integration** (Figure 3):
   - Patch-seq profiling of human cortical neurons links *ADORA1* expression to intrinsic hyperexcitability.
   - Reduced *ADORA1* correlates with altered action potential dynamics and seizure burden in patients.

4. **Thalamocortical network involvement** (Figure 7):
   - In a kainate mouse model, adenosine-related abnormalities extend from focal epileptogenic zones to the thalamocortical network.
   - ANT-DBS evokes rapid, neuron-associated adenosine release, acutely suppresses seizure phenotypes, and—with prolonged stimulation—upregulates *Adora1* and reshapes excitatory neuronal transcriptional programs.

---

## 🏗️ Repository Structure

All 55 analysis scripts are organized by publication figure:

```text
2024_EpiAdo/
├── Figure1/                    # Human snRNA-seq atlas & adenosine metabolism landscape
│   ├── Fig1A_snRNA_UMAP_RadialHeatmap.R
│   ├── Fig1B_Marker_BubblePlot.R
│   ├── Fig1C_KEGG_Metabolic_Pathway_Heatmap.R
│   ├── Fig1D_Adenosine_Metabolism_Gene_Heatmap.R
│   ├── Fig1E_GluN_Ado_KEGG_Spearman_Barplot.R
│   └── Fig1F_Purine_Ado_Correlation_Across_CellTypes.R
│
├── Figure2/                    # ADORA1 downregulation and candidate target identification
│   ├── Fig2A_GluN_ADORA1_Expression_Violin.R
│   ├── Fig2B_ADORA1_Density_OddsRatio.R
│   ├── Fig2C_ADORA1_Correlated_Pathways_GSEA.R
│   ├── Fig2D_Candidate_Gene_Heatmap_Correlation_CohenD.R
│   └── Fig2E_Candidate_Gene_Criteria_Ranking.R
│
├── Figure3/                    # Patch-seq electrophysiology-transcriptomics & clinical correlates
│   ├── Fig3B_PatchSeq_UMAP_Mapping.R
│   ├── Fig3C_PatchSeq_Marker_BubblePlot.R
│   ├── Fig3D_GluN_EFeatures_Epilepsy_vs_Control.R
│   ├── Fig3E_PatchSeq_ActionPotential_Traces.R
│   ├── Fig3F_GluN_EFeatures_ADORA1_Low_vs_High.R
│   └── Fig3G_Clinical_SeizureBurden_ADORA1_Correlation.R
│
├── Figure7/                    # Mouse DBS model: adenosine release and transcriptional rescue
│   ├── Fig7G_snRNA_UMAP_Mouse_DBS.R
│   ├── Fig7H_Vglut2_DEG_Volcano_and_GSEA.R
│   ├── Fig7I_Adora1_Expression_Violin_FC.R
│   ├── Fig7J_Ado_Pathways_Expression_Violin.R
│   ├── Fig7K_Adora1_Correlated_Pathways_GSEA.R
│   ├── Fig7L_Dysfunctional_Pathways_Recovery.R
│   └── Fig7M_Ion_Channels_Expression_Violin.R
│
├── SuppFig1/                   # snRNA-seq QC, cell composition, metabolic reprogramming overview
├── SuppFig2/                   # ADORA1 stability across subtypes, LMM, leave-one-out, astrocytes
├── SuppFig3/                   # Stereo-seq spatial transcriptomics QC & landscape
├── SuppFig4/                   # ADORA1^Low Glu.N subclusters and cross-cohort reproducibility
├── SuppFig5/                   # Purine pathway, metabolic-electrophysiology correlations
├── SuppFig6/                   # Candidate marker validation across 18 brain regions
│   └── glun_s6_cache.rds       # Lightweight pre-extracted single-cell cache (1.1 MB)
├── SuppFig7/                   # Patch-seq metadata UMAP & all 12 E-feature projections
│   └── ps_s7_cache.rds         # Lightweight pre-extracted Patch-seq cache (0.02 MB)
├── SuppFig8/                   # Threshold sensitivity, label shuffling, gene permutation, GABA controls
│   └── ps_s8_cache.rds         # Lightweight pre-extracted cache (0.15 MB)
│
├── install_packages.R          # One-command R package installation
├── .gitignore
├── LICENSE
└── README.md
```

---

## 📊 Data Modalities

| Data Type | Source | Cells / Samples | Figure(s) |
| :--- | :--- | :--- | :--- |
| snRNA-seq | Human epileptogenic brain (multi-region, multi-cohort) | 421,714 nuclei | Fig1, Fig2, SuppFig1–6 |
| Stereo-seq (Spatial Transcriptomics) | Human cortical sections | Multiple slides | SuppFig3 |
| Patch-seq (electrophysiology + RNA) | Human cortical neurons | ~380 neurons | Fig3, SuppFig7–8 |
| snRNA-seq | Mouse ANT kainate/DBS model | 31,001 nuclei | Fig7, SuppFig (mouse) |

---

## ⚡ Getting Started

### Requirements
- **R** ≥ 4.2.0 (validated on R 4.5.2)
- **Bioconductor** ≥ 3.17

### Install Dependencies

```R
source("install_packages.R")
```

#### Core Packages

| Category | Packages |
| :--- | :--- |
| Single-cell | `Seurat`, `qs` |
| Differential expression | `DESeq2`, `limma`, `lme4`, `lmerTest` |
| Pathway analysis | `clusterProfiler`, `fgsea`, `enrichplot` |
| Visualization | `ggplot2`, `ggpubr`, `pheatmap`, `ComplexHeatmap`, `cowplot`, `patchwork` |
| Statistics | `effsize`, `rstatix`, `car`, `PMCMRplus` |

---

## 🔬 Reproducing the Figures

Each script is self-contained and can be run from its figure subdirectory. Scripts load data from a local data store (not included due to privacy; see **Data Availability** below), or from lightweight `.rds` caches where available.

```bash
# Example: Reproduce Figure 2B (ADORA1 density & Odds Ratio)
cd Figure2
Rscript Fig2B_ADORA1_Density_OddsRatio.R

# Example: Reproduce SuppFig 8A (threshold sensitivity, uses included cache)
cd SuppFig8
Rscript FigS8A_Threshold_Sensitivity.R
```

Each script:
1. Loads required input data (see `[DATA REQUIRED]` comments for data specification)
2. Exports publication-formatted **Source Data** CSV files to `../sourcedata/`
3. Outputs high-resolution **check plots** to `../sourcedata/check_plots/`

### Cache-Accelerated Scripts

Scripts in `SuppFig5–8` include pre-extracted lightweight `.rds` cache files and can be run immediately **without the raw Seurat objects**, completing in seconds:

| Cache File | Size | Used By |
| :--- | :--- | :--- |
| `SuppFig5/glun_s5_cache.rds` | 1.4 MB | FigS5A–E |
| `SuppFig6/glun_s6_cache.rds` | 1.1 MB | FigS6A |
| `SuppFig7/ps_s7_cache.rds` | 0.02 MB | FigS7A–B |
| `SuppFig8/ps_s8_cache.rds` | 0.15 MB | FigS8A–DE |

---

## 🔒 Data Availability & Privacy

Raw single-nucleus RNA-seq data, Patch-seq recordings, and clinical metadata are subject to patient privacy regulations and are not directly included in this repository. Lines loading patient-level data have been replaced with `# [DATA REQUIRED]` placeholders throughout the scripts.

Upon publication, processed data (including cell-by-gene matrices, UMAP coordinates, and pseudobulk summaries) will be deposited in a public repository (e.g., GEO / Zenodo) in accordance with Science Translational Medicine data deposition policies.

---

## ⚠️ Analysis Notes & Data Filtering Rationale

The following pre-analysis filtering steps are applied and are explicitly commented in the relevant scripts:

| Script | Filter | Rationale |
| :--- | :--- | :--- |
| `Fig2B_ADORA1_Density_OddsRatio.R` | Cohorts `"Tran"`, `"YC"`, `"LWS"` excluded from density plot | These cohorts use different sequencing platforms; CPM-scale expression distributions are not directly comparable. They are included in pseudobulk analyses (Fig2A, Fig2C) with batch correction. |
| `Fig3D_GluN_EFeatures_Epilepsy_vs_Control.R` | Control neurons filtered to Resting Membrane Potential < −50 mV | Standard electrophysiology quality criterion; excludes cells with compromised membrane integrity. Epilepsy neurons are included irrespective of RMP. |
| `Fig3G_Clinical_SeizureBurden_ADORA1_Correlation.R` | Restricted to `Cohort == "HS"` | Only the HS (Huashan Hospital) cohort has complete paired clinical seizure burden metadata; this is a single-cohort analysis by design. |
| All `Fig1A`, `FigS1A` | `celltype_coarse != "LowQuality"` | Standard Seurat doublet/low-quality cell removal step; thresholds defined in Methods. |

---

## 📜 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

---

## 🤝 Acknowledgments

This research was supported by grants from the National Natural Science Foundation of China, the Shanghai Municipal Commission of Health and Family Planning, and the Tianqiao and Chrissy Chen Institute for Neuroscience.

Bioinformatics analysis code was developed with the assistance of **Google Gemini 3.7 Flash** and the **Antigravity** agentic coding system.
