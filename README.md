# Multi-Omics Integration for Hepatocellular Carcinoma Patient Stratification

A graph-based bioinformatics framework for integrating multi-omics data (mRNA-seq and miRNA-seq) to identify molecularly distinct patient subgroups in hepatocellular carcinoma (HCC) using Similarity Network Fusion (SNF).

## Overview

This project implements a complete bioinformatics pipeline for patient stratification in liver hepatocellular carcinoma using data from The Cancer Genome Atlas (TCGA-LIHC). The framework integrates multiple omics layers through similarity network fusion, constructs a patient graph as an integrated representation of molecular relationships, and identifies clinically relevant subgroups through spectral clustering.

### Key Findings

- **366 patients** from TCGA-LIHC were stratified into **6 molecular subtypes** with distinct prognoses
- **Cluster 1** emerges as the most aggressive subtype with poorest survival outcomes
- **Cluster 4** shows the most favorable prognosis
- Cluster membership provides **independent prognostic value** even after adjusting for age, gender, and tumor stage (multivariate Cox regression)
- The molecular stratification captures biological information beyond conventional clinical variables

## Pipeline Overview

```
Raw Data (mRNA + miRNA + Clinical)
                ↓
┌───────────────────────────────────────┐
│         PREPROCESSING                 │
├───────────────────────────────────────┤
│ • Log2 transformation                 │
│ • IQR-based feature selection         │
│ • Differential expression analysis    │
│ • Sample intersection alignment       │
└───────────────────────────────────────┘
                ↓
┌───────────────────────────────────────┐
│      NETWORK ANALYSIS (SNF)           │
├───────────────────────────────────────┤
│ • Spearman correlation distance       │
│ • K-NN affinity matrices              │
│ • Multi-modal network fusion          │
│ • Spectral clustering                 │
│ • Parameter optimization              │
└───────────────────────────────────────┘
                ↓
┌───────────────────────────────────────┐
│     POST-CLUSTERING ANALYSIS          │
├───────────────────────────────────────┤
│ • Kaplan-Meier survival analysis      │
│ • Cox proportional hazards regression │
│ • Gene ontology enrichment            │
│ • KEGG pathway analysis               │
│ • Clinical association tests          │
└───────────────────────────────────────┘
                ↓
    Patient Clusters + Biomarkers + Survival Signatures
```

## Project Structure

```
thesis-project/
├── main.R                              # Main pipeline orchestration
├── data/
│   ├── mRNAseq_LIHC/                   # mRNA expression matrices
│   ├── miRNAseq_LIHC/                  # miRNA expression matrices
│   └── clinical_LIHC/                  # Clinical metadata
├── src/
│   ├── constants.R                     # Global constants and paths
│   ├── utils.R                         # Utility functions
│   ├── missing_values_analysis.R       # Missing data visualization
│   ├── plot_network.R                  # Network visualization
│   ├── network_analysis_SNF.R          # SNF implementation
│   ├── grid_search_SNF.R               # Parameter optimization
│   ├── preprocessing/
│   │   ├── deg_pipeline.R              # Differential gene expression
│   │   ├── dem_pipeline.R              # Differential miRNA expression
│   │   ├── preprocessing_utils.R       # Data filtering utilities
│   │   └── clinical_data_matrix.R      # Clinical data processing
│   └── post_clustering_analysis/
│       ├── survival_analysis.R         # Survival analysis
│       ├── snf_enrichment_analysis.R   # GO/KEGG enrichment
│       ├── snf_cluster_deg_analysis.R  # Cluster-specific DEGs
│       └── snf_clinical_association.R  # Clinical associations
├── DEG/                                # Differential expression results
│   ├── results/
│   └── plot/
└── plots/                              # Output visualizations
```

## Requirements

### R Dependencies

**Core Data Manipulation:**
- dplyr, tidyverse, reshape2, tidyr

**Statistical Analysis:**
- limma, edgeR, survival

**Bioinformatics:**
- SNFtool (Similarity Network Fusion)
- clusterProfiler (Enrichment analysis)
- org.Hs.eg.db (Human gene annotation)
- ComplexHeatmap

**Network & Clustering:**
- igraph
- mclust
- cluster
- caret

**Visualization:**
- ggplot2, pheatmap, ggrepel
- RColorBrewer, gridExtra
- survminer, forestplot
- ggpubr, ggsignif, corrplot

**Missing Value Analysis:**
- VIM, Hmisc

**Clinical Analysis:**
- tableone

### Installation

Install required R packages:

```r
# CRAN packages
install.packages(c("dplyr", "tidyverse", "reshape2", "tidyr", "ggplot2",
                   "pheatmap", "ggrepel", "RColorBrewer", "gridExtra",
                   "igraph", "mclust", "cluster", "caret", "survival",
                   "survminer", "forestplot", "ggpubr", "ggsignif",
                   "corrplot", "VIM", "Hmisc", "tableone", "stringr"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "edgeR", "SNFtool", "clusterProfiler",
                       "org.Hs.eg.db", "ComplexHeatmap"))
```

## Data

The input data required to reproduce the analyses is available as a [release archive (v1-data)](https://github.com/ppbevilacqua/hcc-snf-molecular-subtypes/releases/tag/v1-data). Download `data.zip` and extract it into the project root so that the `data/` directory is created with the expected structure. The analysis scripts assume all data files are available in this location and can be executed without further configuration.

The pipeline uses TCGA-LIHC (Liver Hepatocellular Carcinoma) data:

| Data Type | Description | Samples |
|-----------|-------------|---------|
| mRNA-seq | Gene expression (60,660 genes) | 424 |
| miRNA-seq | miRNA expression | 424 |
| Clinical | Patient metadata (80+ variables) | 377 |

**Final cohort:** 366 patients with complete multi-omics data

## Usage

### Running the Complete Pipeline

```r
source("main.R")
```

### Step-by-Step Execution

1. **Data Preprocessing**
```r
source("src/utils.R")
source("src/constants.R")
source("src/preprocessing/deg_pipeline.R")
source("src/preprocessing/dem_pipeline.R")
```

2. **Network Analysis**
```r
source("src/network_analysis_SNF.R")
# Parameters: K=10, sigma=1, iterations=20, clusters=6
```

3. **Parameter Optimization**
```r
source("src/grid_search_SNF.R")
# Tests: K=[10-120], sigma=[0.4-1.2], T=[20-200]
```

4. **Post-Clustering Analysis**
```r
source("src/post_clustering_analysis/survival_analysis.R")
source("src/post_clustering_analysis/snf_enrichment_analysis.R")
source("src/post_clustering_analysis/snf_clinical_association.R")
```

## Methods

### Differential Expression Analysis
- Zero expression removal and log2 transformation
- IQR-based feature selection (64th percentile threshold)
- Fold-change filtering (|log2FC| > 1)
- FDR-adjusted p-value < 0.01

### Similarity Network Fusion
- **Distance metric:** 1 - Spearman correlation
- **Affinity matrix:** K-nearest neighbors approach
- **Fusion:** Iterative network fusion of mRNA and miRNA networks
- **Clustering:** Spectral clustering with eigen-gap optimization

### Cluster Validation
- Silhouette analysis
- Within-cluster similarity
- Modularity metrics

### Clinical Characterization
- Kaplan-Meier survival curves with log-rank test
- Univariate and multivariate Cox proportional hazards regression
- Clinical feature association tests (chi-square, Fisher's exact)

## Outputs

### Visualizations (`plots/`)
- `fused_network_clusters.png` - Network visualization with cluster coloring
- `similarity_heatmap_*.png` - Similarity matrices with cluster annotations
- `silhouette_analysis.png` - Cluster quality assessment
- `consensus_clustering.png` - Clustering stability analysis

### Results (`DEG/results/`)
- `DEG.txt` - Differentially expressed genes (95 genes: 71 up, 24 down)
- `cancer_matrix_DEG.txt` - DEG-filtered expression matrix

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| K | 10 | Number of nearest neighbors |
| sigma | 1 | Affinity kernel width |
| T | 20 | SNF iterations |
| n_clusters | 6 | Number of patient clusters |
| IQR threshold | 0.64 | Feature selection percentile |
| log2FC | 1 | Fold-change cutoff |
| FDR | 0.01 | Adjusted p-value threshold |

## References

This project implements the Similarity Network Fusion (SNF) algorithm. If you use this pipeline, please cite the following:

### SNF Algorithm

> Wang, B., Mezlini, A.M., Demir, F., Fiume, M., Tu, Z., Brudno, M., Haibe-Kains, B., & Goldenberg, A. (2014). **Similarity network fusion for aggregating data types on a genomic scale.** *Nature Methods*, 11(3), 333–337. https://doi.org/10.1038/nmeth.2810

### SNFtool R Package

> Wang, B., Mezlini, A.M., Demir, F., Fiume, M., Tu, Z., Brudno, M., Haibe-Kains, B., & Goldenberg, A. (2021). *SNFtool: Similarity Network Fusion* (Version 2.3.1). CRAN. https://doi.org/10.32614/CRAN.package.SNFtool

## License

This project is developed for academic research purposes. The SNFtool package is licensed under GPL-2 | GPL-3.

## Acknowledgments

- [TCGA Research Network](https://www.cancer.gov/tcga) for providing the LIHC dataset
- Wang et al. for developing the [Similarity Network Fusion](https://www.nature.com/articles/nmeth.2810) algorithm
- The [SNFtool](https://CRAN.R-project.org/package=SNFtool) package maintainers (Benjamin Brew and the Goldenberg Lab)