library(stringr)
library(pheatmap)

# Network and plot
library(caret)
library(igraph)
library(mclust)
library(cluster)
library(ggrepel)
library(SNFtool)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)

### K-Nearest Neighbors Imputation
library(VIM)

### Imputation quality
library(Hmisc)

# Missing values analysis
library(dplyr)
library(reshape2)
library(tidyr)

# Survival analysis
library(survival)
library(survminer)
library(tableone)
library(forestplot)
library(ComplexHeatmap)

# Ensure dplyr functions take precedence
select <- dplyr::select
filter <- dplyr::filter
arrange <- dplyr::arrange

# Cluster analysis and enrichment analysis
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("limma")
  BiocManager::install("edgeR")
  BiocManager::install("org.Hs.eg.db")
  BiocManager::install("clusterProfiler")
}

suppressPackageStartupMessages({
  library(limma)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggpubr)
  library(ggsignif)
  library(corrplot)
})

# Set theme for plots
theme_set(theme_bw())

# NOTE: a global set.seed() used to live here. Seeding is now applied
# locally before each stochastic block (spectral clustering loops, grid
# search, cross-validation) using the SEED constant from constants.R.


# Filter to get only cancer samples
get_cancer_samples <- function(data_matrix) {
  cancer_cols <- grep(CANCER_PZ_REGEX, colnames(data_matrix), value = TRUE)
  cancer_matrix <- data_matrix[, cancer_cols]

  # Extract unique patient IDs
  patient_ids <- unique(str_extract(colnames(cancer_matrix), PZ_TCGA_REGEX))

  # For patients with multiple cancer samples, take the first one
  final_cols <- sapply(patient_ids, function(pid) {
    patient_samples <- grep(pid, colnames(cancer_matrix), value = TRUE)
    patient_samples[1]  # Take first sample if multiple exist
  })

  return(cancer_matrix[, final_cols])
}

get_var_label <- function(var) {
  var_labels <- c(
    "age_at_index" = "Age at Diagnosis",
    "gender" = "Gender",
    "hepatitis_b" = "Hepatitis B Status",
    "ajcc_pathologic_stage" = "AJCC Pathologic Stage",
    "ajcc_pathologic_t" = "AJCC Pathologic T",
    "ajcc_pathologic_m" = "AJCC Pathologic M"
  )
  
  ifelse(var %in% names(var_labels), var_labels[var], var)
}

