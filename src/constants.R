PNG_WIDTH <- 800
PNG_HEIGHT <- 550

########################### DATA PATHS ###########################
mRNAseq_data_matrix_path <- './data/mRNAseq_LIHC/matrix_RNAseq_LIHC.txt'
miRNAseq_data_matrix_path <- './data/miRNAseq_LIHC/matrix_lihc_miRNASeq.txt'

matrix_DEG_path <- 'DEG/results/matrix_DEG_FC_filtered.txt'
matrix_DEG_mRNA_cancer_path <- 'DEG/results/cancer_matrix_DEG.txt'

clinical_data_matrix_path <- './data/clinical_LIHC/clinical_lihc.txt'
clinical_data_integrated_matrix_path <- './data/clinical_LIHC/clinical_lihc_integration.txt'
##################################################################

PZ_TCGA_REGEX <- "TCGA-\\w+-\\w+"
PZ_SAMPLE_TCGA_REGEX <- "TCGA-\\w+-\\w+-\\w+"

# Reg-ex for patient types (Normal and Cancer)
NORMAL_PZ_REGEX <- "TCGA-\\w+-\\w+-1\\d"
CANCER_PZ_REGEX <- "TCGA-\\w+-\\w+-0\\d"

# Set TRUE to re-run the SNF parameter grid search.
# When FALSE the pipeline uses the hardcoded optimal values in main.R
# (K=10, sigma=1, T=20, n_clusters=6) which are the output of a previous
# grid search run (see plots/snf_top_parameters.csv).
ENABLE_GRID_SEARCH <- FALSE

# Master seed for reproducibility. Re-applied before each stochastic block.
SEED <- 123

dirPlots <- './plots/'

if (!dir.exists(dirPlots)) {
  dir.create(dirPlots)
}

dirDEG <- './DEG/'

if (!dir.exists(dirDEG)) {
  dir.create(dirDEG)
}

dirResults <- paste0(dirDEG, 'results/')

if (!dir.exists(dirResults)) {
  dir.create(dirResults)
}

dirPlotDEG <- './DEG/plot/'

if (!dir.exists(dirPlotDEG)) {
  dir.create(dirPlotDEG)
}





