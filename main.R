setwd("./")
rm(list = ls())
options(stringsAsFactors = F)

source("src/utils.R")
source("src/constants.R")
source("src/plot_network.R")
source("src/missing_values_analysis.R")


# Pre-processing
source("src/preprocessing/deg_pipeline.R")
source("src/preprocessing/dem_pipeline.R")
source("src/preprocessing/preprocessing_utils.R")
source("src/preprocessing/clinical_data_matrix.R")

# Network fusion
source("src/grid_search_SNF.R")
source("src/network_analysis_SNF.R")

# Output evaluation
source("src/post_clustering_analysis/survival_analysis.R")
source("src/post_clustering_analysis/snf_enrichment_analysis.R")
source("src/post_clustering_analysis/snf_clinical_association.R")
source("src/post_clustering_analysis/snf_cluster_deg_analysis.R")
source("src/post_clustering_analysis/modified_deg_for_log2_data.R")
source("src/post_clustering_analysis/tcga_lihc_complete_analysis.R")

################################ 1. Input data ################################

# Read mRNA-seq and miRNA-seq
mRNA_seq_matrix <- get_mRNAseq_DEG_dataframe(mRNAseq_data_matrix_path)
miRNA_seq_matrix <- get_miRNAseq_dataframe(miRNAseq_data_matrix_path)

# Filtering to have the same patients
filtered_matrix <- intersect_pz_mRNA_miRNA(mRNA_seq_matrix,
                                           miRNA_seq_matrix)

mRNA_matrix_cancer <- filtered_matrix$mRNA_matrix_final
miRNA_matrix_cancer <- filtered_matrix$miRNA_matrix_final
common_pz_omics <- filtered_matrix$common_pz_cancer

# Ensure clinical matrix matches the cancer patients
cancer_patient_ids <- unique(str_extract(colnames(mRNA_matrix_cancer), PZ_TCGA_REGEX))

final_clinical_matrix <- compute_final_clinical_matrix(common_pz_omics = cancer_patient_ids)

# Verify all matrices have the same patients in the same order
stopifnot(identical(
  str_extract(colnames(mRNA_matrix_cancer), PZ_TCGA_REGEX),
  str_extract(colnames(miRNA_matrix_cancer), PZ_TCGA_REGEX)
))
stopifnot(identical(
  str_extract(colnames(mRNA_matrix_cancer), PZ_TCGA_REGEX),
  rownames(final_clinical_matrix)
))

############################# 2. Network analysis #############################
distance <- "correlation"

net_analysis_results <- network_analysis_SNF(mRNA_matrix_cancer,
                                             miRNA_matrix_cancer,
                                             distance,
                                             K = 10,
                                             sigma = 1,
                                             n_clusters = 6,
                                             iterations = 20,
                                             dir_plots = dirPlots)

W_fused <- net_analysis_results$W_fused
final_clusters <- net_analysis_results$final_clusters
patient_to_cluster <- net_analysis_results$patient_to_cluster

dist_mRNA <- net_analysis_results$dist_mRNA
dist_miRNA <- net_analysis_results$dist_miRNA

########################### Parameter Grid Search for SNF ###########################

# Define parameter ranges
K_range <- c(10, 15, 20, 50, 60, 80, 100, 120)
alpha_range <- seq(0.4, 1.2, by = 0.2)
T_range <- c(20, 25, 30, 50, 100, 200)

grid_search_SNF(dir_plots = dirPlots,
                dist_mRNA,
                dist_miRNA,
                K_range,
                alpha_range,
                T_range)

############################## Post-clustering analysis #############################

################################# Survival analysis #################################

clinical_legacy <- read.table(clinical_data_matrix_path,
                              header = T, sep = "\t", check.names = F,
                              row.names = 1, quote = "")

clinical_legacy <- clear_clinical_data_columns(clinical_legacy)

columns_to_remove <- c(
  'synchronous_malignancy',
  'days_to_diagnosis',
  'primary_diagnosis',
  'year_of_diagnosis',
  'prior_treatment',
  'ajcc_staging_system_edition',
  'morphology',
  'race',
  'ethnicity',
  'days_to_birth',
  'year_of_birth',
  'year_of_death',
  'bcr_patient_barcode')

clinical_legacy <- clinical_legacy[, setdiff(colnames(clinical_legacy), columns_to_remove)]

clinical_updated <- read.table(clinical_data_integrated_matrix_path,
                               header = T, sep = "\t", check.names = F,
                               row.names = 1, quote = "")

clinical_updated <- clear_clinical_data_columns(clinical_updated)

columns_to_remove <- c(
  'asian_race',
  'caucasian_race')

clinical_updated <- clinical_updated[, setdiff(colnames(clinical_updated), columns_to_remove)]

summarize_features(clinical_legacy, "summary_clinical_legacy.csv")
summarize_features(clinical_updated, "summary_clinical_updated.csv")

# Run the survival analysis
survival_results <- run_complete_hcc_suvival_analysis(clinical_legacy,
                                                      clinical_updated,
                                                      patient_to_cluster)

# Run complete final analysis
final_clinical_matrix$snf_cluster <- final_clusters

results <- run_tcga_lihc_snf_analysis(
  mRNA_matrix_cancer = mRNA_matrix_cancer,
  miRNA_matrix_cancer = miRNA_matrix_cancer,
  final_clinical_matrix = final_clinical_matrix,
  output_base_dir = "TCGA_LIHC_SNF_Results_final"
)

#####################################################################################
