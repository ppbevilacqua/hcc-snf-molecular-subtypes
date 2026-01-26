# Helper function to filter columns by patient IDs
filter_columns_by_patients <- function(matrix, patient_ids) {

  # Extract patient IDs from column names
  col_patient_ids <- str_extract(colnames(matrix), PZ_SAMPLE_TCGA_REGEX)

  # Find columns that match the desired patient IDs
  matching_cols <- which(col_patient_ids %in% patient_ids)

  # Return filtered matrix
  return(matrix[, matching_cols, drop = FALSE])
}

intersect_pz_mRNA_miRNA <- function(mRNA_seq_matrix,
                                    miRNA_seq_matrix) {

  # Get common patients from mRNA and miRNA
  pz_mRNA <- unique(str_extract(colnames(mRNA_seq_matrix), PZ_SAMPLE_TCGA_REGEX))

  pz_miRNA <- unique(str_extract(colnames(miRNA_seq_matrix), PZ_SAMPLE_TCGA_REGEX))

  common_pz <- intersect(pz_mRNA, pz_miRNA)

  # Filter columns based on common patients
  mRNA_matrix <- filter_columns_by_patients(mRNA_seq_matrix, common_pz)
  miRNA_matrix <- filter_columns_by_patients(miRNA_seq_matrix, common_pz)

  # Edit patient names to use only TCGA barcode (removing sample type suffix)
  colnames(mRNA_matrix) <- str_extract(colnames(mRNA_matrix), PZ_TCGA_REGEX)
  colnames(miRNA_matrix) <- str_extract(colnames(miRNA_matrix), PZ_TCGA_REGEX)

  # Order columns by patient ID for consistent ordering
  mRNA_matrix <- mRNA_matrix[, order(colnames(mRNA_matrix)), drop = FALSE]
  miRNA_matrix <- miRNA_matrix[, order(colnames(miRNA_matrix)), drop = FALSE]

  return(list(
    mRNA_matrix_final = mRNA_matrix,
    miRNA_matrix_final = miRNA_matrix,
    common_pz_cancer = colnames(mRNA_matrix)
  ))
}
