get_miRNAseq_dataframe <- function(
  data_miRNAseq_matrix_path
) {
  tmp <- read.table(data_miRNAseq_matrix_path, header = T,
                    sep = "\t", check.names = F,
                    row.names = 1, quote = "")

  pz <- colnames(tmp)

  pzC <- grep(CANCER_PZ_REGEX, pz, value = TRUE)

  # Remove duplicates based on patient ID
  pzC <- pzC[!duplicated(str_extract(pzC, PZ_TCGA_REGEX))]

  # Get ALL cancer samples for final output
  dataC_all <- tmp[, pzC]

  # Remove zero means data using paired samples for consistency
  overall_mean <- rowMeans(dataC_all)
  ind <- which(overall_mean == 0)

  if (length(ind) > 0) {
    dataC_all <- dataC_all[-ind,]  # Apply same filtering to all cancer samples
  }

  # Apply logarithmic transformation
  dataC_all <- log2(dataC_all + 1)

  stopifnot(all(!is.na(dataC_all)))

  # Returns: filtered and log2-transformed cancer-only miRNA expression matrix
  return(dataC_all)
}