clear_clinical_data_columns <- function(
  clinical_data
) {
  # Remove columns with all NA values
  ind <- unlist(lapply(clinical_data, function(col) { !all(is.na(col)) }))
  clinical_data <- clinical_data[, ind]

  # Remove single value columns
  ind <- unlist(lapply(clinical_data, function(col) { length(unique(col)) > 1 }))
  clinical_data <- clinical_data[, ind]

  # Remove id-column (i.e. diagnosis_id...)
  ind <- !grepl("_id$", colnames(clinical_data))
  clinical_data <- clinical_data[, ind]
}


# Merge clinical data (legacy and updated)
merge_clinical_datas <- function(clinical_data_legacy,
                                 clinical_data_updated,
                                 common_pz_omics) {

  # Columns to remove in clinical data updated 
  columns_to_remove_updated <- c(
    'asian_race',
    'caucasian_race',
    'obesity_class_1',
    'viral_hepatitis_serology'
  )

  clinical_data_updated <- clinical_data_updated[, setdiff(colnames(clinical_data_updated), columns_to_remove_updated)]

  # Clear clinical_data_updated row names (patients are with normal/cancer info)
  rownames(clinical_data_updated) <- str_extract(rownames(clinical_data_updated), PZ_TCGA_REGEX)

  clinical_data_legacy <- clinical_data_legacy[rownames(clinical_data_legacy) %in% common_pz_omics,]
  clinical_data_updated <- clinical_data_updated[rownames(clinical_data_updated) %in% common_pz_omics,]

  missing_rows <- setdiff(rownames(clinical_data_legacy), rownames(clinical_data_updated))

  # Create a new data frame with the missing rows filled with "Not reported"
  missing_data <- data.frame(matrix("not reported",
                                    nrow = length(missing_rows),
                                    ncol = ncol(clinical_data_updated)))

  colnames(missing_data) <- colnames(clinical_data_updated)
  rownames(missing_data) <- missing_rows

  # Add the missing rows to clinical_data_updated
  clinical_data_updated <- rbind(clinical_data_updated, missing_data)

  clinical_data_legacy <- clinical_data_legacy[common_pz_omics,]
  clinical_data_updated <- clinical_data_updated[common_pz_omics,]

  final_clinical <- cbind(clinical_data_updated, clinical_data_legacy)
  final_clinical <- clear_clinical_data_columns(final_clinical)

  # Standardize all missing representations
  final_clinical[final_clinical == "not reported"] <- NA
  final_clinical[final_clinical == "Not Reported"] <- NA
  final_clinical[final_clinical == ""] <- NA

  # Run the analysis
  analyze_missing_values(
    clinical_data = final_clinical,
    threshold = 0.7,
    save_plots = TRUE,
    output_dir = paste0(dirPlots, "missing_analysis")
  )

  return(final_clinical_matrix = final_clinical)
}


# Fix data types
fix_clinical_data_types <- function(data) {

  # Fix BMI - convert to numeric
  if ("bmi" %in% colnames(data)) {
    data$bmi <- as.numeric(as.character(data$bmi))
  }

  # Fix age - convert to numeric
  if ("age_at_index" %in% colnames(data)) {
    data$age_at_index <- as.numeric(as.character(data$age_at_index))
  }

  # Ensure all other character columns that should be factors are converted
  char_cols <- sapply(data, is.character)
  if (any(char_cols)) {
    char_col_names <- names(data)[char_cols]
    cat("Converting character columns to factors:", char_col_names, "\n")

    for (col in char_col_names) {
      data[[col]] <- factor(data[[col]])
    }
  }

  return(data)
}


compute_final_clinical_matrix <- function(common_pz_omics,
                                          clinical_legacy_path = clinical_data_matrix_path,
                                          clinical_updated_path = clinical_data_integrated_matrix_path) {

  ################### Read legacy clinical data ################### 
  clinical_legacy_data <- read.table(clinical_legacy_path,
                                     header = T,
                                     sep = "\t",
                                     check.names = F,
                                     row.names = 1,
                                     quote = "")

  clinical_legacy_data <- clear_clinical_data_columns(clinical_legacy_data)

  columns_to_remove <- c(
    'synchronous_malignancy',
    'days_to_last_follow_up',
    'age_at_diagnosis',
    'primary_diagnosis',
    'year_of_diagnosis',
    'ajcc_staging_system_edition',
    'morphology',
    'ethnicity',
    'vital_status',
    'days_to_birth',
    'year_of_birth',
    'year_of_death',
    'days_to_death',
    'bcr_patient_barcode')

  # Remove selected columns
  clinical_legacy <- clinical_legacy_data[, setdiff(colnames(clinical_legacy_data), columns_to_remove)]

  ################### Read updated clinical data ###################
  clinical_updated <- read.table(clinical_updated_path,
                                 header = T,
                                 sep = "\t",
                                 check.names = F,
                                 row.names = 1,
                                 quote = "")


  ###################### Merge clinical data #######################
  final_clinical_matrix <- merge_clinical_datas(clinical_legacy,
                                                clinical_updated,
                                                common_pz_omics)

  #################### Normalize clinical data #####################


  # Convert variables to appropriate R data types
  final_clinical_matrix$obesity_class_2 <- factor(final_clinical_matrix$obesity_class_2,
                                                  exclude = NA,
                                                  levels = c("Normal", "Obese"))

  final_clinical_matrix$ajcc_pathologic_stage <- factor(final_clinical_matrix$ajcc_pathologic_stage,
                                                        exclude = NA,
                                                        levels = c("Stage I",
                                                                   "Stage II",
                                                                   "Stage III",
                                                                   "Stage IIIA",
                                                                   "Stage IIIC",
                                                                   "Stage IV"),
                                                        ordered = TRUE)

  final_clinical_matrix$ajcc_pathologic_t <- factor(final_clinical_matrix$ajcc_pathologic_t,
                                                    exclude = NA,
                                                    levels = c("T1", "T2", "T3", "T3a", "T3b", "T4"),
                                                    ordered = TRUE)

  final_clinical_matrix$ajcc_pathologic_n <- factor(final_clinical_matrix$ajcc_pathologic_n,
                                                    exclude = NA,
                                                    levels = c("N0", "N1", "NX"),
                                                    ordered = TRUE)

  final_clinical_matrix$ajcc_pathologic_m <- factor(final_clinical_matrix$ajcc_pathologic_m,
                                                    exclude = NA,
                                                    levels = c("M0", "MX", "M1"),
                                                    ordered = TRUE)

  final_clinical_matrix$gender <- factor(final_clinical_matrix$gender,
                                         exclude = NA,
                                         levels = c("female", "male"),
                                         ordered = TRUE)

  # Standardize binary columns

  # Define the columns to standardize
  columns_to_fix <- c("hepatitis_c",
                      "hepatitis_b",
                      "alcoholic_liver",
                      "prior_malignancy",
                      "treatments_radiation_treatment_or_therapy",
                      "treatments_pharmaceutical_treatment_or_therapy")

  # Apply standardization to each specified column
  for (col in columns_to_fix) {
    if (col %in% names(final_clinical_matrix)) {
      # Convert to character to handle any factor issues
      final_clinical_matrix[[col]] <- as.character(final_clinical_matrix[[col]])

      # Standardize values
      final_clinical_matrix[[col]] <- ifelse(final_clinical_matrix[[col]] %in% c("No", "no"), "No", "Yes")
    }
  }

  # Convert binary variables to factors
  binary_vars <- c("cirrhosis",
                   "hepatitis_c",
                   "hepatitis_b",
                   "family_history",
                   "alcoholic_liver",
                   "prior_malignancy",
                   "treatments_radiation_treatment_or_therapy",
                   "treatments_pharmaceutical_treatment_or_therapy")

  for (var in binary_vars) {
    final_clinical_matrix[[var]] <- factor(final_clinical_matrix[[var]],
                                           levels = c("No", "Yes"))
  }

  # Reorder clinical data by patients
  final_clinical_matrix <- final_clinical_matrix[cancer_patient_ids,]

  cat("Clinical matrix:", dim(final_clinical_matrix), "\n")

  # Remove columns not used in similarity matrix computation
  columns_to_remove <- c("race")

  final_clinical_matrix <- final_clinical_matrix[, setdiff(colnames(final_clinical_matrix), columns_to_remove)]

  # Apply the fix
  final_clinical_matrix <- fix_clinical_data_types(final_clinical_matrix)
  final_clinical_matrix <- final_clinical_matrix[cancer_patient_ids,]

  return(final_clinical_matrix)
}
