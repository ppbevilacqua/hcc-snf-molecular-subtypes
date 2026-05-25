########################### Base survival analysis ###########################

# ============================================================================
# 1. DATA PREPROCESSING FUNCTIONS
# ============================================================================

# Function to prepare and merge clinical data with clusters
prepare_clinical_data <- function(clinical_legacy, clinical_updated, patient_to_cluster) {

  # Filter clinical legacy with clustered patients
  clinical_legacy <- clinical_legacy[rownames(patient_to_cluster),]

  # Adjust clinical updated patient ids
  rownames(clinical_updated) <- str_extract(rownames(clinical_updated), PZ_TCGA_REGEX)

  # Create a new data frame with the missing rows filled with NA
  missing_rows <- setdiff(rownames(clinical_legacy), rownames(clinical_updated))

  missing_data <- data.frame(matrix(NA,
                                    nrow = length(missing_rows),
                                    ncol = ncol(clinical_updated)))

  colnames(missing_data) <- colnames(clinical_updated)
  rownames(missing_data) <- missing_rows

  # Add the missing rows to clinical_updated
  clinical_updated <- rbind(clinical_updated, missing_data)

  # Get patient IDs from rownames
  legacy_patients <- rownames(clinical_legacy)
  cluster_patients <- rownames(patient_to_cluster)

  # Find common patients
  common_patients <- intersect(legacy_patients, cluster_patients)

  # Merge clinical data for common patients
  clinical_combined <- cbind(
    clinical_legacy[common_patients, , drop = FALSE],
    clinical_updated[common_patients, , drop = FALSE]
  )

  # Remove duplicate columns if any
  clinical_combined <- clinical_combined[, !duplicated(colnames(clinical_combined))]

  # Add patient IDs as a column
  clinical_combined$patient_id <- rownames(clinical_combined)

  # Add SNF cluster assignments
  if (nrow(patient_to_cluster) == nrow(clinical_combined)) {
    cluster_df <- data.frame(
      patient_id = rownames(patient_to_cluster),
      snf_cluster = patient_to_cluster$cluster,
      stringsAsFactors = FALSE
    )
  } else {
    stop("Mismatch between number of clusters and patients. Please check cluster assignment.")
  }

  # Merge with cluster information
  clinical_combined <- merge(clinical_combined, cluster_df,
                             by = "patient_id", all.x = TRUE)

  # Convert cluster to factor with meaningful labels
  clinical_combined$snf_cluster <- factor(clinical_combined$snf_cluster,
                                          levels = sort(unique(final_clusters)),
                                          labels = paste0("C", sort(unique(final_clusters))))

  # Set rownames back
  rownames(clinical_combined) <- clinical_combined$patient_id

  return(clinical_combined)
}

# Function to check and report missing data
check_missing_data <- function(data) {

  # Replace "not reported" with NA for consistent checking
  data[data == "not reported"] <- NA

  # Calculate missing data statistics (variables)
  missing_summary <- data.frame(
    variable = colnames(data),
    n_missing = sapply(data, function(x) sum(is.na(x))),
    pct_missing = sapply(data, function(x) round(100 * sum(is.na(x)) / length(x), 2)),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(pct_missing))

  # Mark variables to exclude (>30% missing)
  missing_summary$keep <- missing_summary$pct_missing <= 30

  cat("\n=== Missing Data Summary ===\n")
  cat(sprintf("Total variables: %d\n", nrow(missing_summary)))
  cat(sprintf("Variables with >30%% missing: %d\n", sum(!missing_summary$keep)))

  if (sum(!missing_summary$keep) > 0) {
    cat("\nVariables to be excluded due to high missingness:\n")
    print(missing_summary[!missing_summary$keep, c("variable", "pct_missing")])
  }

  # ----  Patient-level survival eligibility check ----
  invalid_patients <- with(data, {
    # Case 1: vital_status missing
    is_na_vital <- is.na(vital_status)

    # Case 2: Dead but no days_to_death
    dead_no_time <- vital_status == "Dead" & is.na(days_to_death)

    # Case 3: Alive but no days_to_last_follow_up
    alive_no_time <- vital_status == "Alive" & is.na(days_to_last_follow_up)

    # Combine conditions
    is_na_vital | dead_no_time | alive_no_time
  })

  n_invalid <- sum(invalid_patients, na.rm = TRUE)

  cat("\n=== Survival Analysis Eligibility ===\n")
  cat(sprintf("Total patients: %d\n", nrow(data)))
  cat(sprintf("Patients to exclude from survival analysis: %d\n", n_invalid))

  if (n_invalid > 0) {
    cat("\nExample patient IDs to exclude:\n")
    if ("patient_id" %in% colnames(data)) {
      print(head(data$patient_id[invalid_patients], 10))
    } else {
      print(head(which(invalid_patients), 10))
    }
  }

  return(list(
    missing_summary = missing_summary,
    invalid_patients = invalid_patients
  ))
}

# Function to create survival variables
create_survival_data <- function(data) {

  # Replace "not reported" with NA for consistency
  data[data == "not reported"] <- NA

  # Required columns for survival
  required_cols <- c("vital_status", "days_to_death", "days_to_last_follow_up")

  if (!all(required_cols %in% colnames(data))) {
    stop("Missing required survival columns: ",
         paste(required_cols, collapse = ", "))
  }

  # Ensure proper data types
  data <- data %>%
    mutate(
      vital_status = as.character(vital_status),
      days_to_death = suppressWarnings(as.numeric(as.character(days_to_death))),
      days_to_last_follow_up = suppressWarnings(as.numeric(as.character(days_to_last_follow_up))),

      # Event indicator: 1 = death, 0 = alive
      os_event = case_when(
        tolower(vital_status) %in% c("dead", "deceased", "1") ~ 1,
        tolower(vital_status) %in% c("alive", "living", "0") ~ 0,
        TRUE ~ NA_real_
      ),

      # Survival time in months
      os_time = case_when(
        os_event == 1 ~ days_to_death / 30.44,
        os_event == 0 ~ days_to_last_follow_up / 30.44,
        TRUE ~ NA_real_
      ),

      # Survival time in years
      os_years = os_time / 12
    )

  # Final check: remove any rows still missing survival time/event
  n_before <- nrow(data)
  data <- data %>%
    filter(!is.na(os_time) & !is.na(os_event) & os_time > 0)
  n_after <- nrow(data)

  cat("\n=== Survival Variables Created ===\n")
  cat(sprintf("Removed %d patients with incomplete survival variables after filtering\n",
              n_before - n_after))
  cat(sprintf("Final cohort size: %d patients\n", n_after))

  return(data)
}

# Function to process clinical variables
process_clinical_variables <- function(data, missing_summary) {

  # Replace "not reported" with NA
  data[data == "not reported"] <- NA

  # Keep only variables with acceptable missing data
  vars_to_keep <- missing_summary %>%
    filter(keep) %>%
    pull(variable)

  # Always keep essential variables
  essential_vars <- c("patient_id", "snf_cluster", "os_time", "os_event", "os_years")
  vars_to_keep <- unique(c(vars_to_keep, essential_vars[essential_vars %in% names(data)]))

  data_processed <- data[, vars_to_keep, drop = FALSE]

  # Process age variable if present
  age_cols <- c("age_at_index", "age_at_initial_pathologic_diagnosis", "age")
  age_match <- age_cols[age_cols %in% colnames(data_processed)]
  if (length(age_match) > 0) {
    age_col <- age_match[1]
    data_processed[[age_col]] <- as.numeric(as.character(data_processed[[age_col]]))
    data_processed$age_group <- cut(
      data_processed[[age_col]],
      breaks = c(0, 50, 60, 70, Inf),
      labels = c("<50", "50-60", "60-70", ">70"),
      right = FALSE
    )
  }

  # Process stage if available
  stage_cols <- c("ajcc_pathologic_stage", "pathologic_stage", "tumor_stage")
  stage_match <- stage_cols[stage_cols %in% colnames(data_processed)]
  if (length(stage_match) > 0) {
    stage_col <- stage_match[1]
    data_processed$stage_group <- case_when(
      grepl("^Stage I$", data_processed[[stage_col]]) ~ "Stage I",
      grepl("^Stage II$", data_processed[[stage_col]]) ~ "Stage II",
      grepl("^Stage IV$", data_processed[[stage_col]]) ~ "Stage IV",
      grepl("^Stage III", data_processed[[stage_col]]) ~ "Stage III",
      TRUE ~ "Unknown"
    )
  }

  # Process gender if available
  if ("gender" %in% colnames(data_processed)) {
    data_processed$gender <- tolower(as.character(data_processed$gender))
  }

  # Process hepatitis and other liver conditions if available
  hep_vars <- c("hepatitis_b", "hepatitis_c", "cirrhosis", "alcoholic_liver", "viral_hepatitis_serology")
  for (var in hep_vars) {
    if (var %in% colnames(data_processed)) {
      data_processed[[var]] <- factor(data_processed[[var]])
    }
  }

  return(data_processed)
}

# ============================================================================
# 2. MAIN SURVIVAL ANALYSIS FUNCTIONS
# ============================================================================

# Kaplan-Meier analysis
perform_km_analysis <- function(data) {

  # Ensure cluster is factor
  data$snf_cluster <- as.factor(data$snf_cluster)

  # Check cluster distribution
  cat("\nCluster distribution:\n")
  print(table(data$snf_cluster))

  # Fit KM model
  km_fit <- survfit(Surv(os_time, os_event) ~ snf_cluster, data = data)

  # Log-rank test
  log_rank <- survdiff(Surv(os_time, os_event) ~ snf_cluster, data = data)
  p_value <- pchisq(log_rank$chisq, df = length(log_rank$n) - 1, lower.tail = FALSE)

  # Define color palette
  n_clusters <- length(levels(data$snf_cluster))
  if (n_clusters <= 8) {
    palette <- RColorBrewer::brewer.pal(n_clusters, "Dark2")
  } else if (n_clusters <= 20) {
    palette <- c(
      RColorBrewer::brewer.pal(8, "Dark2"),
      RColorBrewer::brewer.pal(n_clusters - 8, "Set3")
    )
  } else {
    warning("Number of clusters exceeds palette capacity; colors will be recycled.")
    palette <- rep(RColorBrewer::brewer.pal(12, "Set3"), length.out = n_clusters)
  }

  # Create KM plot
  km_plot <- ggsurvplot(
    km_fit,
    data = data,
    pval = TRUE,
    pval.method = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.3,
    palette = palette,
    legend.title = "Cluster",
    legend.labs = levels(data$snf_cluster),
    xlab = "Time (months)",
    ylab = "Overall Survival Probability",
    title = "Kaplan-Meier Survival Curves by SNF Cluster",
    break.time.by = 12,
    ggtheme = theme_bw(),
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    conf.int = TRUE,
    conf.int.alpha = 0.15,
    surv.median.line = "hv",
    tables.theme = theme_cleantable()
  )

  # Median survival summary
  km_summary <- summary(km_fit)
  median_surv <- km_summary$table[, "median"]

  # Replace NA with "NR" (Not Reached)
  median_surv_out <- ifelse(is.na(median_surv), "NR", round(median_surv, 1))

  cluster_names <- gsub("snf_cluster=", "", rownames(km_summary$table))

  surv_summary <- data.frame(
    Cluster = cluster_names,
    N = km_summary$table[, "records"],
    Events = km_summary$table[, "events"],
    Median_Survival = median_surv_out,
    CI_Lower = round(km_summary$table[, "0.95LCL"], 1),
    CI_Upper = round(km_summary$table[, "0.95UCL"], 1),
    stringsAsFactors = FALSE
  )

  return(list(
    fit = km_fit,
    plot = km_plot,
    log_rank_p = p_value,
    summary_table = surv_summary
  ))
}

# Cox regression analysis
perform_cox_analysis <- function(data) {

  # Check events per cluster
  event_table <- table(data$snf_cluster, data$os_event)
  cat("\nEvents per cluster:\n")
  print(event_table)

  if (any(rowSums(event_table[, "1", drop = FALSE]) == 0)) {
    cat("\nWarning: One or more clusters have zero events. Cox model may be unstable.\n")
  }

  # Force snf_cluster to factor
  data$snf_cluster <- factor(data$snf_cluster)

  # Univariate Cox
  cox_uni <- coxph(Surv(os_time, os_event) ~ snf_cluster, data = data)

  # Potential covariates
  potential_covars <- c("age_at_diagnosis", "age_at_initial_pathologic_diagnosis",
                        "gender", "stage_group", "hepatitis_b", "hepatitis_c",
                        "cirrhosis", "alcoholic_liver", "bmi")

  available_covars <- potential_covars[potential_covars %in% names(data)]
  cat("\nAvailable covariates for multivariate analysis:",
      paste(available_covars, collapse = ", "), "\n")

  # Convert categorical vars to factor
  categorical_covars <- c("gender", "stage_group", "hepatitis_b", "hepatitis_c",
                          "cirrhosis", "alcoholic_liver")
  for (covar in categorical_covars) {
    if (covar %in% colnames(data)) {
      data[[covar]] <- as.factor(data[[covar]])
    }
  }

  # Filter covariates by missingness
  covar_missing <- sapply(available_covars, function(x) sum(is.na(data[[x]])))
  covars_to_use <- available_covars[covar_missing < nrow(data) * 0.3]

  results <- list(univariate = cox_uni, multivariate = NULL)

  # Fit multivariate model if covariates available
  if (length(covars_to_use) > 0) {
    formula_str <- paste("Surv(os_time, os_event) ~ snf_cluster +",
                         paste(covars_to_use, collapse = " + "))

    tryCatch({
      cox_multi <- coxph(as.formula(formula_str), data = data)
      results$multivariate <- cox_multi
      cat("\nMultivariate model includes:", paste(covars_to_use, collapse = ", "), "\n")
    }, error = function(e) {
      cat("\nWarning: Could not fit multivariate model:", e$message, "\n")
    })
  } else {
    cat("\nNo covariates passed missingness filter — skipping multivariate model.\n")
  }

  return(results)
}

# Create forest plot
create_forest_plot <- function(cox_model, title = "Forest Plot") {

  # Extract model information
  coef_summary <- summary(cox_model)

  # Prepare data for forest plot
  forest_data <- data.frame(
    Variable = rownames(coef_summary$coefficients),
    HR = exp(coef_summary$coefficients[, "coef"]),
    Lower = exp(confint(cox_model)[, 1]),
    Upper = exp(confint(cox_model)[, 2]),
    P_value = coef_summary$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )

  # Clean variable names
  forest_data$Variable <- gsub("snf_cluster", "", forest_data$Variable)
  forest_data$Variable <- gsub("_", " ", forest_data$Variable)

  # Add reference line text
  forest_data$HR_text <- sprintf("%.2f (%.2f-%.2f)",
                                 forest_data$HR,
                                 forest_data$Lower,
                                 forest_data$Upper)
  forest_data$P_text <- ifelse(forest_data$P_value < 0.001, "<0.001",
                               sprintf("%.3f", forest_data$P_value))

  # Create plot
  p <- ggplot(forest_data, aes(x = HR, y = reorder(Variable, HR))) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.5) +
    scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10)) +
    labs(x = "Hazard Ratio (95% CI)", y = "", title = title) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  # Add text annotations
  max_x <- max(forest_data$Upper, na.rm = TRUE) * 1.5
  p <- p +
    geom_text(aes(x = max_x, label = HR_text), hjust = 0, size = 3) +
    geom_text(aes(x = max_x * 1.5, label = paste("p =", P_text)), hjust = 0, size = 3) +
    expand_limits(x = max_x * 2)

  return(list(plot = p, data = forest_data))
}

# Clinical characterization
compare_clinical_features <- function(data) {

  # Exclude variables we don't want in the comparison
  exclude_vars <- c("patient_id", "os_time", "os_event", "os_years",
                    "days_to_diagnosis", "days_to_death", "days_to_last_follow_up")

  vars_to_analyze <- setdiff(names(data), exclude_vars)

  vars_continuous <- c()
  vars_categorical <- c()

  for (var in vars_to_analyze) {
    if (var == "snf_cluster") next

    if (is.numeric(data[[var]])) {
      # Continuous vs categorical numeric
      if (length(unique(data[[var]][!is.na(data[[var]])])) > 10) {
        vars_continuous <- c(vars_continuous, var)
      } else {
        vars_categorical <- c(vars_categorical, var)
      }
    } else {
      vars_categorical <- c(vars_categorical, var)
    }
  }

  # Convert numeric-coded categorical vars to factor
  for (cat_var in vars_categorical) {
    if (!is.factor(data[[cat_var]])) {
      data[[cat_var]] <- as.factor(data[[cat_var]])
    }
  }

  cat("\nContinuous variables:", paste(vars_continuous, collapse = ", "), "\n")
  cat("Categorical variables:", paste(vars_categorical, collapse = ", "), "\n")

  # Ensure snf_cluster is factor
  data$snf_cluster <- as.factor(data$snf_cluster)

  # Check if we have >1 cluster level
  if (length(levels(data$snf_cluster)) < 2) {
    cat("\nWarning: Only one SNF cluster present — skipping comparison.\n")
    return(NULL)
  }

  # Create comparison table
  if (length(c(vars_continuous, vars_categorical)) > 0) {
    table1 <- CreateTableOne(
      vars = c(vars_continuous, vars_categorical),
      strata = "snf_cluster",
      data = data,
      test = TRUE
    )

    table1_print <- print(table1,
                          quote = FALSE,
                          noSpaces = TRUE,
                          printToggle = FALSE,
                          showAllLevels = TRUE)

    return(list(
      table = table1,
      summary = table1_print,
      continuous_vars = vars_continuous,
      categorical_vars = vars_categorical
    ))
  } else {
    cat("\nNo clinical variables available for comparison\n")
    return(NULL)
  }
}

# ============================================================================
# 3. MAIN ANALYSIS PIPELINE
# ============================================================================

run_hcc_survival_analysis <- function(clinical_legacy, clinical_updated, patient_to_cluster) {

  cat("SURVIVAL ANALYSIS FOR SNF-DERIVED HCC SUBTYPES\n")

  # Initialize results list
  results <- list()

  # Step 1: Check clinical data
  cat("\n1. CHECK CLINICAL DATA\n")

  cat(sprintf("Legacy data: %d patients, %d variables\n",
              nrow(clinical_legacy), ncol(clinical_legacy)))
  cat(sprintf("Updated data: %d patients, %d variables\n",
              nrow(clinical_updated), ncol(clinical_updated)))

  # Step 2: Prepare and merge data
  cat("\n2. MERGING DATA AND ADDING CLUSTERS\n")

  clinical_data <- prepare_clinical_data(clinical_legacy, clinical_updated, patient_to_cluster)

  cat(sprintf("Combined data: %d patients, %d variables\n",
              nrow(clinical_data), ncol(clinical_data)))

  cat(sprintf("Number of clusters: %d\n", length(unique(clinical_data$snf_cluster))))

  # Step 3: Check missing data
  cat("\n3. MISSING DATA ASSESSMENT\n")

  missing_summary_results <- check_missing_data(clinical_data)

  clinical_data <- clinical_data[!missing_summary_results$invalid_patients,]

  missing_summary <- missing_summary_results$missing_summary
  results$missing_summary <- missing_summary

  # Step 4: Create survival data
  cat("\n4. CREATING SURVIVAL VARIABLES\n")

  survival_data <- create_survival_data(clinical_data)
  survival_data <- process_clinical_variables(survival_data, missing_summary)

  results$survival_data <- survival_data

  # Print survival summary
  cat("\nSurvival data summary:\n")
  cat(sprintf("  Total events: %d (%.1f%%)\n",
              sum(survival_data$os_event),
              100 * mean(survival_data$os_event)))
  cat(sprintf("  Median follow-up: %.1f months\n",
              median(survival_data$os_time[survival_data$os_event == 0])))
  cat(sprintf("  Median survival: %.1f months\n",
              median(survival_data$os_time[survival_data$os_event == 1])))

  # Step 5: Kaplan-Meier analysis
  cat("\n5. KAPLAN-MEIER ANALYSIS\n")

  km_results <- perform_km_analysis(survival_data)
  results$km_results <- km_results

  cat(sprintf("\nLog-rank test p-value: %.4f\n", km_results$log_rank_p))

  if (km_results$log_rank_p < 0.01) {
    cat("=> Statistically significant survival differences between clusters\n")
  } else {
    cat("=> No statistically significant survival differences\n")
  }

  cat("\nMedian survival by cluster:\n")
  print(km_results$summary_table)

  # Step 6: Cox regression
  cat("\n6. COX REGRESSION ANALYSIS\n")

  cox_results <- perform_cox_analysis(survival_data)
  results$cox_results <- cox_results

  cat("\nUnivariate Cox model:\n")
  print(summary(cox_results$univariate))

  if (!is.null(cox_results$multivariate)) {
    cat("\nMultivariate Cox model:\n")
    print(summary(cox_results$multivariate))
  }

  # Step 7: Forest plots
  cat("\n7. GENERATING FOREST PLOTS\n")

  forest_uni <- create_forest_plot(cox_results$univariate,
                                   "Univariate Cox Regression - SNF Clusters")
  results$forest_plots <- list(univariate = forest_uni)

  if (!is.null(cox_results$multivariate)) {
    forest_multi <- create_forest_plot(cox_results$multivariate,
                                       "Multivariate Cox Regression")
    results$forest_plots$multivariate <- forest_multi
  }

  # Step 8: Clinical characterization
  cat("\n8. CLINICAL CHARACTERIZATION\n")

  clinical_comparison <- compare_clinical_features(survival_data)
  results$clinical_comparison <- clinical_comparison

  if (!is.null(clinical_comparison)) {
    cat("\nClinical features comparison completed\n")
  }

  # Step 9: Check proportional hazards assumption
  cat("\n9. MODEL DIAGNOSTICS\n")

  ph_test <- cox.zph(cox_results$univariate)
  results$ph_test <- ph_test

  cat("\nProportional hazards test:\n")
  print(ph_test)

  cat("ANALYSIS COMPLETE!\n")

  return(results)
}

# ============================================================================
# 4. VISUALIZATION AND REPORTING FUNCTIONS
# ============================================================================

# Function to save all results
save_results <- function(results, output_dir = "HCC_survival_results") {

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save KM plot - Handle ggsurvplot object properly
  tryCatch({
    pdf(file.path(output_dir, "kaplan_meier_curves.pdf"), width = 10, height = 8)

    if ("ggsurvplot" %in% class(results$km_results$plot)) {
      # For ggsurvplot objects, print the entire object
      print(results$km_results$plot, newpage = FALSE)
    } else {
      print(results$km_results$plot)
    }
    dev.off()
  }, error = function(e) {
    dev.off()
    cat("Note: Could not save KM plot in standard format, trying alternative...\n")

    # Alternative method: save plot and table separately
    if ("ggsurvplot" %in% class(results$km_results$plot)) {
      ggsave(filename = file.path(output_dir, "kaplan_meier_plot.pdf"),
             plot = results$km_results$plot$plot,
             width = 10, height = 6)

      if (!is.null(results$km_results$plot$table)) {
        ggsave(filename = file.path(output_dir, "risk_table.pdf"),
               plot = results$km_results$plot$table,
               width = 10, height = 3)
      }
    }
  })

  # Save forest plots
  if (!is.null(results$forest_plots$univariate)) {
    tryCatch({
      # Get the data
      forest_data <- results$forest_plots$univariate$data

      # Create a clean forest plot
      uni_plot_clean <- ggplot(forest_data, aes(x = HR, y = reorder(Variable, HR))) +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.5) +
        scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8)) +
        labs(x = "Hazard Ratio (95% CI)",
             y = "",
             title = "Univariate Cox Regression - SNF Clusters") +
        theme_bw() +
        theme(
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")
        )

      pdf(file.path(output_dir, "forest_plot_univariate.pdf"), width = 8, height = 6)
      print(uni_plot_clean)
      dev.off()

      # Save the data table
      forest_table <- forest_data %>%
        mutate(
          HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, Lower, Upper),
          P_value = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value))
        ) %>%
        select(Variable, HR_CI, P_value)

      write.csv(forest_table,
                file.path(output_dir, "forest_plot_univariate_data.csv"),
                row.names = FALSE)

    }, error = function(e) {
      cat("Warning: Could not create univariate forest plot:", e$message, "\n")
    })
  }

  if (!is.null(results$forest_plots$multivariate)) {
    tryCatch({
      # Get the multivariate forest plot data
      forest_data_multi <- results$forest_plots$multivariate$data

      # Calculate appropriate x-axis limits
      x_min <- min(0.1, min(forest_data_multi$Lower, na.rm = TRUE) * 0.9)
      x_max <- max(10, max(forest_data_multi$Upper, na.rm = TRUE) * 1.1)

      # Create a clean forest plot
      multi_plot_clean <- ggplot(forest_data_multi, aes(x = HR, y = reorder(Variable, HR))) +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.5) +
        scale_x_log10(
          breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16),
          limits = c(x_min, x_max)
        ) +
        labs(x = "Hazard Ratio (95% CI)",
             y = "",
             title = "Multivariate Cox Regression") +
        theme_bw() +
        theme(
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(hjust = 1),
          plot.margin = unit(c(1, 1, 1, 1), "cm")
        )

      # Save the clean plot
      pdf(file.path(output_dir, "forest_plot_multivariate.pdf"), width = 10, height = 8)
      print(multi_plot_clean)
      dev.off()

      # Create a data table with HR and p-values
      table_data <- forest_data_multi %>%
        mutate(
          `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR, Lower, Upper),
          `P-value` = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value))
        ) %>%
        select(Variable, `HR (95% CI)`, `P-value`)

      # Save as CSV
      write.csv(table_data,
                file.path(output_dir, "forest_plot_multivariate_data.csv"),
                row.names = FALSE)

      # Try to create combined plot with table
      if (require(gridExtra, quietly = TRUE)) {
        tryCatch({
          # Create table grob
          table_grob <- tableGrob(
            table_data,
            rows = NULL,
            theme = ttheme_minimal(
              base_size = 9,
              padding = unit(c(3, 3), "mm")
            )
          )

          # Create combined plot
          pdf(file.path(output_dir, "forest_plot_multivariate_with_table.pdf"),
              width = 14, height = max(8, nrow(forest_data_multi) * 0.5 + 2))

          grid.arrange(
            multi_plot_clean,
            table_grob,
            ncol = 2,
            widths = c(3, 2),
            top = textGrob("Multivariate Cox Regression Analysis",
                           gp = gpar(fontsize = 14, fontface = "bold"))
          )
          dev.off()
        }, error = function(e) {
          cat("Note: Could not create combined forest plot with table\n")
        })
      }

    }, error = function(e) {
      cat("Warning: Could not create multivariate forest plot:", e$message, "\n")
    })
  }

  # Save survival summary table
  tryCatch({
    write.csv(results$km_results$summary_table,
              file.path(output_dir, "survival_summary.csv"),
              row.names = FALSE)
  }, error = function(e) {
    cat("Warning: Could not save survival summary table\n")
  })

  # Save clinical comparison if available
  if (!is.null(results$clinical_comparison)) {
    tryCatch({
      write.csv(results$clinical_comparison$summary,
                file.path(output_dir, "clinical_comparison.csv"))
    }, error = function(e) {
      cat("Warning: Could not save clinical comparison table\n")
    })
  }

  # Missing data bar plot
  if (!is.null(results$missing_summary)) {
    miss_df <- results$missing_summary

    if (!("pct_missing" %in% colnames(miss_df))) {
      stop("missing_summary must contain column 'pct_missing'")
    }

    miss_df$variable <- factor(miss_df$variable, levels = rev(miss_df$variable[order(miss_df$pct_missing)]))

    base_size <- 12

    miss_plot <- ggplot(miss_df, aes(x = variable, y = pct_missing)) +
      geom_col(fill = "#d73027") +
      coord_flip() +
      labs(x = "", y = "% missing", title = "Missing Data per Variable") +
      theme_bw(base_size = base_size) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = base_size - 1))

    miss_pdf <- file.path(output_dir, "missing_data.pdf")
    ggplot2::ggsave(miss_pdf, plot = miss_plot, width = 7, height = 5)
  }

  # Save the complete results object
  saveRDS(results, file.path(output_dir, "complete_results.rds"))

  # Create a summary text file
  sink(file.path(output_dir, "analysis_summary.txt"))

  cat("HCC SNF SURVIVAL ANALYSIS SUMMARY\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")

  cat("Dataset Summary:\n")
  cat(sprintf("  Total patients: %d\n", nrow(results$survival_data)))
  cat(sprintf("  Number of clusters: %d\n", length(unique(results$survival_data$snf_cluster))))
  cat(sprintf("  Total events: %d (%.1f%%)\n",
              sum(results$survival_data$os_event),
              100 * mean(results$survival_data$os_event)))

  cat("\nSurvival Analysis Results:\n")
  cat(sprintf("  Log-rank p-value: %.4f\n", results$km_results$log_rank_p))

  if (results$km_results$log_rank_p < 0.01) {
    cat("  Result: Statistically significant differences between clusters\n")
  } else {
    cat("  Result: No statistically significant differences between clusters\n")
  }

  cat("\nMedian Survival by Cluster:\n")
  print(results$km_results$summary_table)

  cat("\nFiles saved:\n")
  files_saved <- list.files(output_dir, pattern = "\\.(pdf|csv|txt|rds)$")
  for (file in files_saved) {
    cat(sprintf("  - %s\n", file))
  }

  sink()

  cat(sprintf("\nResults successfully saved to: %s\n", output_dir))

  # Return the output directory path
  invisible(output_dir)
}

summarize_features <- function(data_matrix, output_file = "feature_summary.csv") {
  # Ensure the input is a data frame
  df <- as.data.frame(data_matrix, stringsAsFactors = FALSE)

  summary_list <- list()

  for (feature in names(df)) {
    col_data <- df[[feature]]

    # Detect if numeric
    suppressWarnings({
      as_num <- as.numeric(col_data)
    })

    if (is.numeric(col_data) || all(suppressWarnings(!is.na(as.numeric(na.omit(col_data)))))) {
      col_data_num <- as_num
      n <- length(col_data_num)
      missing_pct <- sum(is.na(col_data_num)) / n * 100

      if (all(is.na(col_data_num))) {
        min_val <- NA
        max_val <- NA
        median_val <- NA
      } else {
        min_val <- min(col_data_num, na.rm = TRUE)
        max_val <- max(col_data_num, na.rm = TRUE)
        median_val <- median(col_data_num, na.rm = TRUE)
      }

      summary_list[[feature]] <- data.frame(
        Feature = feature,
        Type = "Numeric",
        MissingPct = round(missing_pct, 2),
        Min = min_val,
        Max = max_val,
        Median = median_val,
        Details = NA
      )

    } else {
      # Treat as categorical
      col_data_fac <- as.factor(col_data)
      value_counts <- as.data.frame(table(col_data_fac, useNA = "ifany"))
      counts_str <- paste0(value_counts$col_data_fac, ":", value_counts$Freq, collapse = "; ")

      summary_list[[feature]] <- data.frame(
        Feature = feature,
        Type = "Categorical",
        MissingPct = NA,
        Min = NA,
        Max = NA,
        Median = NA,
        Details = counts_str
      )
    }
  }

  summary_df <- do.call(rbind, summary_list)

  write.csv(summary_df, output_file, row.names = FALSE)
}

######################### Complete survival analysis #########################

# ============================================================================
# 1: STRATIFIED SURVIVAL ANALYSIS
# ============================================================================

# Function to perform stratified analysis by clinical factors
stratified_survival_analysis <- function(survival_data, stratify_by = "stage_group") {

  if (!stratify_by %in% colnames(survival_data)) {
    cat(sprintf("Variable '%s' not found in data\n", stratify_by))
    return(NULL)
  }

  # Remove missing values for stratification variable
  data_subset <- survival_data %>%
    filter(!is.na(!!sym(stratify_by)))

  # Get unique strata
  strata <- unique(data_subset[[stratify_by]])

  # Results list
  stratified_results <- list()

  cat(sprintf("\nStratified Analysis by %s\n", stratify_by))
  cat(paste(rep("-", 40), collapse = ""), "\n")

  for (stratum in strata) {
    # Subset data
    stratum_data <- data_subset %>%
      filter(!!sym(stratify_by) == stratum)

    if (nrow(stratum_data) < 10) {
      cat(sprintf("  %s: Skipping - too few patients (n=%d)\n", stratum, nrow(stratum_data)))
      next
    }

    # Perform KM analysis
    km_fit <- survfit(Surv(os_time, os_event) ~ snf_cluster, data = stratum_data)
    log_rank <- survdiff(Surv(os_time, os_event) ~ snf_cluster, data = stratum_data)
    p_value <- 1 - pchisq(log_rank$chisq, df = length(log_rank$n) - 1)

    cat(sprintf("  %s: N = %d, Events = %d, p-value = %.4f\n",
                stratum, nrow(stratum_data), sum(stratum_data$os_event), p_value))

    stratified_results[[stratum]] <- list(
      data = stratum_data,
      km_fit = km_fit,
      p_value = p_value,
      n = nrow(stratum_data),
      events = sum(stratum_data$os_event)
    )
  }

  return(stratified_results)
}

# ============================================================================
# 2: IDENTIFY BEST AND WORST PERFORMING CLUSTERS
# ============================================================================

identify_extreme_clusters <- function(km_results, survival_data) {

  # Get median survival times
  median_surv <- km_results$summary_table

  # Identify best and worst clusters
  best_idx <- which.max(median_surv$Median_Survival)
  worst_idx <- which.min(median_surv$Median_Survival)

  best_cluster <- median_surv$Cluster[best_idx]
  worst_cluster <- median_surv$Cluster[worst_idx]

  best_cluster <- gsub("snf_cluster=", "", best_cluster)
  worst_cluster <- gsub("snf_cluster=", "", worst_cluster)

  cat("\n=== Cluster Performance Summary ===\n")
  cat(sprintf("Best survival: %s (median = %s months)\n",
              best_cluster, median_surv$Median_Survival[best_idx]))
  cat(sprintf("Worst survival: %s (median = %s months)\n",
              worst_cluster, median_surv$Median_Survival[worst_idx]))

  # Compare characteristics between best and worst
  extreme_data <- survival_data %>%
    filter(snf_cluster %in% c(best_cluster, worst_cluster)) %>%
    mutate(performance = ifelse(snf_cluster == best_cluster, "Best", "Worst"))

  # Statistical comparison
  comparison_vars <- c("age_at_diagnosis", "age_at_initial_pathologic_diagnosis",
                       "gender", "stage_group", "hepatitis_b", "hepatitis_c",
                       "cirrhosis", "alcoholic_liver")

  comparison_results <- list()

  for (var in comparison_vars) {
    if (var %in% colnames(extreme_data)) {
      if (is.numeric(extreme_data[[var]])) {
        # Wilcoxon test for continuous variables
        tryCatch({
          test <- wilcox.test(as.formula(paste(var, "~ performance")),
                              data = extreme_data)
          comparison_results[[var]] <- list(
            type = "continuous",
            p_value = test$p.value,
            best_median = median(extreme_data[[var]][extreme_data$performance == "Best"],
                                 na.rm = TRUE),
            worst_median = median(extreme_data[[var]][extreme_data$performance == "Worst"],
                                  na.rm = TRUE)
          )
        }, error = function(e) {
          cat(sprintf("  Could not compare %s\n", var))
        })
      } else {
        # Chi-square test for categorical variables
        tryCatch({
          tab <- table(extreme_data[[var]], extreme_data$performance)
          if (min(dim(tab)) > 1) {
            test <- chisq.test(tab)
            comparison_results[[var]] <- list(
              type = "categorical",
              p_value = test$p.value,
              table = tab
            )
          }
        }, error = function(e) {
          cat(sprintf("  Could not compare %s\n", var))
        })
      }
    }
  }

  return(list(
    best_cluster = best_cluster,
    worst_cluster = worst_cluster,
    comparison = comparison_results,
    data = extreme_data
  ))
}

# ============================================================================
# 3: CREATE CLUSTER PROFILES
# ============================================================================

create_cluster_profiles <- function(survival_data) {

  # Calculate key metrics for each cluster
  cluster_profiles <- survival_data %>%
    group_by(snf_cluster) %>%
    summarise(
      n_patients = n(),
      n_events = sum(os_event),
      event_rate = mean(os_event) * 100,
      median_survival = median(os_time[os_event == 1], na.rm = TRUE),
      median_followup = median(os_time[os_event == 0], na.rm = TRUE),
      .groups = "drop"
    )

  # Add age statistics if available
  age_cols <- c("age_at_diagnosis", "age_at_initial_pathologic_diagnosis")
  age_col <- age_cols[age_cols %in% colnames(survival_data)][1]

  if (!is.na(age_col) && length(age_col) > 0) {
    age_stats <- survival_data %>%
      group_by(snf_cluster) %>%
      summarise(
        mean_age = mean(!!sym(age_col), na.rm = TRUE),
        sd_age = sd(!!sym(age_col), na.rm = TRUE),
        .groups = "drop"
      )
    cluster_profiles <- left_join(cluster_profiles, age_stats, by = "snf_cluster")
  }

  # Add gender distribution if available
  if ("gender" %in% colnames(survival_data)) {
    gender_dist <- survival_data %>%
      group_by(snf_cluster) %>%
      summarise(
        male_pct = sum(tolower(gender) == "male", na.rm = TRUE) / n() * 100,
        .groups = "drop"
      )
    cluster_profiles <- left_join(cluster_profiles, gender_dist, by = "snf_cluster")
  }

  # Add stage distribution if available
  if ("stage_group" %in% colnames(survival_data)) {
    stage_dist <- survival_data %>%
      filter(stage_group %in% c("Stage III", "Stage IV")) %>%
      group_by(snf_cluster) %>%
      summarise(
        n_advanced = n(),
        .groups = "drop"
      )

    stage_pct <- survival_data %>%
      group_by(snf_cluster) %>%
      summarise(n_total = n(), .groups = "drop") %>%
      left_join(stage_dist, by = "snf_cluster") %>%
      mutate(
        advanced_stage_pct = ifelse(is.na(n_advanced), 0, n_advanced / n_total * 100)
      ) %>%
      select(snf_cluster, advanced_stage_pct)

    cluster_profiles <- left_join(cluster_profiles, stage_pct, by = "snf_cluster")
  }

  # Round numeric columns
  cluster_profiles <- cluster_profiles %>%
    mutate(across(where(is.numeric), ~round(., 1)))

  return(cluster_profiles)
}

# ============================================================================
# 4: CROSS-VALIDATION FOR MODEL STABILITY
# ============================================================================

assess_cluster_stability <- function(survival_data, n_splits = 10) {

  # Create folds
  set.seed(SEED)
  folds <- createFolds(survival_data$snf_cluster, k = n_splits)

  concordance_values <- numeric(n_splits)

  cat("\nCross-Validation (", n_splits, " folds)\n", sep = "")
  cat(paste(rep("-", 40), collapse = ""), "\n")

  for (i in 1:n_splits) {
    # Split data
    train_idx <- unlist(folds[-i])
    test_idx <- folds[[i]]

    train_data <- survival_data[train_idx,]
    test_data <- survival_data[test_idx,]

    # Fit Cox model on training data
    cox_train <- coxph(Surv(os_time, os_event) ~ snf_cluster, data = train_data)

    # Predict on test data
    pred_risk <- predict(cox_train, newdata = test_data, type = "risk")

    # Calculate concordance
    conc <- concordance(Surv(test_data$os_time, test_data$os_event) ~ pred_risk)
    concordance_values[i] <- conc$concordance

    cat(sprintf("  Fold %d: C-index = %.3f\n", i, concordance_values[i]))
  }

  cat("\n=== Cross-Validation Results ===\n")
  cat(sprintf("Mean C-index: %.3f\n", mean(concordance_values)))
  cat(sprintf("SD C-index: %.3f\n", sd(concordance_values)))
  cat(sprintf("95%% CI: [%.3f, %.3f]\n",
              mean(concordance_values) - 1.96 * sd(concordance_values),
              mean(concordance_values) + 1.96 * sd(concordance_values)))

  return(concordance_values)
}

# ============================================================================
# 5: PAIRWISE CLUSTER COMPARISONS
# ============================================================================

pairwise_cluster_comparison <- function(survival_data) {

  clusters <- sort(unique(survival_data$snf_cluster))
  n_clusters <- length(clusters)

  # Initialize matrix for p-values
  p_matrix <- matrix(NA, n_clusters, n_clusters)
  rownames(p_matrix) <- clusters
  colnames(p_matrix) <- clusters

  cat("\nPairwise Log-Rank Tests\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")

  for (i in 1:(n_clusters - 1)) {
    for (j in (i + 1):n_clusters) {
      # Subset data for two clusters
      pair_data <- survival_data %>%
        filter(snf_cluster %in% c(clusters[i], clusters[j]))

      # Perform log-rank test
      log_rank <- survdiff(Surv(os_time, os_event) ~ snf_cluster, data = pair_data)
      p_value <- 1 - pchisq(log_rank$chisq, df = 1)

      p_matrix[i, j] <- p_value
      p_matrix[j, i] <- p_value

      if (p_value < 0.01) {
        cat(sprintf("%s vs %s: p = %.4f *\n", clusters[i], clusters[j], p_value))
      }
    }
  }

  # Add diagonal
  diag(p_matrix) <- 1

  return(p_matrix)
}

# ============================================================================
# 6: SURVIVAL RATES AT SPECIFIC TIME POINTS
# ============================================================================

calculate_survival_rates <- function(km_fit, time_points = c(12, 36, 60)) {

  # Get survival estimates at specific time points
  surv_summary <- summary(km_fit, times = time_points)

  # Extract cluster names
  cluster_names <- names(km_fit$strata)
  cluster_names <- gsub("snf_cluster=", "", cluster_names)

  # Create results data frame
  results <- data.frame()

  for (i in seq_along(time_points)) {
    time_data <- data.frame(
      Cluster = cluster_names,
      Time_Months = time_points[i],
      Survival_Rate = surv_summary$surv[surv_summary$time == time_points[i]],
      Lower_CI = surv_summary$lower[surv_summary$time == time_points[i]],
      Upper_CI = surv_summary$upper[surv_summary$time == time_points[i]]
    )
    results <- rbind(results, time_data)
  }

  # Round values
  results <- results %>%
    mutate(across(c(Survival_Rate, Lower_CI, Upper_CI), ~round(. * 100, 1)))

  return(results)
}

# ============================================================================
# 7: MAIN WRAPPER FOR COMPLETE SURVIVAL ANALYSIS
# ============================================================================

run_complete_hcc_suvival_analysis <- function(clinical_legacy,
                                              clinical_updated,
                                              patient_to_cluster,
                                              output_dir = "HCC_complete_survival_results") {

  cat("\n")
  cat("COMPREHENSIVE HCC SURVIVAL ANALYSIS WITH SNF CLUSTERS\n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # STEP 1: Run main survival analysis
  cat("\n>>> STEP 1: MAIN SURVIVAL ANALYSIS\n")
  results <- run_hcc_survival_analysis(clinical_legacy, clinical_updated, patient_to_cluster)

  # STEP 2: Identify extreme clusters
  cat("\n>>> STEP 2: IDENTIFYING EXTREME CLUSTERS\n")
  extreme_clusters <- identify_extreme_clusters(results$km_results, results$survival_data)
  results$extreme_clusters <- extreme_clusters

  # STEP 3: Create cluster profiles
  cat("\n>>> STEP 3: CREATING CLUSTER PROFILES\n")
  cluster_profiles <- create_cluster_profiles(results$survival_data)
  results$cluster_profiles <- cluster_profiles

  cat("\nCluster Profiles:\n")
  print(cluster_profiles)

  # STEP 4: Pairwise comparisons
  cat("\n>>> STEP 4: PAIRWISE CLUSTER COMPARISONS\n")
  pairwise_p <- pairwise_cluster_comparison(results$survival_data)
  results$pairwise_comparisons <- pairwise_p

  # STEP 5: Survival rates at specific time points
  cat("\n>>> STEP 5: CALCULATING SURVIVAL RATES\n")
  survival_rates <- calculate_survival_rates(results$km_results$fit, c(12, 36, 60))
  results$survival_rates <- survival_rates

  cat("\nSurvival Rates (%):\n")
  print(survival_rates)

  # STEP 6: Cross-validation
  if (nrow(results$survival_data) >= 50) {  # Only if enough samples
    cat("\n>>> STEP 6: CROSS-VALIDATION\n")
    cv_results <- assess_cluster_stability(results$survival_data, n_splits = 10)
    results$cross_validation <- cv_results
  }

  # STEP 7: Save all results
  cat("\n>>> STEP 7: SAVING RESULTS\n")

  # Save main results
  save_results(results, output_dir)

  # Save additional analysis results
  write.csv(cluster_profiles,
            file.path(output_dir, "cluster_profiles.csv"),
            row.names = FALSE)

  write.csv(survival_rates,
            file.path(output_dir, "survival_rates.csv"),
            row.names = FALSE)

  write.csv(pairwise_p,
            file.path(output_dir, "pairwise_comparisons.csv"))

  # Create final summary
  create_final_summary(results, output_dir)

  cat("\n")
  cat(sprintf("ANALYSIS COMPLETE! All results saved to: %s\n", output_dir))

  return(results)
}

# ============================================================================
# 8: CREATE FINAL SUMMARY REPORT
# ============================================================================

create_final_summary <- function(results, output_dir) {

  summary_file <- file.path(output_dir, "complete_analysis_summary.txt")

  sink(summary_file)

  cat("COMPREHENSIVE HCC SURVIVAL ANALYSIS - FINAL REPORT\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  cat("1. DATASET OVERVIEW\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  cat(sprintf("  Total patients: %d\n", nrow(results$survival_data)))
  cat(sprintf("  Number of clusters: %d\n", length(unique(results$survival_data$snf_cluster))))
  cat(sprintf("  Total events: %d (%.1f%%)\n",
              sum(results$survival_data$os_event),
              100 * mean(results$survival_data$os_event)))

  cat("\n2. KEY FINDINGS\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  cat(sprintf("  Log-rank p-value: %.4f\n", results$km_results$log_rank_p))

  if (!is.null(results$extreme_clusters)) {
    cat(sprintf("  Best cluster: %s\n", results$extreme_clusters$best_cluster))
    cat(sprintf("  Worst cluster: %s\n", results$extreme_clusters$worst_cluster))
  }

  if (!is.null(results$cross_validation)) {
    cat(sprintf("  Cross-validation C-index: %.3f (SD: %.3f)\n",
                mean(results$cross_validation),
                sd(results$cross_validation)))
  }

  cat("\n3. CLUSTER PROFILES\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  print(results$cluster_profiles)

  cat("\n4. CLINICAL IMPLICATIONS\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")

  sink()

  cat(sprintf("\nFinal summary saved to: %s\n", summary_file))
}
