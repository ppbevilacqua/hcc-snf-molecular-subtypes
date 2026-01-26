# Clinical Feature Association Analysis for SNF clusters

#' Perform association analysis between clusters and clinical features
#'
#' @param clinical_data Data frame with clinical features and snf_cluster column
#' @param cluster_col Name of the cluster column (default: "snf_cluster")
#' @param exclude_cols Columns to exclude from analysis (default: patient ID columns)
#' @param categorical_cols Vector of categorical column names (auto-detected if NULL)
#' @param continuous_cols Vector of continuous column names (auto-detected if NULL)
#' @param output_dir Directory to save results (default: NULL)
#' @param verbose Print progress messages (default: TRUE)
#' @return List with association test results
perform_clinical_association <- function(clinical_data,
                                         cluster_col = "snf_cluster",
                                         exclude_cols = c("patient_id", "barcode", "sample_id"),
                                         categorical_cols = NULL,
                                         continuous_cols = NULL,
                                         output_dir = NULL,
                                         verbose = TRUE) {

  # Set seed for reproducibility
  set.seed(123)

  # Check if cluster column exists
  if (!cluster_col %in% names(clinical_data)) {
    stop(paste("Cluster column", cluster_col, "not found in clinical data"))
  }

  # Convert cluster to factor
  clinical_data[[cluster_col]] <- as.factor(clinical_data[[cluster_col]])

  # Get feature columns (excluding cluster and specified columns)
  feature_cols <- setdiff(names(clinical_data), c(cluster_col, exclude_cols))

  # Auto-detect variable types if not specified
  if (is.null(categorical_cols) || is.null(continuous_cols)) {
    var_types <- detect_variable_types(clinical_data[, feature_cols])

    if (is.null(categorical_cols)) {
      categorical_cols <- var_types$categorical
    }
    if (is.null(continuous_cols)) {
      continuous_cols <- var_types$continuous
    }
  }

  # Ensure columns exist
  categorical_cols <- intersect(categorical_cols, feature_cols)
  continuous_cols <- intersect(continuous_cols, feature_cols)

  if (verbose) {
    cat("Analyzing clinical associations with", length(unique(clinical_data[[cluster_col]])), "clusters\n")
    cat("Categorical variables:", length(categorical_cols), "\n")
    cat("Continuous variables:", length(continuous_cols), "\n\n")
  }

  # Initialize results
  association_results <- list(
    categorical = data.frame(),
    continuous = data.frame(),
    post_hoc = list()
  )

  # Test categorical variables
  if (length(categorical_cols) > 0) {
    if (verbose) cat("Testing categorical variables:\n")

    cat_results <- list()

    for (var in categorical_cols) {
      if (verbose) cat("  -", var, "...")

      # Remove NA values
      data_subset <- clinical_data[!is.na(clinical_data[[var]]), c(cluster_col, var)]

      if (nrow(data_subset) < 10) {
        if (verbose) cat(" skipped (insufficient data)\n")
        next
      }

      # Create contingency table
      cont_table <- table(data_subset[[cluster_col]], data_subset[[var]])

      # Perform appropriate test
      test_result <- perform_categorical_test(cont_table)

      cat_results[[var]] <- data.frame(
        variable = var,
        test = test_result$test,
        statistic = test_result$statistic,
        df = test_result$df,
        p_value = test_result$p_value,
        n_samples = nrow(data_subset),
        n_categories = ncol(cont_table),
        stringsAsFactors = FALSE
      )

      if (verbose) cat(" p =", format(test_result$p_value, digits = 3), "\n")
    }

    if (length(cat_results) > 0) {
      association_results$categorical <- bind_rows(cat_results)

      # Add multiple testing correction
      association_results$categorical$p_adjusted <- p.adjust(
        association_results$categorical$p_value,
        method = "BH"
      )
    }
  }

  # Test continuous variables
  if (length(continuous_cols) > 0) {
    if (verbose) cat("\nTesting continuous variables:\n")

    cont_results <- list()

    for (var in continuous_cols) {
      if (verbose) cat("  -", var, "...")

      # Remove NA values
      data_subset <- clinical_data[!is.na(clinical_data[[var]]), c(cluster_col, var)]

      if (nrow(data_subset) < 10) {
        if (verbose) cat(" skipped (insufficient data)\n")
        next
      }

      # Convert to numeric if needed
      data_subset[[var]] <- as.numeric(data_subset[[var]])

      # Perform Kruskal-Wallis test
      kw_result <- kruskal.test(
        formula = as.formula(paste(var, "~", cluster_col)),
        data = data_subset
      )

      # Calculate effect size (eta-squared)
      eta_squared <- calculate_eta_squared_kw(data_subset[[var]], data_subset[[cluster_col]])

      # Summary statistics per cluster
      summary_stats <- data_subset %>%
        group_by(!!sym(cluster_col)) %>%
        summarise(
          mean = mean(!!sym(var), na.rm = TRUE),
          median = median(!!sym(var), na.rm = TRUE),
          sd = sd(!!sym(var), na.rm = TRUE),
          min = min(!!sym(var), na.rm = TRUE),
          max = max(!!sym(var), na.rm = TRUE),
          n = n(),
          .groups = "drop"
        )

      cont_results[[var]] <- data.frame(
        variable = var,
        test = "Kruskal-Wallis",
        statistic = kw_result$statistic,
        df = kw_result$parameter,
        p_value = kw_result$p.value,
        eta_squared = eta_squared,
        n_samples = nrow(data_subset),
        global_mean = mean(data_subset[[var]], na.rm = TRUE),
        global_median = median(data_subset[[var]], na.rm = TRUE),
        global_sd = sd(data_subset[[var]], na.rm = TRUE),
        stringsAsFactors = FALSE
      )

      # Store summary statistics
      association_results$post_hoc[[paste0(var, "_summary")]] <- summary_stats

      # If significant, perform post-hoc pairwise comparisons
      if (kw_result$p.value < 0.01) {
        pw_result <- pairwise.wilcox.test(
          data_subset[[var]],
          data_subset[[cluster_col]],
          p.adjust.method = "BH"
        )
        association_results$post_hoc[[paste0(var, "_pairwise")]] <- pw_result
      }

      if (verbose) cat(" p =", format(kw_result$p.value, digits = 3), "\n")
    }

    if (length(cont_results) > 0) {
      association_results$continuous <- bind_rows(cont_results)

      # Add multiple testing correction
      association_results$continuous$p_adjusted <- p.adjust(
        association_results$continuous$p_value,
        method = "BH"
      )
    }
  }

  # Create summary table
  all_results <- bind_rows(
    association_results$categorical %>%
      select(variable, test, p_value, p_adjusted) %>%
      mutate(type = "categorical"),
    association_results$continuous %>%
      select(variable, test, p_value, p_adjusted) %>%
      mutate(type = "continuous")
  ) %>%
    arrange(p_adjusted)

  association_results$summary <- all_results

  # Save results if output directory provided
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    # Save main results
    write.csv(association_results$summary,
              file = file.path(output_dir, "clinical_associations_summary.csv"),
              row.names = FALSE)

    if (nrow(association_results$categorical) > 0) {
      write.csv(association_results$categorical,
                file = file.path(output_dir, "clinical_associations_categorical.csv"),
                row.names = FALSE)
    }

    if (nrow(association_results$continuous) > 0) {
      write.csv(association_results$continuous,
                file = file.path(output_dir, "clinical_associations_continuous.csv"),
                row.names = FALSE)
    }

    # Save post-hoc results
    for (name in names(association_results$post_hoc)) {
      if (grepl("_summary$", name)) {
        write.csv(association_results$post_hoc[[name]],
                  file = file.path(output_dir, paste0("posthoc_", name, ".csv")),
                  row.names = FALSE)
      }
    }
  }

  if (verbose) {
    cat("\n=== Significant associations (p.adj < 0.01) ===\n")
    sig_results <- all_results %>% filter(p_adjusted < 0.01)
    if (nrow(sig_results) > 0) {
      print(sig_results)
    } else {
      cat("No significant associations found after multiple testing correction.\n")
    }
  }

  return(association_results)
}

#' Detect variable types in clinical data
#'
#' @param data Data frame with clinical variables
#' @param max_unique_categorical Maximum unique values for categorical (default: 10)
#' @return List with categorical and continuous variable names
detect_variable_types <- function(data, max_unique_categorical = 10) {
  categorical <- c()
  continuous <- c()

  for (col in names(data)) {
    # Skip if all NA
    if (all(is.na(data[[col]]))) {
      next
    }

    # Check if numeric
    if (is.numeric(data[[col]])) {
      n_unique <- length(unique(na.omit(data[[col]])))

      # If few unique values, treat as categorical
      if (n_unique <= max_unique_categorical) {
        categorical <- c(categorical, col)
      } else {
        continuous <- c(continuous, col)
      }
    } else {
      # Non-numeric are categorical
      categorical <- c(categorical, col)
    }
  }

  return(list(categorical = categorical, continuous = continuous))
}

#' Perform appropriate test for categorical variables
#'
#' @param cont_table Contingency table
#' @return List with test results
perform_categorical_test <- function(cont_table) {
  # Check expected frequencies
  expected <- chisq.test(cont_table, simulate.p.value = FALSE)$expected

  # Use Fisher's exact test if any expected frequency < 5
  if (any(expected < 5) || sum(cont_table) < 100) {
    # For large tables, use chi-square with simulation
    if (prod(dim(cont_table)) > 6) {
      test_result <- chisq.test(cont_table, simulate.p.value = TRUE, B = 10000)
      test_name <- "Chi-square (simulated)"
    } else {
      test_result <- fisher.test(cont_table, simulate.p.value = TRUE, B = 10000)
      test_name <- "Fisher's exact"
    }
  } else {
    test_result <- chisq.test(cont_table)
    test_name <- "Chi-square"
  }

  return(list(
    test = test_name,
    statistic = ifelse(is.null(test_result$statistic), NA, as.numeric(test_result$statistic)),
    df = ifelse(is.null(test_result$parameter), NA, as.numeric(test_result$parameter)),
    p_value = test_result$p.value
  ))
}

#' Calculate eta-squared effect size for Kruskal-Wallis test
#'
#' @param values Numeric vector of values
#' @param groups Factor vector of groups
#' @return Eta-squared value
calculate_eta_squared_kw <- function(values, groups) {
  kw_stat <- kruskal.test(values ~ groups)$statistic
  n <- length(values)
  k <- length(unique(groups))

  eta_squared <- (kw_stat - k + 1) / (n - k)
  eta_squared <- max(0, eta_squared)  # Ensure non-negative

  return(eta_squared)
}

#' Create visualization for clinical associations
#'
#' @param clinical_data Clinical data with clusters
#' @param association_results Results from perform_clinical_association
#' @param cluster_col Name of cluster column
#' @param output_dir Directory to save plots
#' @param top_n Number of top associations to plot
#' @return NULL (plots are saved)
plot_clinical_associations <- function(clinical_data,
                                       association_results,
                                       cluster_col = "snf_cluster",
                                       output_dir,
                                       top_n = 10) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Convert cluster to factor for consistent plotting
  clinical_data[[cluster_col]] <- as.factor(clinical_data[[cluster_col]])

  # Get top significant associations
  top_vars <- association_results$summary %>%
    filter(p_adjusted < 0.1) %>%
    head(top_n) %>%
    pull(variable)

  if (length(top_vars) == 0) {
    warning("No significant associations to plot")
    return(NULL)
  }

  # Plot each variable
  for (var in top_vars) {
    var_info <- association_results$summary %>%
      filter(variable == var)

    # Skip if variable not in data
    if (!var %in% names(clinical_data)) {
      next
    }

    if (var_info$type == "continuous") {
      # Boxplot for continuous variables
      p <- ggplot(clinical_data %>% filter(!is.na(!!sym(var))),
                  aes(x = !!sym(cluster_col), y = !!sym(var), fill = !!sym(cluster_col))) +
        geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
        stat_compare_means(method = "kruskal.test", label.y.npc = 0.95) +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(
          title = paste(get_var_label(var), "across SNF Clusters"),
          subtitle = paste("Kruskal-Wallis p =",
                           format(var_info$p_value, digits = 3),
                           "| Adjusted p =",
                           format(var_info$p_adjusted, digits = 3)),
          x = "SNF Cluster",
          y = get_var_label(var)
        ) +
        scale_fill_brewer(palette = "Set1")

    } else {
      # Stacked barplot for categorical variables
      # Create proportion table
      prop_table <- clinical_data %>%
        filter(!is.na(!!sym(var))) %>%
        group_by(!!sym(cluster_col), !!sym(var)) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(!!sym(cluster_col)) %>%
        mutate(prop = n / sum(n))

      p <- ggplot(prop_table,
                  aes(x = !!sym(cluster_col), y = prop, fill = !!sym(var))) +
        geom_bar(stat = "identity", position = "stack") +
        theme_minimal() +
        labs(
          title = paste(get_var_label(var), "Distribution across SNF Clusters"),
          subtitle = paste(var_info$test, "p =",
                           format(var_info$p_value, digits = 3),
                           "| Adjusted p =",
                           format(var_info$p_adjusted, digits = 3)),
          x = "SNF Cluster",
          y = "Proportion",
          fill = var
        ) +
        scale_fill_brewer(palette = "Set2") +
        scale_y_continuous(labels = scales::percent)
    }

    # Save plot
    ggsave(filename = file.path(output_dir, paste0("clinical_", var, ".pdf")),
           plot = p,
           width = 8,
           height = 6,
           device = "pdf")
  }

  # Create summary heatmap of p-values
  create_association_heatmap(association_results, output_dir)
}

#' Create heatmap of association p-values
#'
#' @param association_results Results from perform_clinical_association
#' @param output_dir Directory to save plot
#' @return NULL
create_association_heatmap <- function(association_results, output_dir) {

  # Prepare data for heatmap
  p_values <- association_results$summary %>%
    select(variable, p_adjusted) %>%
    mutate(
      neg_log_p = -log10(p_adjusted),
      neg_log_p = pmin(neg_log_p, 10)  # Cap at 10 for visualization
    )

  # Create bar plot of -log10(p-values)
  p <- ggplot(p_values, aes(x = reorder(variable, neg_log_p), y = neg_log_p)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red") +
    coord_flip() +
    theme_minimal() +
    labs(
      title = "Clinical Variable Associations with SNF Clusters",
      x = "Clinical Variable",
      y = "-log10(Adjusted P-value)",
      caption = "Red line indicates p.adj = 0.01"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold")
    )

  ggsave(filename = file.path(output_dir, "clinical_associations_summary_plot.pdf"),
         plot = p,
         width = 8,
         height = 6 + nrow(p_values) * 0.2,
         device = "pdf",
         limitsize = FALSE)
}