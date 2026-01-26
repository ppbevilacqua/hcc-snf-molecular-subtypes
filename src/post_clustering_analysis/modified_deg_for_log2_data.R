
#' @param expression_matrix_log2 Already log2-transformed expression matrix
#' @param cluster_assignments Named vector of cluster assignments
#' @param adj_pval_threshold Adjusted p-value threshold (default: 0.01)
#' @param logfc_threshold Log fold-change threshold (default: 1)
#' @param filter_degs Boolean to filter DEGs (default: TRUE)
#' @param output_dir Directory to save results (default: NULL)
#' @param verbose Print progress messages (default: TRUE)
#' @return List of data frames with DEGs for each cluster
perform_cluster_deg_analysis_log2 <- function(expression_matrix_log2,
                                              cluster_assignments,
                                              adj_pval_threshold = 0.01,
                                              logfc_threshold = 1,
                                              filter_degs = TRUE,
                                              output_dir = NULL,
                                              verbose = TRUE) {

  # Validate inputs
  if (!is.matrix(expression_matrix_log2) && !is.data.frame(expression_matrix_log2)) {
    stop("expression_matrix_log2 must be a matrix or data frame")
  }

  # Convert to matrix if needed
  if (is.data.frame(expression_matrix_log2)) {
    expression_matrix_log2 <- as.matrix(expression_matrix_log2)
  }

  # IMPORTANT: Verify that data appears to be log2-transformed
  data_range <- range(expression_matrix_log2, na.rm = TRUE)
  if (data_range[2] > 100) {
    warning("Data appears to NOT be log2-transformed (max value > 100). Please verify!")
  }

  if (verbose) {
    cat("Data range: [", data_range[1], ",", data_range[2], "]\n")
    cat("Assuming data is already log2-transformed\n\n")
  }

  # Check patient ID matching
  common_patients <- intersect(names(cluster_assignments), colnames(expression_matrix_log2))
  if (length(common_patients) == 0) {
    stop("No matching patient IDs between expression matrix and cluster assignments")
  }

  if (verbose) {
    cat("Found", length(common_patients), "patients with both expression and cluster data\n")
  }

  # Subset and align data
  expression_matrix_log2 <- expression_matrix_log2[, common_patients]
  cluster_assignments <- cluster_assignments[common_patients]

  # Get unique clusters
  unique_clusters <- sort(unique(cluster_assignments))

  if (verbose) {
    cat("Analyzing", length(unique_clusters), "clusters:", paste(unique_clusters, collapse = ", "), "\n")
  }

  # Initialize results list
  deg_results <- list()

  # Perform one-vs-all comparison for each cluster
  for (cluster_id in unique_clusters) {

    if (verbose) {
      cat("\n--- Analyzing Cluster", cluster_id, "---\n")
    }

    # Create binary design for one-vs-all comparison
    cluster_binary <- ifelse(cluster_assignments == cluster_id, 1, 0)

    # Count samples per group
    n_cluster <- sum(cluster_binary == 1)
    n_others <- sum(cluster_binary == 0)

    if (verbose) {
      cat("Cluster", cluster_id, "samples:", n_cluster, "\n")
      cat("Other clusters samples:", n_others, "\n")
    }

    # Check minimum sample size
    if (n_cluster < 3 || n_others < 3) {
      warning(paste("Cluster", cluster_id, "or comparison group has fewer than 3 samples. Skipping."))
      next
    }

    # Create design matrix
    design <- model.matrix(~cluster_binary)
    colnames(design) <- c("Intercept", paste0("Cluster", cluster_id, "_vs_Others"))

    # IMPORTANT: Data is already log2-transformed, so use it directly
    # Fit linear model using limma on log2 data
    fit <- lmFit(expression_matrix_log2, design)
    fit <- eBayes(fit)

    # Extract results for the cluster comparison
    deg_table <- topTable(fit,
                          coef = 2,  # Compare cluster vs others
                          number = Inf,
                          sort.by = "P")

    # Add gene names as a column
    deg_table$gene <- rownames(deg_table)

    # Calculate average expression for cluster and others
    cluster_samples <- names(cluster_assignments)[cluster_assignments == cluster_id]
    other_samples <- names(cluster_assignments)[cluster_assignments != cluster_id]

    deg_table$avg_expr_cluster <- rowMeans(expression_matrix_log2[deg_table$gene, cluster_samples, drop = FALSE])
    deg_table$avg_expr_others <- rowMeans(expression_matrix_log2[deg_table$gene, other_samples, drop = FALSE])

    # Add direction of regulation
    deg_table$direction <- ifelse(deg_table$logFC > 0, "UP", "DOWN")

    # Reorganize columns
    deg_table <- deg_table %>%
      dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val, B,
                    avg_expr_cluster, avg_expr_others, direction) %>%
      dplyr::arrange(adj.P.Val, desc(abs(logFC)))

    # Apply filtering if requested
    if (filter_degs) {
      deg_table_filtered <- deg_table %>%
        filter(adj.P.Val < adj_pval_threshold,
               abs(logFC) > logfc_threshold)

      n_up <- sum(deg_table_filtered$direction == "UP")
      n_down <- sum(deg_table_filtered$direction == "DOWN")

      if (verbose) {
        cat("DEGs found (filtered):", nrow(deg_table_filtered), "\n")
        cat("  - Upregulated:", n_up, "\n")
        cat("  - Downregulated:", n_down, "\n")
      }

      # Store both filtered and unfiltered results
      deg_results[[paste0("cluster_", cluster_id)]] <- list(
        filtered = deg_table_filtered,
        unfiltered = deg_table,
        summary = data.frame(
          cluster = cluster_id,
          n_samples = n_cluster,
          n_degs_total = nrow(deg_table),
          n_degs_filtered = nrow(deg_table_filtered),
          n_up = n_up,
          n_down = n_down,
          adj_pval_threshold = adj_pval_threshold,
          logfc_threshold = logfc_threshold
        )
      )
    } else {
      n_up <- sum(deg_table$direction == "UP" &
                    deg_table$adj.P.Val < adj_pval_threshold &
                    abs(deg_table$logFC) > logfc_threshold)
      n_down <- sum(deg_table$direction == "DOWN" &
                      deg_table$adj.P.Val < adj_pval_threshold &
                      abs(deg_table$logFC) > logfc_threshold)

      if (verbose) {
        cat("Total genes analyzed:", nrow(deg_table), "\n")
        cat("Significant DEGs (at thresholds):",
            sum(deg_table$adj.P.Val < adj_pval_threshold & abs(deg_table$logFC) > logfc_threshold), "\n")
      }

      deg_results[[paste0("cluster_", cluster_id)]] <- list(
        filtered = NULL,
        unfiltered = deg_table,
        summary = data.frame(
          cluster = cluster_id,
          n_samples = n_cluster,
          n_degs_total = nrow(deg_table),
          n_degs_filtered = sum(deg_table$adj.P.Val < adj_pval_threshold &
                                  abs(deg_table$logFC) > logfc_threshold),
          n_up = n_up,
          n_down = n_down,
          adj_pval_threshold = adj_pval_threshold,
          logfc_threshold = logfc_threshold
        )
      )
    }

    # Save results if output directory is provided
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }

      # Save filtered results (if filtering was applied)
      if (filter_degs) {
        write.csv(deg_table_filtered,
                  file = file.path(output_dir, paste0("DEGs_cluster_", cluster_id, "_filtered.csv")),
                  row.names = FALSE)
      }

      # Save unfiltered results
      write.csv(deg_table,
                file = file.path(output_dir, paste0("DEGs_cluster_", cluster_id, "_all.csv")),
                row.names = FALSE)
    }
  }

  # Create and save summary table
  if (length(deg_results) > 0) {
    summary_table <- bind_rows(lapply(deg_results, function(x) x$summary))

    if (verbose) {
      cat("\n=== Overall Summary ===\n")
      print(summary_table)
    }

    if (!is.null(output_dir)) {
      write.csv(summary_table,
                file = file.path(output_dir, "DEG_analysis_summary.csv"),
                row.names = FALSE)
    }
  }

  return(deg_results)
}

run_additional_cluster_analysis_log2 <- function(expression_matrix_log2,
                                                 clinical_data,
                                                 output_dir,
                                                 cluster_col = "snf_cluster",
                                                 run_deg = TRUE,
                                                 run_enrichment = TRUE,
                                                 run_clinical = TRUE,
                                                 deg_params = list(),
                                                 enrichment_params = list(),
                                                 clinical_params = list(),
                                                 verbose = TRUE) {

  # Create main output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Create timestamp for this run
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # Log file
  log_file <- file.path(output_dir, paste0("analysis_log_", timestamp, ".txt"))

  # Start logging
  if (verbose) {
    sink(log_file, split = TRUE)
    cat("===========================================\n")
    cat("SNF Post-Clustering Analysis Pipeline\n")
    cat("Using LOG2-TRANSFORMED Data\n")
    cat("===========================================\n")
    cat("Start time:", as.character(Sys.time()), "\n")
    cat("Output directory:", output_dir, "\n\n")
  }

  # Validate inputs
  if (!cluster_col %in% names(clinical_data)) {
    stop(paste("Cluster column", cluster_col, "not found in clinical data"))
  }

  # Extract cluster assignments
  cluster_assignments <- setNames(
    clinical_data[[cluster_col]],
    rownames(clinical_data)
  )

  # Check for patient ID matching
  common_patients <- intersect(names(cluster_assignments), colnames(expression_matrix_log2))

  if (verbose) {
    cat("Data Summary:\n")
    cat("- Genes in expression matrix:", nrow(expression_matrix_log2), "\n")
    cat("- Total patients in expression matrix:", ncol(expression_matrix_log2), "\n")
    cat("- Total patients in clinical data:", nrow(clinical_data), "\n")
    cat("- Matched patients:", length(common_patients), "\n")
    cat("- Number of clusters:", length(unique(cluster_assignments)), "\n")
    cat("- Data is LOG2-TRANSFORMED\n")
    cat("- Cluster distribution:\n")
    print(table(cluster_assignments[common_patients]))
    cat("\n")
  }

  # Initialize results container
  all_results <- list(
    metadata = list(
      timestamp = timestamp,
      n_genes = nrow(expression_matrix_log2),
      n_patients = length(common_patients),
      n_clusters = length(unique(cluster_assignments)),
      cluster_sizes = table(cluster_assignments[common_patients]),
      data_type = "log2-transformed"
    )
  )

  # ============================================================
  # 1. Differential Expression Analysis (Modified for log2 data)
  # ============================================================

  if (run_deg) {
    if (verbose) {
      cat("\n===========================================\n")
      cat("STEP 1: Differential Expression Analysis\n")
      cat("        (Using log2-transformed data)\n")
      cat("===========================================\n")
    }

    # Set default DEG parameters
    deg_params_default <- list(
      adj_pval_threshold = 0.01,
      logfc_threshold = 1,
      filter_degs = TRUE,
      output_dir = file.path(output_dir, "DEG_analysis"),
      verbose = verbose
    )

    # Merge with user parameters
    deg_params <- modifyList(deg_params_default, deg_params)

    # Run modified DEG analysis for log2 data
    deg_results <- do.call(perform_cluster_deg_analysis_log2, c(
      list(
        expression_matrix_log2 = expression_matrix_log2,
        cluster_assignments = cluster_assignments
      ),
      deg_params
    ))

    top_10_degs_per_cluster <- list()

    for (cluster_name in names(deg_results)) {
      # Get filtered DEGs for this cluster
      filtered_degs <- deg_results[[cluster_name]]$filtered

      if (nrow(filtered_degs) > 0) {
        # Sort by adjusted p-value and select top 10
        top_10 <- filtered_degs %>%
          arrange(adj.P.Val) %>%
          head(10) %>%
          mutate(cluster = cluster_name, rank = 1:n()) %>%
          select(cluster, rank, gene, logFC, adj.P.Val, direction, AveExpr)

        top_10_degs_per_cluster[[cluster_name]] <- top_10

      }
    }

    # Combine all clusters into single data frame
    top_10_combined <- bind_rows(top_10_degs_per_cluster)
    write.csv(top_10_combined,  file = file.path(deg_params$output_dir, "top_10_expression.csv"))

    # Create volcano plots
    if (length(deg_results) > 0) {
      if (verbose) cat("\nCreating volcano plots...\n")

      plot_deg_volcano(
        deg_results = deg_results,
        output_dir = deg_params$output_dir,
        adj_pval_threshold = deg_params$adj_pval_threshold,
        logfc_threshold = deg_params$logfc_threshold
      )

      # Create heatmap of top DEGs
      if (verbose) cat("Creating DEG heatmap...\n")

      plot_deg_heatmap(
        expression_matrix = expression_matrix_log2,
        deg_results = deg_results,
        cluster_assignments = cluster_assignments,
        top_n_genes = 50,
        output_file = file.path(deg_params$output_dir, "top_DEGs_heatmap.pdf")
      )
    }

    all_results$deg <- deg_results
  }

  # ============================================================
  # 2. Functional Enrichment Analysis
  # ============================================================

  if (run_enrichment &&
    run_deg &&
    length(all_results$deg) > 0) {
    if (verbose) {
      cat("\n===========================================\n")
      cat("STEP 2: Functional Enrichment Analysis\n")
      cat("===========================================\n")
    }

    # Run enrichment analysis
    enrichment_results <- do.call(perform_enrichment_analysis, c(
      list(deg_results = all_results$deg),
      enrichment_params
    ))

    # Create enrichment plots
    if (length(enrichment_results) > 0) {
      if (verbose) cat("\nCreating enrichment plots...\n")

      plot_enrichment_results(
        enrichment_results = enrichment_results,
        output_dir = enrichment_params$output_dir,
        plot_type = "both",
        top_n = 10
      )

      # Create summary plots
      for (enrich_type in c("GO_BP", "KEGG")) {
        create_enrichment_summary_plot(
          enrichment_results = enrichment_results,
          output_file = file.path(enrichment_params$output_dir,
                                  paste0(enrich_type, "_summary_heatmap.pdf")),
          top_n = 5,
          enrichment_type = enrich_type
        )
      }
    }

    all_results$enrichment <- enrichment_results
  }

  # ============================================================
  # 3. Clinical Feature Association Analysis
  # ============================================================

  if (run_clinical) {
    if (verbose) {
      cat("\n===========================================\n")
      cat("STEP 3: Clinical Feature Association Analysis\n")
      cat("===========================================\n")
    }

    # Run clinical association analysis
    clinical_results <- do.call(perform_clinical_association, c(
      list(
        clinical_data = clinical_data
      ),
      clinical_params
    ))

    # Create clinical plots
    if (nrow(clinical_results$summary) > 0) {
      if (verbose) cat("\nCreating clinical association plots...\n")

      plot_clinical_associations(
        clinical_data = clinical_data,
        association_results = clinical_results,
        cluster_col = cluster_col,
        output_dir = clinical_params$output_dir,
        top_n = 10
      )
    }

    all_results$clinical <- clinical_results
  }

  # End logging
  if (verbose) {
    cat("\n===========================================\n")
    cat("Analysis Complete\n")
    cat("End time:", as.character(Sys.time()), "\n")
    cat("===========================================\n")
    sink()
  }

  return(all_results)
}