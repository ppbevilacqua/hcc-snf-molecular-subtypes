#' Main analysis function for TCGA-LIHC data clustered
#'
#' @param mRNA_matrix_cancer log2-transformed mRNA matrix
#' @param miRNA_matrix_cancer log2-transformed miRNA matrix
#' @param final_clinical_matrix clinical data with snf_cluster column
#' @param output_base_dir Base directory for all outputs
#'
#' @return List of all analysis results
run_tcga_lihc_snf_analysis <- function(mRNA_matrix_cancer,
                                       miRNA_matrix_cancer,
                                       final_clinical_matrix,
                                       output_base_dir = "TCGA_LIHC_SNF_Analysis") {

  # Create timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  cat("=====================================\n")
  cat("TCGA-LIHC SNF Post-Clustering Analysis\n")
  cat("=====================================\n")
  cat("Start time:", as.character(Sys.time()), "\n\n")

  # ============================================================
  # Data Validation
  # ============================================================

  # Check for snf_cluster column
  if (!"snf_cluster" %in% names(final_clinical_matrix)) {
    stop("'snf_cluster' column not found in clinical data!")
  }

  mrna_patients <- str_extract(colnames(mRNA_matrix_cancer), PZ_TCGA_REGEX)
  mirna_patients <- str_extract(colnames(miRNA_matrix_cancer), PZ_TCGA_REGEX)
  clinical_patients <- rownames(final_clinical_matrix)

  # Find common patients
  common_patients <- Reduce(intersect, list(mrna_patients, mirna_patients, clinical_patients))

  # Create cluster assignments
  cluster_assignments <- setNames(
    final_clinical_matrix$snf_cluster,
    rownames(final_clinical_matrix)
  )

  cat("SNF Cluster distribution:\n")
  print(table(cluster_assignments))
  cat("\n")

  # ============================================================
  # 1. mRNA Analysis
  # ============================================================

  cat("\n=====================================\n")
  cat("PART 1: mRNA Expression Analysis\n")
  cat("=====================================\n\n")

  mrna_output_dir <- file.path(output_base_dir, "mRNA_analysis")

  mrna_results <- run_additional_cluster_analysis_log2(
    expression_matrix_log2 = mRNA_matrix_cancer,
    clinical_data = final_clinical_matrix,
    output_dir = mrna_output_dir,
    cluster_col = "snf_cluster",
    run_deg = TRUE,
    run_enrichment = TRUE,
    run_clinical = TRUE,
    deg_params = list(
      adj_pval_threshold = 0.01,
      logfc_threshold = 1,
      filter_degs = TRUE,
      output_dir = file.path(mrna_output_dir, "DEG_analysis"),
      verbose = TRUE
    ),
    enrichment_params = list(
      organism = "hsa",
      pval_cutoff = 0.01,
      qval_cutoff = 0.1,
      min_gs_size = 10,
      max_gs_size = 500,
      output_dir = file.path(mrna_output_dir, "Enrichment_analysis"),
      verbose = TRUE
    ),
    clinical_params = list(
      cluster_col = "snf_cluster",
      exclude_cols = c("days_to_diagnosis"),  # Exclude non-informative columns

      # Categorical variables from your data
      categorical_cols = c(
        "obesity_class_2",
        "family_history",
        "alcoholic_liver",
        "hepatitis_c",
        "hepatitis_b",
        "cirrhosis",
        "ajcc_pathologic_stage",
        "prior_malignancy",
        "prior_treatment",
        "ajcc_pathologic_t",
        "ajcc_pathologic_n",
        "ajcc_pathologic_m",
        "gender",
        "treatments_pharmaceutical_treatment_or_therapy",
        "treatments_radiation_treatment_or_therapy"
      ),

      # Continuous variables
      continuous_cols = c(
        "bmi",
        "age_at_index"
      ),
      output_dir = file.path(mrna_output_dir, "Clinical_associations"),
      verbose = TRUE
    ),
    verbose = TRUE
  )

  # ============================================================
  # 2. miRNA Analysis (without enrichment)
  # ============================================================

  cat("\n=====================================\n")
  cat("PART 2: miRNA Expression Analysis\n")
  cat("=====================================\n\n")

  mirna_output_dir <- file.path(output_base_dir, "miRNA_analysis")

  mirna_results <- run_additional_cluster_analysis_log2(
    expression_matrix_log2 = miRNA_matrix_cancer,
    clinical_data = final_clinical_matrix,
    output_dir = mirna_output_dir,
    cluster_col = "snf_cluster",
    run_deg = TRUE,
    run_enrichment = FALSE,  # miRNA enrichment needs different tools
    run_clinical = FALSE,    # Already done with mRNA
    deg_params = list(
      adj_pval_threshold = 0.01,
      logfc_threshold = 0.5,    # Lower threshold for miRNAs
      filter_degs = TRUE,
      output_dir = file.path(mirna_output_dir, "DEM_analysis"),
      verbose = TRUE
    ),
    verbose = TRUE
  )

  # ============================================================
  # 3. Integrated Analysis
  # ============================================================

  cat("\n=====================================\n")
  cat("PART 3: Integrated Analysis\n")
  cat("=====================================\n\n")

  integrated_dir <- file.path(output_base_dir, "integrated_analysis")
  if (!dir.exists(integrated_dir)) {
    dir.create(integrated_dir, recursive = TRUE)
  }

  # Create integrated summary
  integrated_summary <- create_integrated_summary(
    mrna_results = mrna_results,
    mirna_results = mirna_results,
    clinical_results = mrna_results$clinical,
    output_dir = integrated_dir
  )

  # ============================================================
  # 4. Create Final Report
  # ============================================================

  cat("\nGenerating final integrated report...\n")

  final_report <- generate_final_report(
    mrna_results = mrna_results,
    mirna_results = mirna_results,
    integrated_summary = integrated_summary,
    clinical_data = final_clinical_matrix,
    output_dir = output_base_dir,
    timestamp = timestamp
  )

  # Save all results as RDS
  saveRDS(
    list(
      mrna = mrna_results,
      mirna = mirna_results,
      integrated = integrated_summary,
      metadata = list(
        timestamp = timestamp,
        n_patients = length(common_patients),
        n_mrna_genes = nrow(mRNA_matrix_cancer),
        n_mirnas = nrow(miRNA_matrix_cancer)
      )
    ),
    file = file.path(output_base_dir, paste0("all_results_", timestamp, ".rds"))
  )

  cat("\n=====================================\n")
  cat("Analysis Complete!\n")
  cat("End time:", as.character(Sys.time()), "\n")
  cat("Results saved to:", output_base_dir, "\n")
  cat("=====================================\n")

  return(list(
    mrna = mrna_results,
    mirna = mirna_results,
    integrated = integrated_summary
  ))
}

#' Create integrated summary of mRNA and miRNA results
#'
#' @param mrna_results Results from mRNA analysis
#' @param mirna_results Results from miRNA analysis
#' @param clinical_results Results from clinical analysis
#' @param output_dir Output directory
#' @return Integrated summary list
create_integrated_summary <- function(mrna_results,
                                     mirna_results,
                                     clinical_results,
                                     output_dir) {

  cat("Creating integrated summary...\n")

  # Extract DEG summaries
  mrna_deg_summary <- bind_rows(lapply(mrna_results$deg, function(x) x$summary))
  mirna_deg_summary <- bind_rows(lapply(mirna_results$deg, function(x) x$summary))

  # Add data type column
  mrna_deg_summary$data_type <- "mRNA"
  mirna_deg_summary$data_type <- "miRNA"

  # Combine summaries
  combined_deg_summary <- bind_rows(mrna_deg_summary, mirna_deg_summary)

  # Save combined summary
  write.csv(combined_deg_summary,
           file = file.path(output_dir, "integrated_deg_summary.csv"),
           row.names = FALSE)

  # Create visualization of DEGs across clusters
  p_deg_summary <- ggplot(combined_deg_summary,
                          aes(x = factor(cluster), y = n_degs_filtered, fill = data_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Differentially Expressed Features per Cluster",
         x = "SNF Cluster",
         y = "Number of DEGs/DEMs",
         fill = "Data Type") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1")

  ggsave(filename = file.path(output_dir, "integrated_deg_barplot.pdf"),
         plot = p_deg_summary,
         width = 10,
         height = 6)

  # Extract top clinical associations
  if (!is.null(clinical_results) && nrow(clinical_results$summary) > 0) {
    top_clinical <- clinical_results$summary %>%
      filter(p_adjusted < 0.1) %>%
      arrange(p_adjusted)

    write.csv(top_clinical,
             file = file.path(output_dir, "top_clinical_associations.csv"),
             row.names = FALSE)
  }

  return(list(
    deg_summary = combined_deg_summary,
    clinical_summary = clinical_results$summary
  ))
}

#' Generate final integrated report
#'
#' @param mrna_results mRNA analysis results
#' @param mirna_results miRNA analysis results
#' @param integrated_summary Integrated summary
#' @param clinical_data Clinical data
#' @param output_dir Output directory
#' @param timestamp Analysis timestamp
#' @return Report content as character vector
generate_final_report <- function(mrna_results,
                                 mirna_results,
                                 integrated_summary,
                                 clinical_data,
                                 output_dir,
                                 timestamp) {

  report <- c(
    "================================================================================",
    "TCGA-LIHC SNF POST-CLUSTERING ANALYSIS - FINAL REPORT",
    "================================================================================",
    paste("Generated:", Sys.time()),
    paste("Analysis ID:", timestamp),
    "",
    "--------------------------------------------------------------------------------",
    "STUDY OVERVIEW",
    "--------------------------------------------------------------------------------",
    paste("Total Patients Analyzed:", nrow(clinical_data)),
    paste("Number of SNF Clusters:", length(unique(clinical_data$snf_cluster))),
    "",
    "Cluster Distribution:",
    capture.output(print(table(clinical_data$snf_cluster))),
    "",
    "--------------------------------------------------------------------------------",
    "mRNA EXPRESSION ANALYSIS",
    "--------------------------------------------------------------------------------",
    paste("Total Genes Analyzed:", mrna_results$metadata$n_genes),
    ""
  )

  # Add mRNA DEG summary per cluster
  for (cluster_name in names(mrna_results$deg)) {
    summary_info <- mrna_results$deg[[cluster_name]]$summary
    report <- c(report,
      paste0("  ", toupper(gsub("_", " ", cluster_name)), ":"),
      paste("    - Significant DEGs:", summary_info$n_degs_filtered),
      paste("      Upregulated:", summary_info$n_up),
      paste("      Downregulated:", summary_info$n_down)
    )
  }

  report <- c(report,
    "",
    "--------------------------------------------------------------------------------",
    "miRNA EXPRESSION ANALYSIS",
    "--------------------------------------------------------------------------------",
    paste("Total miRNAs Analyzed:", mirna_results$metadata$n_genes),
    ""
  )

  # Add miRNA DEG summary per cluster
  for (cluster_name in names(mirna_results$deg)) {
    summary_info <- mirna_results$deg[[cluster_name]]$summary
    report <- c(report,
      paste0("  ", toupper(gsub("_", " ", cluster_name)), ":"),
      paste("    - Significant DEMs:", summary_info$n_degs_filtered),
      paste("      Upregulated:", summary_info$n_up),
      paste("      Downregulated:", summary_info$n_down)
    )
  }

  # Add clinical associations
  if (!is.null(mrna_results$clinical)) {
    report <- c(report,
      "",
      "--------------------------------------------------------------------------------",
      "CLINICAL ASSOCIATIONS",
      "--------------------------------------------------------------------------------",
      paste("Total Clinical Variables Tested:", nrow(mrna_results$clinical$summary)),
      "",
      "Top Significant Associations (p.adj < 0.01):"
    )

    sig_clinical <- mrna_results$clinical$summary %>%
      filter(p_adjusted < 0.01) %>%
      arrange(p_adjusted) %>%
      head(10)

    if (nrow(sig_clinical) > 0) {
      for (i in 1:nrow(sig_clinical)) {
        report <- c(report,
          paste("  ", i, ". ", sig_clinical$variable[i],
                " (p.adj = ", format(sig_clinical$p_adjusted[i], digits = 3), ")",
                sep = "")
        )
      }
    } else {
      report <- c(report, "  No significant associations at p.adj < 0.01")
    }
  }

  # Add enrichment summary if available
  if (!is.null(mrna_results$enrichment)) {
    report <- c(report,
      "",
      "--------------------------------------------------------------------------------",
      "FUNCTIONAL ENRICHMENT HIGHLIGHTS",
      "--------------------------------------------------------------------------------"
    )

    for (cluster_name in names(mrna_results$enrichment)) {
      cluster_enrich <- mrna_results$enrichment[[cluster_name]]
      if (length(cluster_enrich) > 0) {
        report <- c(report,
          "",
          paste0("  ", toupper(gsub("_", " ", cluster_name)), ":")
        )

        # Add top GO terms
        if ("GO_BP" %in% names(cluster_enrich) &&
            !is.null(cluster_enrich$GO_BP) &&
            nrow(cluster_enrich$GO_BP@result) > 0) {
          top_go <- head(cluster_enrich$GO_BP@result$Description, 3)
          report <- c(report,
            "    Top GO Biological Processes:",
            paste("      -", top_go)
          )
        }

        # Add top KEGG pathways
        if ("KEGG" %in% names(cluster_enrich) &&
            !is.null(cluster_enrich$KEGG) &&
            nrow(cluster_enrich$KEGG@result) > 0) {
          top_kegg <- head(cluster_enrich$KEGG@result$Description, 3)
          report <- c(report,
            "    Top KEGG Pathways:",
            paste("      -", top_kegg)
          )
        }
      }
    }
  }

  report <- c(report,
    "",
    "--------------------------------------------------------------------------------",
    "OUTPUT FILES GENERATED",
    "--------------------------------------------------------------------------------",
    "  1. DEG analysis results (CSV files per cluster)",
    "  2. Enrichment analysis results (CSV files per cluster)",
    "  3. Clinical association results (CSV files)",
    "  4. Visualization plots (PDF files):",
    "     - Volcano plots",
    "     - Heatmaps",
    "     - Enrichment plots",
    "     - Clinical association plots",
    "  5. Integrated summaries",
    "  6. R data objects (.rds files)",
    "",
    "--------------------------------------------------------------------------------",
    "END OF REPORT",
    "--------------------------------------------------------------------------------"
  )

  # Save report
  report_file <- file.path(output_dir, paste0("final_report_", timestamp, ".txt"))
  writeLines(report, report_file)

  cat("Final report saved to:", report_file, "\n")

  return(report)
}

# ============================================================
# EXECUTE THE ANALYSIS
# ============================================================

# This is the main execution block
# Uncomment and run with your actual data objects

# results <- run_tcga_lihc_snf_analysis(
#   mRNA_matrix_cancer = mRNA_matrix_cancer,
#   miRNA_matrix_cancer = miRNA_matrix_cancer,
#   final_clinical_matrix = final_clinical_matrix,
#   output_base_dir = "TCGA_LIHC_SNF_Results"
# )

# To access specific results:
# mrna_deg_cluster1 <- results$mrna$deg$cluster_1$filtered
# mirna_deg_cluster1 <- results$mirna$deg$cluster_1$filtered
# clinical_associations <- results$mrna$clinical$summary
# enrichment_cluster1 <- results$mrna$enrichment$cluster_1