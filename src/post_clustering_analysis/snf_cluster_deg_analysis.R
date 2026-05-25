plot_deg_volcano <- function(deg_results,
                             output_dir,
                             adj_pval_threshold = 0.01,
                             logfc_threshold = 1,
                             top_genes = 10,
                             feature_label = "DEGs") {

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  for (cluster_name in names(deg_results)) {

    # Get DEG data (use unfiltered for complete volcano plot)
    deg_data <- deg_results[[cluster_name]]$unfiltered

    # Add significance categories
    deg_data <- deg_data %>%
      mutate(
        significance = case_when(
          adj.P.Val >= adj_pval_threshold ~ "Not Significant",
          abs(logFC) <= logfc_threshold ~ "Not Significant",
          logFC > logfc_threshold ~ "Upregulated",
          logFC < -logfc_threshold ~ "Downregulated",
          TRUE ~ "Not Significant"
        )
      )

    # Select top genes to label
    top_up <- deg_data %>%
      filter(significance == "Upregulated") %>%
      arrange(adj.P.Val) %>%
      head(top_genes / 2)

    top_down <- deg_data %>%
      filter(significance == "Downregulated") %>%
      arrange(adj.P.Val) %>%
      head(top_genes / 2)

    genes_to_label <- bind_rows(top_up, top_down)

    # Create volcano plot
    p <- ggplot(deg_data, aes(x = logFC, y = -log10(adj.P.Val))) +
      geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Upregulated" = "red",
                                    "Downregulated" = "blue",
                                    "Not Significant" = "gray60")) +
      geom_hline(yintercept = -log10(adj_pval_threshold),
                 linetype = "dashed", color = "gray40", alpha = 0.7) +
      geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
                 linetype = "dashed", color = "gray40", alpha = 0.7) +
      labs(
        title = paste("Volcano Plot -", gsub("_", " ", cluster_name)),
        subtitle = paste("Adj. p-value <", adj_pval_threshold, "& |logFC| >", logfc_threshold),
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value",
        color = paste(sub("s$", "", feature_label), "Status")
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
      )

    # Add gene labels if ggrepel is available
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = genes_to_label,
        aes(label = gene),
        size = 3,
        max.overlaps = 20,
        segment.size = 0.3,
        segment.alpha = 0.5
      )
    }

    # Save plot
    ggsave(filename = file.path(output_dir, paste0("volcano_", cluster_name, ".pdf")),
           plot = p,
           width = 10,
           height = 8,
           device = "pdf")
  }
}

plot_deg_heatmap <- function(expression_matrix,
                             deg_results,
                             cluster_assignments,
                             output_file,
                             top_n_genes = 50,
                             feature_label = "DEGs") {

  # Collect top DEGs from each cluster
  top_degs <- c()
  for (cluster_name in names(deg_results)) {
    if (!is.null(deg_results[[cluster_name]]$filtered)) {
      cluster_degs <- deg_results[[cluster_name]]$filtered
    } else {
      cluster_degs <- deg_results[[cluster_name]]$unfiltered %>%
        filter(adj.P.Val < 0.01, abs(logFC) > 1)
    }

    if (nrow(cluster_degs) > 0) {
      top_genes <- head(cluster_degs$gene, top_n_genes)
      top_degs <- unique(c(top_degs, top_genes))
    }
  }

  if (length(top_degs) == 0) {
    warning(paste("No significant", feature_label, "found for heatmap"))
    return(NULL)
  }

  # Subset expression matrix
  common_patients <- intersect(names(cluster_assignments), colnames(expression_matrix))
  heatmap_data <- expression_matrix[top_degs, common_patients]

  # Create annotation
  annotation_col <- data.frame(
    Cluster = as.factor(cluster_assignments[common_patients]),
    row.names = common_patients
  )

  # Define colors for clusters
  n_clusters <- length(unique(cluster_assignments))
  cluster_colors <- setNames(
    colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))(n_clusters),
    sort(unique(cluster_assignments))
  )

  annotation_colors <- list(Cluster = cluster_colors)

  # Create heatmap
  pdf(output_file, width = 12, height = 10)
  
  pheatmap::pheatmap(
    heatmap_data,
    scale = "row",
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "complete",
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = paste("Top", top_n_genes, feature_label, "per Cluster"),
    border_color = NA,
    color = colorRampPalette(c("blue", "white", "red"))(100)
  )

  dev.off()

}