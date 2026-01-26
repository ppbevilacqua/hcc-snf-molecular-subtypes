# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tidyverse)
  library(enrichplot)
})

#' Perform GO and KEGG enrichment analysis for cluster DEGs
#'
#' @param deg_results Results from perform_cluster_deg_analysis_log2
#' @param organism Organism for enrichment analysis (default: "hsa" for human)
#' @param pval_cutoff P-value cutoff for enrichment (default: 0.01)
#' @param qval_cutoff Q-value cutoff for enrichment (default: 0.1)
#' @param min_gs_size Minimum gene set size (default: 10)
#' @param max_gs_size Maximum gene set size (default: 500)
#' @param output_dir Directory to save results (default: NULL)
#' @param verbose Print progress messages (default: TRUE)
#' @return List of enrichment results for each cluster
perform_enrichment_analysis <- function(deg_results,
                                       organism = "hsa",
                                       pval_cutoff = 0.01,
                                       qval_cutoff = 0.1,
                                       min_gs_size = 10,
                                       max_gs_size = 500,
                                       output_dir = NULL,
                                       verbose = TRUE) {

  # Set seed for reproducibility
  set.seed(123)

  # Initialize results list
  enrichment_results <- list()

  # Process each cluster
  for (cluster_name in names(deg_results)) {

    if (verbose) {
      cat("\n--- Enrichment Analysis for", cluster_name, "---\n")
    }

    # Get DEGs for this cluster
    if (!is.null(deg_results[[cluster_name]]$filtered)) {
      degs <- deg_results[[cluster_name]]$filtered
    } else {
      degs <- deg_results[[cluster_name]]$unfiltered %>%
        filter(adj.P.Val < 0.01, abs(logFC) > 1)
    }

    if (nrow(degs) < 5) {
      warning(paste("Too few DEGs for", cluster_name, "- skipping enrichment analysis"))
      next
    }

    # Separate up and down regulated genes
    up_genes <- degs %>% filter(direction == "UP") %>% pull(gene)
    down_genes <- degs %>% filter(direction == "DOWN") %>% pull(gene)
    all_degs <- degs$gene

    extract_symbol <- function(gene_ids) {
      sapply(strsplit(gene_ids, "\\|"), function(x) x[length(x)])
    }

    up_symbols <- extract_symbol(up_genes)
    down_symbols <- extract_symbol(down_genes)

    if (verbose) {
      cat("Analyzing", length(all_degs), "DEGs\n")
      cat("  - Upregulated:", length(up_genes), "\n")
      cat("  - Downregulated:", length(down_genes), "\n")
    }

    # Convert gene symbols to Entrez IDs
    gene_mapping <- convert_symbols_to_entrez(all_degs)

    up_entrez <- gene_mapping$entrez[gene_mapping$symbol %in% up_symbols]
    down_entrez <- gene_mapping$entrez[gene_mapping$symbol %in% down_symbols]
    all_entrez <- gene_mapping$entrez

    # Initialize cluster results
    cluster_enrichment <- list()

    # GO Biological Process enrichment
    if (length(all_entrez) >= 5) {
      if (verbose) cat("Running GO Biological Process enrichment...\n")

      go_bp <- enrichGO(
        gene = all_entrez,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = pval_cutoff,
        qvalueCutoff = qval_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        readable = TRUE
      )

      if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
        cluster_enrichment$GO_BP <- go_bp
        if (verbose) cat("  Found", nrow(go_bp@result), "enriched GO BP terms\n")
      }
    }

    # GO enrichment for upregulated genes
    if (length(up_entrez) >= 5) {
      if (verbose) cat("Running GO BP enrichment for upregulated genes...\n")

      go_bp_up <- enrichGO(
        gene = up_entrez,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = pval_cutoff,
        qvalueCutoff = qval_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        readable = TRUE
      )

      if (!is.null(go_bp_up) && nrow(go_bp_up@result) > 0) {
        cluster_enrichment$GO_BP_UP <- go_bp_up
        if (verbose) cat("  Found", nrow(go_bp_up@result), "enriched GO BP terms (UP)\n")
      }
    }

    # GO enrichment for downregulated genes
    if (length(down_entrez) >= 5) {
      if (verbose) cat("Running GO BP enrichment for downregulated genes...\n")

      go_bp_down <- enrichGO(
        gene = down_entrez,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = pval_cutoff,
        qvalueCutoff = qval_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        readable = TRUE
      )

      if (!is.null(go_bp_down) && nrow(go_bp_down@result) > 0) {
        cluster_enrichment$GO_BP_DOWN <- go_bp_down
        if (verbose) cat("  Found", nrow(go_bp_down@result), "enriched GO BP terms (DOWN)\n")
      }
    }

    # KEGG pathway enrichment
    if (length(all_entrez) >= 5) {
      if (verbose) cat("Running KEGG pathway enrichment...\n")

      kegg <- enrichKEGG(
        gene = all_entrez,
        organism = organism,
        pAdjustMethod = "BH",
        pvalueCutoff = pval_cutoff,
        qvalueCutoff = qval_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size
      )

      if (!is.null(kegg) && nrow(kegg@result) > 0) {
        # Convert Entrez IDs to symbols in KEGG results
        kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        cluster_enrichment$KEGG <- kegg
        if (verbose) cat("  Found", nrow(kegg@result), "enriched KEGG pathways\n")
      }
    }

    # KEGG for upregulated genes
    if (length(up_entrez) >= 5) {
      if (verbose) cat("Running KEGG pathway enrichment for upregulated genes...\n")

      kegg_up <- enrichKEGG(
        gene = up_entrez,
        organism = organism,
        pAdjustMethod = "BH",
        pvalueCutoff = pval_cutoff,
        qvalueCutoff = qval_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size
      )

      if (!is.null(kegg_up) && nrow(kegg_up@result) > 0) {
        kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        cluster_enrichment$KEGG_UP <- kegg_up
        if (verbose) cat("  Found", nrow(kegg_up@result), "enriched KEGG pathways (UP)\n")
      }
    }

    # KEGG for downregulated genes
    if (length(down_entrez) >= 5) {
      if (verbose) cat("Running KEGG pathway enrichment for downregulated genes...\n")

      kegg_down <- enrichKEGG(
        gene = down_entrez,
        organism = organism,
        pAdjustMethod = "BH",
        pvalueCutoff = pval_cutoff,
        qvalueCutoff = qval_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size
      )

      if (!is.null(kegg_down) && nrow(kegg_down@result) > 0) {
        kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        cluster_enrichment$KEGG_DOWN <- kegg_down
        if (verbose) cat("  Found", nrow(kegg_down@result), "enriched KEGG pathways (DOWN)\n")
      }
    }

    # Store results
    enrichment_results[[cluster_name]] <- cluster_enrichment

    # Save results if output directory provided
    if (!is.null(output_dir) && length(cluster_enrichment) > 0) {
      cluster_dir <- file.path(output_dir, cluster_name)
      if (!dir.exists(cluster_dir)) {
        dir.create(cluster_dir, recursive = TRUE)
      }

      # Save each enrichment result as CSV
      for (enrich_type in names(cluster_enrichment)) {
        result_df <- as.data.frame(cluster_enrichment[[enrich_type]])
        if (nrow(result_df) > 0) {
          write.csv(result_df,
                   file = file.path(cluster_dir, paste0(enrich_type, "_enrichment.csv")),
                   row.names = FALSE)
        }
      }
    }
  }

  return(enrichment_results)
}

#' Convert gene symbols to Entrez IDs
#'
#' @param gene_symbols Vector of gene symbols
#' @return Data frame with symbol and entrez columns
convert_symbols_to_entrez <- function(gene_symbols) {

  # Remove version numbers if present (e.g., "GENE.1" -> "GENE")
  gene_symbols_clean <- gsub("\\.\\d+$", "", gene_symbols)

  # If genes have format "SYMBOL|ENSEMBL", extract symbol
  if (any(grepl("\\|", gene_symbols_clean))) {
    gene_symbols_clean <- sapply(strsplit(gene_symbols_clean, "\\|"), `[`, 2)
  }

  # Convert to Entrez IDs
  gene_mapping <- bitr(
    gene_symbols_clean,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  # See which genes couldn't be mapped
  all_genes <- gene_symbols_clean
  mapped_genes <- gene_mapping$symbol
  unmapped_genes <- setdiff(all_genes, mapped_genes)
  cat("Number of unmapped genes:", length(unmapped_genes), "\n")
  cat("Examples of unmapped genes:", head(unmapped_genes, 10), "\n")

  colnames(gene_mapping) <- c("symbol", "entrez")

  return(gene_mapping)
}

#' Create enrichment plots for each cluster
#'
#' @param enrichment_results Results from perform_enrichment_analysis
#' @param output_dir Directory to save plots
#' @param plot_type Type of plot: "bar", "dot", or "both" (default: "both")
#' @param top_n Number of top terms to show (default: 10)
#'
#' @return NULL (plots are saved to files)
plot_enrichment_results <- function(enrichment_results,
                                  output_dir,
                                  plot_type = "both",
                                  top_n = 10) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  for (cluster_name in names(enrichment_results)) {
    cluster_enrichment <- enrichment_results[[cluster_name]]

    if (length(cluster_enrichment) == 0) {
      next
    }

    cluster_dir <- file.path(output_dir, cluster_name)
    if (!dir.exists(cluster_dir)) {
      dir.create(cluster_dir, recursive = TRUE)
    }

    # Plot each enrichment type
    for (enrich_type in names(cluster_enrichment)) {
      enrich_obj <- cluster_enrichment[[enrich_type]]

      if (is.null(enrich_obj) || nrow(enrich_obj@result) == 0) {
        next
      }

      # Barplot
      if (plot_type %in% c("bar", "both")) {
        # Check if there are enriched terms to plot
        if (!is.null(enrich_obj) && nrow(enrich_obj@result) > 0) {
          tryCatch({
            p_bar <- barplot(enrich_obj,
                            showCategory = min(top_n, nrow(enrich_obj@result)),
                            title = paste(cluster_name, "-", enrich_type, "Enrichment"))

            ggsave(filename = file.path(cluster_dir, paste0(enrich_type, "_barplot.pdf")),
                   plot = p_bar,
                   width = 10,
                   height = 8,
                   device = "pdf")
          }, error = function(e) {
            message(paste("Could not create barplot for", cluster_name, "-", enrich_type))
          })
        }
      }

      # Dotplot
      if (plot_type %in% c("dot", "both")) {
        # Check if there are enriched terms to plot
        if (!is.null(enrich_obj) && nrow(enrich_obj@result) > 0) {
          tryCatch({
            p_dot <- dotplot(enrich_obj,
                            showCategory = min(top_n, nrow(enrich_obj@result)),
                            title = paste(cluster_name, "-", enrich_type, "Enrichment"))

            ggsave(filename = file.path(cluster_dir, paste0(enrich_type, "_dotplot.pdf")),
                   plot = p_dot,
                   width = 10,
                   height = 8,
                   device = "pdf")
          }, error = function(e) {
            message(paste("Could not create dotplot for", cluster_name, "-", enrich_type))
          })
        }
      }

      # Additional plots for GO terms
      if (grepl("GO", enrich_type)) {
        # Gene-Concept Network
        if (nrow(enrich_obj@result) >= 3) {
          tryCatch({
            p_cnet <- cnetplot(enrich_obj,
                              showCategory = min(5, nrow(enrich_obj@result)),
                              foldChange = NULL,
                              circular = FALSE,
                              colorEdge = TRUE)

            ggsave(filename = file.path(cluster_dir, paste0(enrich_type, "_cnetplot.pdf")),
                   plot = p_cnet,
                   width = 12,
                   height = 10,
                   device = "pdf")
          }, error = function(e) {
            message("Could not create cnetplot for ", enrich_type)
          })
        }
      }
    }
  }
}

#' Create a combined enrichment summary plot across all clusters
#'
#' @param enrichment_results Results from perform_enrichment_analysis
#' @param output_file Path to save the summary plot
#' @param top_n Number of top terms per cluster to include
#' @param enrichment_type Type of enrichment to plot ("GO_BP", "KEGG", etc.)
#' @return ggplot object
create_enrichment_summary_plot <- function(enrichment_results,
                                          output_file = NULL,
                                          top_n = 5,
                                          enrichment_type = "GO_BP") {

  # Collect top enriched terms from each cluster
  summary_data <- data.frame()

  for (cluster_name in names(enrichment_results)) {
    cluster_enrichment <- enrichment_results[[cluster_name]]

    if (enrichment_type %in% names(cluster_enrichment)) {
      enrich_obj <- cluster_enrichment[[enrichment_type]]

      if (!is.null(enrich_obj) && nrow(enrich_obj@result) > 0) {
        top_terms <- head(enrich_obj@result, top_n)
        top_terms$cluster <- gsub("cluster_", "Cluster ", cluster_name)
        top_terms <- top_terms %>%
          select(cluster, Description, p.adjust, Count)
        summary_data <- rbind(summary_data, top_terms)
      }
    }
  }

  if (nrow(summary_data) == 0) {
    warning("No enrichment data found for summary plot")
    return(NULL)
  }

  # Create heatmap-style plot
  p <- ggplot(summary_data, aes(x = cluster, y = Description, fill = -log10(p.adjust))) +
    geom_tile(color = "white") +
    geom_text(aes(label = Count), color = "black", size = 3) +
    scale_fill_gradient(low = "white", high = "darkred",
                       name = "-log10(adj.P)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    labs(
      title = paste("Top", top_n, enrichment_type, "Terms per Cluster"),
      x = "Cluster",
      y = "Enriched Term",
      caption = "Numbers show gene count"
    )

  if (!is.null(output_file)) {
    ggsave(filename = output_file,
           plot = p,
           width = 10 + length(unique(summary_data$cluster)) * 0.5,
           height = 6 + top_n * 0.3,
           device = "pdf",
           limitsize = FALSE)
  }

  return(p)
}

#' Perform Gene Set Enrichment Analysis (GSEA) for ranked genes
#'
#' @param deg_results Results from DEG analysis
#' @param output_dir Output directory for results
#' @param organism Organism identifier (default: "hsa")
#' @param nPerm Number of permutations (default: 10000)
#' @param verbose Print progress messages
#' @return List of GSEA results
perform_gsea_analysis <- function(deg_results,
                                 output_dir = NULL,
                                 organism = "hsa",
                                 nPerm = 10000,
                                 verbose = TRUE) {

  # Set seed for reproducibility
  set.seed(123)

  gsea_results <- list()

  for (cluster_name in names(deg_results)) {
    if (verbose) {
      cat("\n--- GSEA Analysis for", cluster_name, "---\n")
    }

    # Get DEG data
    deg_data <- deg_results[[cluster_name]]$unfiltered

    # Create ranked gene list
    gene_list <- deg_data$logFC
    names(gene_list) <- deg_data$gene

    # Convert symbols to Entrez IDs
    gene_mapping <- convert_symbols_to_entrez(names(gene_list))

    # Create Entrez-based gene list
    gene_list_entrez <- gene_list[gene_mapping$symbol]
    names(gene_list_entrez) <- gene_mapping$entrez

    # Sort by decreasing order
    gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

    # Remove duplicates if any
    gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]

    cluster_gsea <- list()

    # GSEA for GO Biological Process
    if (verbose) cat("Running GSEA for GO BP...\n")

    gsea_go <- gseGO(
      geneList = gene_list_entrez,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      nPerm = nPerm,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.01,
      verbose = FALSE
    )

    if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
      gsea_go <- setReadable(gsea_go, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      cluster_gsea$GSEA_GO_BP <- gsea_go
      if (verbose) cat("  Found", nrow(gsea_go@result), "enriched GO BP gene sets\n")
    }

    # GSEA for KEGG pathways
    if (verbose) cat("Running GSEA for KEGG pathways...\n")

    gsea_kegg <- gseKEGG(
      geneList = gene_list_entrez,
      organism = organism,
      nPerm = nPerm,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.01,
      verbose = FALSE
    )

    if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
      gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      cluster_gsea$GSEA_KEGG <- gsea_kegg
      if (verbose) cat("  Found", nrow(gsea_kegg@result), "enriched KEGG pathways\n")
    }

    gsea_results[[cluster_name]] <- cluster_gsea

    # Save GSEA results
    if (!is.null(output_dir) && length(cluster_gsea) > 0) {
      cluster_dir <- file.path(output_dir, cluster_name)
      if (!dir.exists(cluster_dir)) {
        dir.create(cluster_dir, recursive = TRUE)
      }

      for (gsea_type in names(cluster_gsea)) {
        result_df <- as.data.frame(cluster_gsea[[gsea_type]])
        if (nrow(result_df) > 0) {
          write.csv(result_df,
                   file = file.path(cluster_dir, paste0(gsea_type, ".csv")),
                   row.names = FALSE)

          # Create enrichment plots
          # Running score plot for top pathways
          if (nrow(result_df) >= 3) {
            pdf(file.path(cluster_dir, paste0(gsea_type, "_running_score.pdf")),
                width = 10, height = 8)
            for (i in 1:min(3, nrow(result_df))) {
              print(gseaplot2(cluster_gsea[[gsea_type]],
                            geneSetID = i,
                            title = result_df$Description[i]))
            }
            dev.off()
          }
        }
      }
    }
  }

  return(gsea_results)
}