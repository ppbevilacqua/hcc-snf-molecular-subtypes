network_analysis_SNF <- function(mRNA_matrix_cancer,
                                 miRNA_matrix_cancer,
                                 K,                     # number of nearest neighbors (~20% of sample size)
                                 sigma,                 # sigma in affinity matrix (variance for local model)
                                 iterations,            # number of SNF iterations
                                 n_clusters = 0,
                                 distance = "correlation",
                                 dir_plots = dirPlots) {

  ############################## Data preparation ###############################

  # Transpose omics data to have patients as rows (required for SNF)
  mRNA_data <- t(mRNA_matrix_cancer)
  miRNA_data <- t(miRNA_matrix_cancer)

  rownames(mRNA_data) <- str_extract(rownames(mRNA_data), PZ_TCGA_REGEX)
  rownames(miRNA_data) <- str_extract(rownames(miRNA_data), PZ_TCGA_REGEX)

  cat("\nData dimensions after transpose:\n")
  cat("mRNA:", dim(mRNA_data), "\n")      # patients x genes
  cat("miRNA:", dim(miRNA_data), "\n")    # patients x genes

  ############################ 1. Normalization ############################

  # Standard normalization for each data type
  mRNA_normalized <- standardNormalization(mRNA_data)
  miRNA_normalized <- standardNormalization(miRNA_data)

  ######################## 2. Distance calculation #########################

  # Calculate pairwise distances (euclidean or correlation)
  if (distance == "euclidean") {
    # Using Euclidean distance for omics data
    dist_mRNA <- (dist2(as.matrix(mRNA_normalized), as.matrix(mRNA_normalized)))^(1 / 2)
    dist_miRNA <- (dist2(as.matrix(miRNA_normalized), as.matrix(miRNA_normalized)))^(1 / 2)

  } else if (distance == "correlation") {
    # Using correlation distance for omics data (1 - correlation)
    dist_mRNA <- as.matrix(1 - cor(t(mRNA_normalized), method = "spearman"))
    dist_miRNA <- as.matrix(1 - cor(t(miRNA_normalized), method = "spearman"))

    # Force sigma to be 1 for correlation
    sigma <- 1
  }

  # Force diagonal to be 0
  diag(dist_mRNA) <- 0
  diag(dist_miRNA) <- 0

  cat("\nDistance matrix dimensions:\n")
  cat("mRNA distance:", dim(dist_mRNA), "\n")
  cat("miRNA distance:", dim(dist_miRNA), "\n")

  ####################### 3. Affinity matrix creation #######################

  # Create affinity matrices
  W_mRNA <- affinityMatrix(dist_mRNA, K, sigma)
  W_miRNA <- affinityMatrix(dist_miRNA, K, sigma)

  ########################### 4. Network fusion ############################

  # Perform SNF
  cat("\nPerforming Similarity Network Fusion...\n")

  W_fused <- SNF(list(W_mRNA, W_miRNA), K = K, t = iterations)

  ################ 5. Optimal cluster number determination #################

  if (n_clusters == 0) {
      # Method 1: Eigen-gap based estimation
      estimated_res <- estimateNumberOfClustersGivenGraph(W_fused, NUMC = 2:10)
      optimal_k <- estimated_res$`Eigen-gap 2nd best`

      # Force optimal_k in [3, 7]
      optimal_k <- min(optimal_k, 7)
      optimal_k <- max(optimal_k, 3)
  } else {
    optimal_k <- n_clusters
  }

  # Method 2: Silhouette analysis
  silhouette_scores <- numeric(19)

  for (k in 2:20) {
    clusters_temp <- spectralClustering(W_fused, k)
    sil <- silhouette(clusters_temp, dist = 1 - W_fused)
    silhouette_scores[k - 1] <- mean(sil[, 3])
  }

  # Plot silhouette scores
  plot_silhouette_path <- paste0(dir_plots, "silhouette_analysis.png")

  png(plot_silhouette_path, width = 800, height = 600)
  plot(2:20,
       silhouette_scores,
       type = "b",
       pch = 19,
       xlab = "Number of clusters",
       ylab = "Average silhouette width",
       main = "Silhouette Analysis for Optimal Cluster Number")
  abline(v = optimal_k, col = "red", lty = 2)
  dev.off()

  # Use the optimal number of clusters
  final_clusters <- spectralClustering(W_fused, optimal_k)

  patient_to_cluster <- data.frame(
    cluster = final_clusters
  )

  rownames(patient_to_cluster) <- rownames(W_fused)

  plot_network_fused(W_fused,
                     final_clusters,
                     optimal_k,
                     dir_plots)

  ########################### 6. Enhanced Visualizations ############################

  # Define blue color gradient for heatmaps
  blue_gradient <- colorRampPalette(c("white", "lightblue", "dodgerblue", "darkblue"))(100)

  # 1. Basic cluster visualization with blue gradient
  png(paste0(dir_plots, "fused_network_clusters.png"), width = 800, height = 800)

  W_display <- W_fused
  diag(W_display) <- 0  # Remove diagonal for better visualization

  # Order patients by cluster
  cluster_order <- order(final_clusters)
  W_ordered <- W_display[cluster_order, cluster_order]

  # Plot with blue gradient
  image(1:nrow(W_ordered),
        1:ncol(W_ordered),
        as.matrix(W_ordered),
        col = blue_gradient,
        xlab = "Patients",
        ylab = "Patients",
        main = paste("SNF Fused Network - HCC Patient Similarity (K =", optimal_k, "clusters)"),
        axes = FALSE)

  dev.off()

  # 2. Individual network visualizations with blue gradient
  png(paste0(dir_plots, "individual_networks_enhanced.png"), width = 1600, height = 800)
  par(mfrow = c(1, 2))

  # Function to plot network with blue gradient
  plot_network_blue <- function(W, clusters, title_text) {
    W_display <- W
    diag(W_display) <- 0  # Remove diagonal

    # Order by clusters
    cluster_order <- order(clusters)
    W_ordered <- W_display[cluster_order, cluster_order]

    # Plot
    image(1:nrow(W_ordered),
          1:ncol(W_ordered),
          as.matrix(W_ordered),
          col = blue_gradient,
          xlab = "", ylab = "",
          main = title_text,
          axes = FALSE)
  }

  # Plot each network
  plot_network_blue(W_mRNA, final_clusters, "mRNA Expression Network")
  plot_network_blue(W_miRNA, final_clusters, "miRNA Expression Network")

  dev.off()

  # 3. Heatmap visualization with cluster annotation and ordering
  # Order patients by cluster assignment for better visualization
  cluster_order <- order(final_clusters)
  W_fused_ordered <- W_fused[cluster_order, cluster_order]

  # Remove diag values to increase heatmap readability
  diag(W_fused_ordered) <- 0

  # Create annotation dataframe with ordered clusters
  cluster_annotation <- data.frame(
    Cluster = as.factor(final_clusters[cluster_order]),
    row.names = rownames(W_fused)[cluster_order]
  )

  # Define colors for annotation
  ann_colors <- list(
    Cluster = setNames(brewer.pal(max(3, optimal_k), "Set1")[1:optimal_k], 1:optimal_k)
  )

  # Create the heatmap with ordered data
  png(paste0(dir_plots, "similarity_heatmap_annotated_ordered.png"), width = 1400, height = 1400, res = 300)
  print(pheatmap(W_fused_ordered,
           cluster_rows = FALSE,  # Don't cluster since already ordered
           cluster_cols = FALSE,  # Don't cluster since already ordered
           annotation_row = cluster_annotation,
           annotation_colors = ann_colors,
           show_rownames = FALSE,
           show_colnames = FALSE,
           color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
           main = "Patient Similarity Matrix Ordered by Clusters",
           name = "Similarity",
           fontsize = 8,
           fontsize_row = 6,
           fontsize_col = 6))
  dev.off()

  png(paste0(dir_plots, "similarity_heatmap_with_boundaries.png"), width = 1400, height = 1400, res = 300)

  # Calculate cluster boundaries
  cluster_sizes <- table(final_clusters[cluster_order])
  cluster_boundaries <- cumsum(cluster_sizes)

  # Create the heatmap
  print(pheatmap(W_fused_ordered,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                annotation_row = cluster_annotation,
                annotation_colors = ann_colors,
                show_rownames = FALSE,
                show_colnames = FALSE,
                color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
                main = "Patient Similarity Matrix with Cluster Boundaries",
                name = "Similarity",
                fontsize = 8,
                fontsize_row = 6,
                fontsize_col = 6,
                gaps_row = cluster_boundaries[-length(cluster_boundaries)],
                gaps_col = cluster_boundaries[-length(cluster_boundaries)],
                border_color = NA,
                annotation_legend = TRUE,
                legend = TRUE))

  dev.off()

  # 5. Consensus clustering for stability assessment
  png(paste0(dir_plots, "consensus_clustering.png"), width = 1200, height = 800)

  par(mfrow = c(2, 3))
  for (k in 2:7) {
    clusters_k <- spectralClustering(W_fused, k)
    displayClusters(W_fused, clusters_k)
    title(main = paste("K =", k, "clusters"))
  }

  dev.off()

  ########################### 7. Cluster Quality Metrics ############################

  # Calculate silhouette coefficient
  sil <- silhouette(final_clusters, dist = 1 - W_fused)
  avg_sil <- mean(sil[, 3])

  cat("\n=== Clustering Results ===\n")
  cat("Optimal number of clusters:", optimal_k, "\n")
  cat("Average silhouette width:", round(avg_sil, 3), "\n")

  # Silhouette plot
  png(paste0(dir_plots, "silhouette_plot.png"), width = 800, height = 600)
  plot(sil,
       col = brewer.pal(max(3, optimal_k), "Set1")[1:optimal_k],
       main = paste("Silhouette Plot (Average width:", round(avg_sil, 3), ")"))
  dev.off()

  # Calculate concordance between individual networks
  concordance_matrix <- concordanceNetworkNMI(list(W_mRNA, W_miRNA), optimal_k)
  cat("\nNetwork Concordance (NMI):\n")
  print(round(concordance_matrix, 3))

  return(list(
    W_fused = W_fused,
    dist_mRNA = dist_mRNA,
    dist_miRNA = dist_miRNA,
    final_clusters = final_clusters,
    patient_to_cluster = patient_to_cluster
  ))
}