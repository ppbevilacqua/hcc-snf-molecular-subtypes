calculate_modularity <- function(W, clusters) {
  # Simple modularity calculation
  m <- sum(W) / 2
  Q <- 0
  for (c in unique(clusters)) {
    nodes_c <- which(clusters == c)
    if (length(nodes_c) > 1) {
      W_c <- W[nodes_c, nodes_c]
      e_c <- sum(W_c) / (2 * m)
      a_c <- sum(W[nodes_c,]) / (2 * m)
      Q <- Q + e_c - a_c^2
    }
  }
  return(Q)
}

grid_search_SNF <- function(dir_plots,
                            dist_mRNA,
                            dist_miRNA,
                            K_range,
                            alpha_range,
                            T_range) {

  # Create results storage
  param_results <- data.frame()

  cat("\nStarting parameter grid search...\n")
  cat("Total combinations to test:", length(K_range) *
    length(alpha_range) *
    length(T_range), "\n\n")

  # Grid search
  iteration <- 0
  total_iterations <- length(K_range) *
    length(alpha_range) *
    length(T_range) *
    4 # number of tested clusters

  for (K in K_range) {
    for (alpha in alpha_range) {
      for (T_iter in T_range) {

        cat(sprintf("\nTesting K=%d, alpha=%.1f, T=%d\n", K, alpha, T_iter))

        # Create affinity matrices with current parameters
        W_mRNA_temp <- affinityMatrix(dist_mRNA, K = K, sigma = alpha)
        W_miRNA_temp <- affinityMatrix(dist_miRNA, K = K, sigma = alpha)

        # Perform SNF
        W_fused_temp <- SNF(list(W_mRNA_temp, W_miRNA_temp),
                            K = K, t = T_iter)

        # Test different cluster numbers (3:6)
        for (n_clust in 3:6) {

          iteration <- iteration + 1
          cat(sprintf("  Testing %d clusters (Progress: %d/%d)\r",
                      n_clust, iteration, total_iterations))

          # Get clusters
          clusters <- spectralClustering(W_fused_temp, n_clust)

          # Skip if any cluster is too small
          if (min(table(clusters)) < 3) {
            next
          }

          # Calculate metrics

          # 1. Silhouette score
          sil <- silhouette(clusters, dist = 1 - W_fused_temp)
          avg_silhouette <- mean(sil[, 3])

          # 2. Within-cluster similarity (cohesion)
          within_sim <- 0
          for (c in unique(clusters)) {
            idx <- which(clusters == c)
            if (length(idx) > 1) {
              within_sim <- within_sim + mean(W_fused_temp[idx, idx][upper.tri(W_fused_temp[idx, idx])])
            }
          }
          avg_within_sim <- within_sim / length(unique(clusters))

          # 3. Between-cluster similarity (separation)
          between_sim <- 0
          count <- 0
          for (c1 in unique(clusters)) {
            for (c2 in unique(clusters)) {
              if (c1 < c2) {
                idx1 <- which(clusters == c1)
                idx2 <- which(clusters == c2)
                between_sim <- between_sim + mean(W_fused_temp[idx1, idx2])
                count <- count + 1
              }
            }
          }
          avg_between_sim <- if (count > 0) between_sim / count else 0

          # 4. Separation ratio
          separation_ratio <- if (avg_between_sim > 0) avg_within_sim / avg_between_sim else 0

          # 5. Modularity
          modularity <- calculate_modularity(W_fused_temp, clusters)

          # 6. Concordance with individual networks
          concordance_matrix <- concordanceNetworkNMI(
            list(W_mRNA_temp, W_miRNA_temp),
            n_clust
          )
          avg_concordance <- mean(concordance_matrix[upper.tri(concordance_matrix)])

          # 7. Cluster size statistics
          cluster_sizes <- table(clusters)
          size_imbalance <- sd(cluster_sizes) / mean(cluster_sizes)  # Coefficient of variation

          # Store results
          param_results <- rbind(param_results, data.frame(
            K = K,
            alpha = alpha,
            T = T_iter,
            n_clusters = n_clust,
            silhouette = avg_silhouette,
            within_similarity = avg_within_sim,
            between_similarity = avg_between_sim,
            separation_ratio = separation_ratio,
            modularity = modularity,
            concordance = avg_concordance,
            min_cluster_size = min(cluster_sizes),
            max_cluster_size = max(cluster_sizes),
            size_imbalance = size_imbalance
          ))
        }
      }
    }
  }

  cat("\n\nGrid search completed!\n")

  ########################### Analyze Results ###########################

  # Remove any NA values
  param_results <- na.omit(param_results)

  # Create composite score focusing on key metrics
  param_results$composite_score <-
    0.4 * pmax(param_results$silhouette, 0) +  # Use max(0, silhouette) to handle negative values
      0.3 * (param_results$separation_ratio / max(param_results$separation_ratio)) +
      0.2 * param_results$concordance +
      0.1 * (1 - param_results$size_imbalance / max(param_results$size_imbalance))

  # Get top 10 parameter combinations
  top_params <- param_results[order(param_results$composite_score, decreasing = TRUE),][1:10,]

  cat("\nTop 10 parameter combinations:\n")
  print(top_params[, c("K", "alpha", "T", "n_clusters", "silhouette",
                       "separation_ratio", "composite_score")])

  # Find best parameters by different criteria
  best_by_silhouette <- param_results[which.max(param_results$silhouette),]
  best_by_separation <- param_results[which.max(param_results$separation_ratio),]
  best_by_composite <- param_results[which.max(param_results$composite_score),]

  cat("\n\nBest parameters by different criteria:")
  cat("\n1. By Silhouette Score:")
  cat(sprintf("\n   K=%d, alpha=%.1f, T=%d, clusters=%d (silhouette=%.3f)",
              best_by_silhouette$K, best_by_silhouette$alpha,
              best_by_silhouette$T, best_by_silhouette$n_clusters,
              best_by_silhouette$silhouette))

  cat("\n\n2. By Separation Ratio:")
  cat(sprintf("\n   K=%d, alpha=%.1f, T=%d, clusters=%d (ratio=%.3f)",
              best_by_separation$K, best_by_separation$alpha,
              best_by_separation$T, best_by_separation$n_clusters,
              best_by_separation$separation_ratio))

  cat("\n\n3. By Composite Score:")
  cat(sprintf("\n   K=%d, alpha=%.1f, T=%d, clusters=%d (score=%.3f)",
              best_by_composite$K, best_by_composite$alpha,
              best_by_composite$T, best_by_composite$n_clusters,
              best_by_composite$composite_score))

  ########################### Visualize Results ###########################

  # Plot 1: Silhouette score heatmap
  p1 <- ggplot(param_results, aes(x = factor(K), y = factor(alpha), fill = silhouette)) +
    geom_tile() +
    facet_grid(n_clusters ~ T, labeller = label_both) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
    theme_minimal() +
    labs(title = "Silhouette Score Heatmap",
         x = "K (neighbors)", y = "Alpha", fill = "Silhouette")

  # Plot 2: Separation ratio vs silhouette
  p2 <- ggplot(param_results, aes(x = silhouette, y = separation_ratio,
                                  color = factor(n_clusters), size = concordance)) +
    geom_point(alpha = 0.6) +
    theme_minimal() +
    labs(title = "Clustering Quality Metrics",
         x = "Silhouette Score", y = "Separation Ratio",
         color = "Clusters", size = "Concordance")

  # Plot 3: Parameter importance
  param_summary <- aggregate(composite_score ~ K + alpha + T, data = param_results, mean)
  p3 <- ggplot(param_summary, aes(x = factor(K), y = composite_score,
                                  group = factor(alpha), color = factor(alpha))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    facet_wrap(~T, labeller = label_both) +
    theme_minimal() +
    labs(title = "Average Composite Score by Parameters",
         x = "K (neighbors)", y = "Composite Score", color = "Alpha")

  # Save plots
  png(paste0(dir_plots, "snf_parameter_optimization.png"), width = 1400, height = 1000)
  grid.arrange(p1, p2, p3, ncol = 2)
  dev.off()

  ########################### Apply Best Parameters ###########################

  # Use best parameters based on composite score
  best_params <- best_by_composite

  cat("\n\nApplying best parameters:")
  cat(sprintf("\nK = %d, alpha = %.1f, T = %d, n_clusters = %d",
              best_params$K, best_params$alpha, best_params$T, best_params$n_clusters))

  # Apply best parameters for final clustering
  W_mRNA_final <- affinityMatrix(dist_mRNA, K = best_params$K, sigma = best_params$alpha)
  W_miRNA_final <- affinityMatrix(dist_miRNA, K = best_params$K, sigma = best_params$alpha)

  W_fused_final <- SNF(list(W_mRNA_final, W_miRNA_final),
                       K = best_params$K, t = best_params$T)

  final_clusters <- spectralClustering(W_fused_final, best_params$n_clusters)

  cat("\n\nFinal clustering results:")
  cat("\nCluster sizes:", table(final_clusters))
  cat("\nFinal silhouette score:", best_params$silhouette)
  cat("\n")

  # Save results
  write.csv(param_results, paste0(dirPlots, "snf_parameter_search_results.csv"), row.names = FALSE)
  write.csv(top_params, paste0(dirPlots, "snf_top_parameters.csv"), row.names = FALSE)
}