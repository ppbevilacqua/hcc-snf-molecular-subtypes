plot_network_fused <- function(W_fused,
                               final_clusters,
                               optimal_k,
                               dir_plots) {

  # Enhanced Network graph visualization with improved clustering and edge scaling
  png(paste0(dir_plots, "network_graph_visualization_enhanced.png"), width = 1400, height = 1400, res = 300)

  # Create adjacency matrix with smarter thresholding
  W_display <- W_fused
  diag(W_display) <- 0  # Remove self-loops

  # Create igraph object
  g <- graph_from_adjacency_matrix(W_display, mode = "undirected", weighted = TRUE)

  # Keep only top 10% edges by weight
  edge_weights_all <- E(g)$weight
  weight_threshold <- quantile(edge_weights_all, 0.90)
  edges_to_remove <- which(edge_weights_all < weight_threshold)
  g <- delete_edges(g, edges_to_remove)

  clusters_for_graph <- final_clusters

  # Enhanced color palette
  if (optimal_k <= 8) {
    cluster_colors <- brewer.pal(max(3, optimal_k), "Set1")
  } else {
    # For more clusters, use a combination of palettes
    cluster_colors <- rainbow(optimal_k, s = 0.7, v = 0.8)
  }

  # Set node properties
  V(g)$color <- cluster_colors[clusters_for_graph]
  V(g)$cluster <- clusters_for_graph
  V(g)$size <- 7
  # White frame for all clusters, dark gray for yellow (cluster 6)
  frame_colors <- ifelse(clusters_for_graph == 6, "#404040", "white")
  V(g)$frame.color <- frame_colors
  V(g)$frame.width <- 0.4

  # Enhanced edge properties with better scaling
  edge_weights <- E(g)$weight
  min_weight <- min(edge_weights)
  max_weight <- max(edge_weights)

  # Normalize weights to a reasonable range for visualization
  normalized_weights <- (edge_weights - min_weight) / (max_weight - min_weight)
  E(g)$width <- 0.1 + normalized_weights * 0.8  # Width from 0.1 to 0.9 (thinner edges)

  # Color edges based on strength 
  edge_colors <- colorRampPalette(c("#E0E0E0", "#404040"))(100)
  edge_color_indices <- pmax(1, pmin(100, round(normalized_weights * 99 + 1)))
  E(g)$color <- edge_colors[edge_color_indices]

  # Calculate optimal cluster centers in a circle pattern
  cluster_separation_factor <- 12.0
  angles <- seq(0, 2 * pi, length.out = optimal_k + 1)[1:optimal_k]
  radius <- cluster_separation_factor
  cluster_centers <- cbind(radius * cos(angles), radius * sin(angles))

  # Initialize layout matrix
  n_nodes <- vcount(g)
  final_layout <- matrix(0, nrow = n_nodes, ncol = 2)

  # Assign nodes to cluster centers and add some spread within clusters
  for (clust in 1:optimal_k) {
    cluster_nodes <- which(clusters_for_graph == clust)
    n_cluster_nodes <- length(cluster_nodes)

    if (n_cluster_nodes > 0) {
      center <- cluster_centers[clust,]

      if (n_cluster_nodes == 1) {
        # Single node goes exactly to cluster center
        final_layout[cluster_nodes,] <- center
      } else {
        # Multiple nodes: randomly distribute within cluster area
        spread_radius <- min(3.5, 2.5 + 0.15 * sqrt(n_cluster_nodes))

        for (i in 1:n_cluster_nodes) {
          node_idx <- cluster_nodes[i]
          # Random position within a disk (not on perimeter)
          r <- spread_radius * sqrt(runif(1))  # sqrt for uniform distribution in disk
          theta <- runif(1, 0, 2 * pi)

          final_layout[node_idx, 1] <- center[1] + r * cos(theta)
          final_layout[node_idx, 2] <- center[2] + r * sin(theta)
        }
      }
    }
  }

  # Fine-tune layout with small adjustments to improve edge routing
  # Apply small random perturbations to avoid overlapping nodes while 
  # maintaining cluster structure
  for (clust in 1:optimal_k) {
    cluster_nodes <- which(clusters_for_graph == clust)
    if (length(cluster_nodes) > 1) {

      # Apply minimal repulsion within cluster
      for (i in 1:length(cluster_nodes)) {
        for (j in 1:length(cluster_nodes)) {
          if (i != j) {
            node1 <- cluster_nodes[i]
            node2 <- cluster_nodes[j]

            # Calculate distance
            dx <- final_layout[node1, 1] - final_layout[node2, 1]
            dy <- final_layout[node1, 2] - final_layout[node2, 2]
            dist <- sqrt(dx^2 + dy^2)

            # If nodes are too close, push them apart slightly
            if (dist < 0.35 && dist > 0) {
              push_strength <- 0.15
              final_layout[node1, 1] <- final_layout[node1, 1] + push_strength * dx / dist
              final_layout[node1, 2] <- final_layout[node1, 2] + push_strength * dy / dist
            }
          }
        }
      }
    }
  }

  # Scale layout to use more area while leaving space for legend
  layout_range_x <- range(final_layout[, 1])
  layout_range_y <- range(final_layout[, 2])

  # Scale to use larger area 
  target_range <- 1.0
  scale_factor <- target_range / max(diff(layout_range_x), diff(layout_range_y))
  final_layout <- final_layout * scale_factor

  # Center the layout
  center_x <- mean(range(final_layout[, 1]))
  center_y <- mean(range(final_layout[, 2]))
  final_layout[, 1] <- final_layout[, 1] - center_x
  final_layout[, 2] <- final_layout[, 2] - center_y

  # Create the plot with optimized margins
  par(mar = c(1, 1, 2.5, 6))  # Reduced top and bottom margins, adjusted right margin

  # Create list of node indices for each cluster
  cluster_groups <- lapply(1:optimal_k, function(k) which(clusters_for_graph == k))

  # Normalize layout to [-1, 1] range to match plot coordinates
  layout_normalized <- norm_coords(final_layout, xmin = -1, xmax = 1, ymin = -1, ymax = 1)

  # Plot network WITHOUT mark.groups first (only edges and nodes)
  plot(g,
       layout = layout_normalized,
       rescale = FALSE,  # Prevent igraph from rescaling
       vertex.label = NA,
       vertex.label.cex = 0.4,
       vertex.label.color = "black",
       #edge.curved = 0.1,
       xlim = c(-1.1, 1.25),
       ylim = c(-1.1, 1.1),
       main = "Patient Similarity Network (SNF)",
       asp = 1)

  # Draw cluster clouds (fill + border) ON TOP of edges
  for (clust in 1:optimal_k) {
    cluster_nodes <- cluster_groups[[clust]]
    if (length(cluster_nodes) >= 3) {
      cluster_coords <- layout_normalized[cluster_nodes, ]
      hull_indices <- chull(cluster_coords)
      hull_coords <- cluster_coords[hull_indices, ]
      # Calculate centroid from ALL cluster nodes
      centroid_x <- mean(cluster_coords[, 1])
      centroid_y <- mean(cluster_coords[, 2])
      # Expand hull outward from centroid + add padding
      padding <- 0.05  # Absolute padding for node size
      expanded_x <- centroid_x + 1.15 * (hull_coords[, 1] - centroid_x)
      expanded_y <- centroid_y + 1.15 * (hull_coords[, 2] - centroid_y)
      # Add outward padding based on direction from centroid
      for (i in 1:length(expanded_x)) {
        dir_x <- expanded_x[i] - centroid_x
        dir_y <- expanded_y[i] - centroid_y
        dir_len <- sqrt(dir_x^2 + dir_y^2)
        if (dir_len > 0) {
          expanded_x[i] <- expanded_x[i] + padding * (dir_x / dir_len)
          expanded_y[i] <- expanded_y[i] + padding * (dir_y / dir_len)
        }
      }
      # Draw filled polygon with border
      polygon(expanded_x, expanded_y,
              col = adjustcolor(cluster_colors[clust], alpha.f = 0.2),
              border = cluster_colors[clust],
              lwd = 1.5)
    } else if (length(cluster_nodes) == 2) {
      # For 2 nodes, draw ellipse around them
      coords <- layout_normalized[cluster_nodes, ]
      center <- colMeans(coords)
      rad <- max(0.08, dist(coords) / 2 + 0.05)
      symbols(center[1], center[2], circles = rad, add = TRUE,
              bg = adjustcolor(cluster_colors[clust], alpha.f = 0.2),
              fg = cluster_colors[clust], lwd = 1.5, inches = FALSE)
    } else if (length(cluster_nodes) == 1) {
      # For single node, draw circle around it
      coords <- layout_normalized[cluster_nodes, ]
      symbols(coords[1], coords[2], circles = 0.05, add = TRUE,
              bg = adjustcolor(cluster_colors[clust], alpha.f = 0.2),
              fg = cluster_colors[clust], lwd = 1.5, inches = FALSE)
    }
  }

  # Replot nodes on top using same normalized layout
  plot(g,
       layout = layout_normalized,
       rescale = FALSE,
       vertex.label = NA,
       edge.color = NA,  # Hide edges (already drawn)
       edge.width = 0,
       xlim = c(-1.1, 1.25),
       ylim = c(-1.1, 1.1),
       add = TRUE,       # Add to existing plot
       asp = 1)

  # Cluster legend
  legend(x = 1.12, y = 1.1,
         legend = paste("Cluster", 1:optimal_k),
         fill = cluster_colors[1:optimal_k],
         border = "black",
         bty = "n",
         cex = 0.9,
         xpd = TRUE)

  dev.off()


  #######################################################################################

  # Additional visualization: Edge weight distribution
  png(paste0(dir_plots, "edge_weight_distribution.png"), width = 2000, height = 1800, res = 300)
  par(mfrow = c(1, 2))

  # Histogram of edge weights
  hist(E(g)$weight,
       main = "Distribution of Edge Weights\nin Network",
       xlab = "Similarity Values",
       ylab = "Frequency",
       col = "lightblue",
       border = "darkblue",
       breaks = 20,
       cex.main = 0.9,  # Reduced title size
       cex.lab = 0.8,   # Reduced axis labels
       cex.axis = 0.7)  # Reduced axis text

  # Box plot of edge weights by connection type
  # Classify edges as within-cluster or between-cluster
  edge_list <- get.edgelist(g, names = FALSE)  # Get numeric indices instead of names
  edge_types <- character(nrow(edge_list))
  for (i in 1:nrow(edge_list)) {
    node1 <- edge_list[i, 1]
    node2 <- edge_list[i, 2]
    if (!is.na(node1) &&
      !is.na(node2) &&
      node1 <= length(clusters_for_graph) &&
      node2 <= length(clusters_for_graph)) {
      if (clusters_for_graph[node1] == clusters_for_graph[node2]) {
        edge_types[i] <- "Within Cluster"
      } else {
        edge_types[i] <- "Between Clusters"
      }
    } else {
      edge_types[i] <- "Unknown"
    }
  }

  # Remove "Unknown" edges if any exist
  valid_edges <- edge_types != "Unknown"
  if (sum(valid_edges) > 0) {
    boxplot(E(g)$weight[valid_edges] ~ edge_types[valid_edges],
            main = "Edge Weights by Connection Type",
            xlab = "",
            ylab = "Similarity Values",
            col = c("lightgreen", "lightcoral"),
            notch = TRUE,
            cex.main = 0.9,
            cex.lab = 0.8,
            cex.axis = 0.7)
  } else {
    plot.new()
    text(0.5, 0.5, "No valid edges for comparison", cex = 1.2, adj = 0.5)
  }

  dev.off()
}