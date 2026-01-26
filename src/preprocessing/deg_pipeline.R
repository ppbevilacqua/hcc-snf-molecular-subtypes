apply_differential_expression_analysis <- function(
  data_matrix_path,       # input data matrix path
  prc_IQR,                # percentile IQR
  thr_fc,                 # threshold full change (2 standard for TCGA)
  thr_pval,               # threshold p-value
  paired,                 # statistical paired test
  generate_plot = FALSE   # boolean to set plot generation
) {

  ################### CONSTANTS ###################

  if (generate_plot) {
    filename_heatmap <- paste0(dirPlotDEG, "heatmap.pdf")
  }

  filename_list_normal <- paste0(dirResults, "normal.txt")
  filename_list_cancer <- paste0(dirResults, "tumor.txt")
  filename_DEG <- paste0(dirResults, "DEG.txt")
  #################################################


  ################# Pre-processing #################

  # STEP 1. Importing data
  tmp <- read.table(data_matrix_path, header = T, sep = "\t", check.names = F,
                    row.names = 1, quote = "", nrows = 10)

  classes <- sapply(tmp, class)

  tmp <- read.table(data_matrix_path, header = T, sep = "\t", check.names = F,
                    row.names = 1, quote = "", colClasses = classes)

  genes <- rownames(tmp) # 60660
  pz <- colnames(tmp)    # 424

  pzN <- grep(NORMAL_PZ_REGEX, pz, value = TRUE) # 50
  pzC <- grep(CANCER_PZ_REGEX, pz, value = TRUE) # 374

  pzN <- pzN[!duplicated(str_extract(pzN, PZ_TCGA_REGEX))] # 50
  pzC <- pzC[!duplicated(str_extract(pzC, PZ_TCGA_REGEX))] # 371

  nameN <- str_extract(pzN, PZ_TCGA_REGEX)
  nameC <- str_extract(pzC, PZ_TCGA_REGEX)

  common_pz <- intersect(nameN, nameC) # 50

  pzN_com <- unlist(lapply(common_pz, function(x) { grep(x, pzN, value = TRUE) }))
  pzC_com <- unlist(lapply(common_pz, function(x) { grep(x, pzC, value = TRUE) }))

  dataN <- tmp[, pzN_com] # 60660   50
  dataC <- tmp[, pzC_com] # 60660   50

  data <- cbind(dataN, dataC)
  pz_com <- colnames(data)

  rm(tmp, pz, nameN, nameC, pzN, pzC, common_pz, classes, pz_com)

  ##################################################

  # STEP 2. Analysis

  # STEP 2.0: Overall mean
  overall_mean <- rowMeans(data)
  ind <- which(overall_mean == 0) # 247

  dataN <- dataN[-ind,]           # 799  49
  dataC <- dataC[-ind,]
  data <- data[-ind,]
  genes <- genes[-ind]

  rm(ind)

  # STEP 2.1: Logarithmic transformation (to reduce data variability)
  dataN <- log2(dataN + 1)
  dataC <- log2(dataC + 1)
  data <- log2(data + 1)

  # STEP 2.2 Pre-processing IQR

  # IQR helps us to see the variation of data around the means
  variation <- apply(data, 1, IQR)

  # IQR filtering
  thr_prc <- quantile(variation, prc_IQR) # 31% 0.04011265
  ind <- which(variation <= thr_prc)

  dataN <- dataN[-ind,]           # 551  49
  dataC <- dataC[-ind,]
  data <- data[-ind,]
  genes <- genes[-ind]

  rm(ind)

  if (generate_plot) {
    png(file = paste0(dirPlotDEG, "IQR_histo.png"), width = PNG_WIDTH, height = PNG_HEIGHT)

    # IQR distribution
    hist(variation,
         main = "IQR frequency distribution",
         breaks = 100,
         xlab = "IQR value",
         ylab = "Frequency",
         col = "gold")

    # Plot threshold
    abline(v = thr_prc, lty = 2, lwd = 4, col = "grey")
    dev.off()
  }

  ##################################################

  # STEP 2.3 Filtering
  logFC <- rowMeans(dataC) - rowMeans(dataN)

  if (generate_plot) {

    png(file = paste0(dirPlotDEG, "FC_log_histo.png"), width = PNG_WIDTH, height = PNG_HEIGHT)

    hist(logFC,
         main = "FC (logarithmic) frequency distribution",
         breaks = 100,
         xlab = "log FC",
         ylab = "Frequency",
         col = "gold")

    abline(v = c(-log2(thr_fc), log2(thr_fc)), lty = 2, lwd = 4, col = 'gray')
    dev.off()

  }

  # logFC filtering - discard what's in the middle of the hist
  ind <- which(abs(logFC) < log2(thr_fc)) # 456

  # check if the threshold allows us to get values
  if (length(ind) > 0) {
    dataN <- dataN[-ind,]           # 95
    dataC <- dataC[-ind,]
    data <- data[-ind,]
    genes <- genes[-ind]
    logFC <- logFC[-ind]
  }

  rm(ind)

  # p-value computation - parametric test with normal distribution -> t-test
  N <- ncol(dataN)    # 50
  M <- ncol(dataC)    # 50

  pval <- apply(data, 1, function(x) {
    t.test(x[1:N], x[(N + 1):(N + M)], paired = paired)$p.value
  })

  # adjustment p-value
  pval_adj <- p.adjust(pval, method = "fdr")

  # p-value filtering
  ind <- which(pval_adj > thr_pval)  # 0

  if (length(ind) > 0) {
    dataN <- dataN[-ind,]
    dataC <- dataC[-ind,]
    data <- data[-ind,]
    genes <- genes[-ind]
    logFC <- logFC[-ind]
    pval <- pval[-ind]
    pval_adj <- pval_adj[-ind]
  }

  rm(ind)

  # STEP 3. Exporting results

  direction <- ifelse(logFC > 0, "UP", "DOWN") # 24 71

  DEG <- data.frame(str_split_fixed(genes, "\\|", 2))
  colnames(DEG) <- c("geneSymbol", "ensembl_id")

  # list of DEG that are results of the analysis for test in BRCA
  results <- data.frame(genes = DEG$geneSymbol,
                        ensable_id = DEG$ensembl_id,
                        pvalue = pval,
                        pval_adj = pval_adj,
                        logFC = logFC,
                        direction = direction) # 1275

  # order results from up regulated to down ones
  results <- results[order(results$logFC, decreasing = T),]

  write.table(results, file = filename_DEG, row.names = F,
              sep = "\t", quote = F)

  write.table(data, file = matrix_DEG_path, row.names = T,
              col.names = NA, sep = "\t", quote = F)

  write.table(pzN_com, filename_list_normal,
              sep = "\t", col.names = F, row.names = F, quote = F)

  write.table(pzC_com, filename_list_cancer,
              sep = "\t", col.names = F, row.names = F, quote = F)

  # STEP 4. Plot
  if (generate_plot) {

    # volcano plot
    png(file = paste0(dirPlotDEG, "volcano.png"), width = PNG_WIDTH, height = PNG_HEIGHT)

    plot(logFC, -log10(pval_adj),
         main = "Volcano plot",
         xlim = c(-8, 8),
         ylim = c(0, 25),
         xlab = "log2 fold change",
         ylab = "-log10 p-value",
         pch = 20, col = "deepskyblue4", cex = 0.8)

    abline(h = -log10(thr_pval), lty = 2, lwd = 2, col = "brown2")
    abline(v = c(-log2(thr_fc), log2(thr_fc)), lty = 2, lwd = 2, col = "grey")
    dev.off()

    # box plot
    ind <- which.max(logFC) # most up-regulated gene or min
    gene_id <- genes[ind]

    df = data.frame(normal = t(dataN[ind,]),
                    cancer = t(dataC[ind,]),
                    row.names = NULL)

    colnames(df) <- c(paste(gene_id, "normal"),
                      paste(gene_id, "cancer"))

    png(file = paste0(dirPlotDEG, "boxplot.png"), width = PNG_WIDTH, height = PNG_HEIGHT)

    boxplot(df,
            main = paste0(gene_id, " adjusted pvalue = ",
                          format(pval_adj[ind], digits = 2)),
            notch = T,
            ylab = "Gene expression value",
            xlab = "Condition",
            col = c("gold", "deepskyblue4"),
            pars = list(boxwex = 0.3, staplewex = 0.6))

    dev.off()

    # pie chart
    count <- table(results$direction)

    png(file = paste0(dirPlotDEG, "pie_chart.png"), width = PNG_WIDTH, height = PNG_HEIGHT)

    pie(count,
        main = "FC UP and DOWN regulating gene count",
        labels = paste0(names(count), " ", round(100 * count / sum(count), 2), "%"),
        col = c("brown2", "darkseagreen"))
    dev.off()

    # heatmap
    test <- grepl(NORMAL_PZ_REGEX, colnames(data))

    condition <- ifelse(test, "normal", "cancer")

    annotation <- data.frame(condition = condition)
    rownames(annotation) <- colnames(data)

    vector_color <- c("green", "orange")
    names(vector_color) <- unique(condition)

    annotation_colors <- list(condition = vector_color)

    # heatmap
    pheatmap(data, scale = "row",
             border_color = NA,
             cluster_cols = T,
             cluster_rows = T,
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             clustering_method = "complete",
             annotation_col = annotation,
             annotation_colors = annotation_colors,
             color = colorRampPalette(colors = c("blue",
                                                 "blue3", "black", "yellow3", "yellow"))(100),
             show_rownames = F,
             show_colnames = F,
             cutree_cols = 2,
             cutree_rows = 2,
             width = 10, height = 10,
             filename = filename_heatmap)
  }

  return(data)
}


apply_DEG_analysis_pvalue_only <- function(
  data_matrix_path,       # input data matrix path
  prc_IQR,                # percentile IQR
  thr_pval,               # threshold p-value
  paired,                 # statistical paired test
  thr_fc = 1,             # threshold fold-change
  generate_plot = FALSE   # boolean to set plot generation
) {

  ################### CONSTANTS ###################

  if (generate_plot) {
    filename_heatmap <- paste0(dirPlotDEG, "heatmap.pdf")
  }

  filename_list_normal <- paste0(dirResults, "normal.txt")
  filename_list_cancer <- paste0(dirResults, "tumor.txt")
  filename_DEG <- paste0(dirResults, "DEG.txt")
  filename_cancer_matrix <- paste0(dirResults, "cancer_matrix_DEG.txt")
  #################################################

  ################# Pre-processing #################

  # STEP 1. Importing data
  tmp <- read.table(data_matrix_path, header = T, sep = "\t", check.names = F,
                    row.names = 1, quote = "", nrows = 10)

  classes <- sapply(tmp, class)

  tmp <- read.table(data_matrix_path, header = T, sep = "\t", check.names = F,
                    row.names = 1, quote = "", colClasses = classes)

  genes <- rownames(tmp)
  pz <- colnames(tmp)

  pzN <- grep(NORMAL_PZ_REGEX, pz, value = TRUE)
  pzC <- grep(CANCER_PZ_REGEX, pz, value = TRUE)

  # Remove duplicates based on patient ID
  pzN <- pzN[!duplicated(str_extract(pzN, PZ_TCGA_REGEX))]
  pzC <- pzC[!duplicated(str_extract(pzC, PZ_TCGA_REGEX))]

  nameN <- str_extract(pzN, PZ_TCGA_REGEX)
  nameC <- str_extract(pzC, PZ_TCGA_REGEX)

  # Find common patients between normal and cancer (for paired analysis)
  common_pz <- intersect(nameN, nameC)

  pzN_com <- unlist(lapply(common_pz, function(x) { grep(x, pzN, value = TRUE) }))

  # Keep ALL cancer patients, not just common ones
  dataN <- tmp[, pzN_com] # Normal samples from common patients
  dataC <- tmp[, pzC]     # ALL cancer samples

  # For statistical testing, we'll use only the paired samples
  dataC_paired <- tmp[, unlist(lapply(common_pz, function(x) { grep(x, pzC, value = TRUE) }))]

  # Create combined data for filtering (using paired samples for statistics)
  data_for_stats <- cbind(dataN, dataC_paired)

  rm(tmp, pz, nameN, nameC, pzN, common_pz, classes)

  ##################################################

  # STEP 2. Analysis

  # STEP 2.0: Overall mean (using all data)
  all_data <- cbind(dataN, dataC)
  overall_mean <- rowMeans(all_data)
  ind <- which(overall_mean == 0) # genes with zero expression

  if (length(ind) > 0) {
    dataN <- dataN[-ind,]
    dataC <- dataC[-ind,]
    dataC_paired <- dataC_paired[-ind,]
    data_for_stats <- data_for_stats[-ind,]
    all_data <- all_data[-ind,]
    genes <- genes[-ind]
  }

  rm(ind)

  # STEP 2.1: Logarithmic transformation (to reduce data variability)
  dataN <- log2(dataN + 1)
  dataC <- log2(dataC + 1)
  dataC_paired <- log2(dataC_paired + 1)
  data_for_stats <- log2(data_for_stats + 1)
  all_data <- log2(all_data + 1)

  # STEP 2.2 Pre-processing IQR (using all data)
  variation <- apply(all_data, 1, IQR)

  # IQR filtering
  thr_prc <- quantile(variation, prc_IQR)
  ind <- which(variation <= thr_prc)

  if (length(ind) > 0) {
    dataN <- dataN[-ind,]
    dataC <- dataC[-ind,]
    dataC_paired <- dataC_paired[-ind,]
    data_for_stats <- data_for_stats[-ind,]
    genes <- genes[-ind]
    variation <- variation[-ind]
  }

  rm(ind)

  if (generate_plot) {
    png(file = paste0(dirPlotDEG, "IQR_histo.png"), width = PNG_WIDTH, height = PNG_HEIGHT)

    # IQR distribution
    hist(variation,
         main = "IQR frequency distribution",
         breaks = 100,
         xlab = "IQR value",
         ylab = "Frequency",
         col = "gold")

    # Plot threshold
    abline(v = thr_prc, lty = 2, lwd = 4, col = "grey")
    dev.off()
  }

  ##################################################

  # STEP 2.3 Statistical testing (using paired samples only)
  # Calculate logFC for plotting purposes
  logFC <- rowMeans(dataC_paired) - rowMeans(dataN)

  if (generate_plot) {
    png(file = paste0(dirPlotDEG, "FC_log_histo.png"), width = PNG_WIDTH, height = PNG_HEIGHT)

    hist(logFC,
         main = "FC (logarithmic) frequency distribution",
         breaks = 100,
         xlab = "log FC",
         ylab = "Frequency",
         col = "gold")

    if (thr_fc > 1) {
      abline(v = c(-log2(thr_fc), log2(thr_fc)), lty = 2, lwd = 4, col = 'gray')
    }
    dev.off()
  }

  # Fold-change filtering
  if (thr_fc > 1) {
    # logFC filtering - discard what's in the middle of the hist
    ind <- which(abs(logFC) < log2(thr_fc))

    # check if the threshold allows us to get values
    if (length(ind) > 0) {
      dataN <- dataN[-ind,]
      dataC <- dataC[-ind,]
      dataC_paired <- dataC_paired[-ind,]
      data_for_stats <- data_for_stats[-ind,]
      genes <- genes[-ind]
      logFC <- logFC[-ind]
    }

    rm(ind)
  }

  # p-value computation - parametric test with normal distribution -> t-test
  N <- ncol(dataN)            # number of normal samples
  M <- ncol(dataC_paired)     # number of paired cancer samples

  pval <- apply(data_for_stats, 1, function(x) {
    t.test(x[1:N], x[(N + 1):(N + M)], paired = paired)$p.value
  })

  # adjustment p-value
  pval_adj <- p.adjust(pval, method = "fdr")

  # ONLY p-value filtering
  ind <- which(pval_adj > thr_pval)

  if (length(ind) > 0) {
    dataN <- dataN[-ind,]
    dataC <- dataC[-ind,]
    dataC_paired <- dataC_paired[-ind,]
    genes <- genes[-ind]
    logFC <- logFC[-ind]
    pval <- pval[-ind]
    pval_adj <- pval_adj[-ind]
  }

  rm(ind)

  # STEP 3. Exporting results

  direction <- ifelse(logFC > 0, "UP", "DOWN")

  DEG <- data.frame(str_split_fixed(genes, "\\|", 2))
  colnames(DEG) <- c("geneSymbol", "ensembl_id")

  # list of DEG that are results of the analysis
  results <- data.frame(genes = DEG$geneSymbol,
                        ensable_id = DEG$ensembl_id,
                        pvalue = pval,
                        pval_adj = pval_adj,
                        logFC = logFC,
                        direction = direction)

  # order results from up regulated to down ones
  results <- results[order(results$logFC, decreasing = T),]

  write.table(results, file = filename_DEG, row.names = F,
              sep = "\t", quote = F)

  # Export matrix with ALL cancer samples and DEGs
  write.table(dataC, file = filename_cancer_matrix, row.names = T,
              col.names = NA, sep = "\t", quote = F)

  write.table(pzN_com, filename_list_normal,
              sep = "\t", col.names = F, row.names = F, quote = F)

  write.table(pzC, filename_list_cancer,  # ALL cancer samples
              sep = "\t", col.names = F, row.names = F, quote = F)

  # STEP 4. Plot
  if (generate_plot) {

    # volcano plot (using paired samples for statistics)
    png(file = paste0(dirPlotDEG, "volcano.png"), width = PNG_WIDTH, height = PNG_HEIGHT)

    plot(logFC, -log10(pval_adj),
         main = "Volcano plot",
         xlim = c(-8, 8),
         ylim = c(0, 25),
         xlab = "log2 fold change",
         ylab = "-log10 p-value",
         pch = 20, col = "deepskyblue4", cex = 0.8)

    abline(h = -log10(thr_pval), lty = 2, lwd = 2, col = "brown2")

    if (thr_fc > 1) {
      abline(v = c(-log2(thr_fc), log2(thr_fc)), lty = 2, lwd = 2, col = "grey")
    }

    dev.off()

    # box plot
    ind <- which.max(logFC) # most differentially expressed gene
    gene_id <- genes[ind]

    # Use paired samples for boxplot
    df <- data.frame(normal = t(dataN[ind,]),
                     cancer = t(dataC_paired[ind,]),
                     row.names = NULL)

    colnames(df) <- c(paste(gene_id, "normal"),
                      paste(gene_id, "cancer"))

    png(file = paste0(dirPlotDEG, "boxplot.png"), width = PNG_WIDTH, height = PNG_HEIGHT)

    boxplot(df,
            main = paste0(gene_id, " adjusted pvalue = ",
                          format(pval_adj[ind], digits = 2)),
            notch = T,
            ylab = "Gene expression value",
            xlab = "Condition",
            col = c("gold", "deepskyblue4"),
            pars = list(boxwex = 0.3, staplewex = 0.6))

    dev.off()

    # pie chart
    count <- table(results$direction)

    png(file = paste0(dirPlotDEG, "pie_chart.png"), width = PNG_WIDTH, height = PNG_HEIGHT)

    pie(count,
        main = "UP and DOWN regulating gene count",
        labels = paste0(names(count), " ", round(100 * count / sum(count), 2), "%"),
        col = c("brown2", "darkseagreen"))
    dev.off()

    # heatmap (using ALL cancer samples)
    # Create annotation for all samples
    all_samples_data <- cbind(dataN, dataC)
    test_normal <- grepl(NORMAL_PZ_REGEX, colnames(all_samples_data))
    condition <- ifelse(test_normal, "normal", "cancer")

    annotation <- data.frame(condition = condition)
    rownames(annotation) <- colnames(all_samples_data)

    vector_color <- c("green", "orange")
    names(vector_color) <- unique(condition)

    annotation_colors <- list(condition = vector_color)

    # heatmap with all samples
    pheatmap(all_samples_data, scale = "row",
             border_color = NA,
             cluster_cols = T,
             cluster_rows = T,
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             clustering_method = "complete",
             annotation_col = annotation,
             annotation_colors = annotation_colors,
             color = colorRampPalette(colors = c("blue",
                                                 "blue3", "black", "yellow3", "yellow"))(100),
             show_rownames = F,
             show_colnames = F,
             cutree_cols = 2,
             cutree_rows = 2,
             width = 10, height = 10,
             filename = filename_heatmap)
  }

  # Return matrix with DEGs and ALL cancer samples
  return(dataC)
}


get_mRNAseq_DEG_dataframe <- function(data_mRNAseq_matrix) {

  if (file.exists(matrix_DEG_mRNA_cancer_path)) {

    data <- read.table(matrix_DEG_mRNA_cancer_path, header = T,
                       sep = "\t", check.names = F,
                       row.names = 1, quote = "", nrows = 1)

    classes <- sapply(data, class)

    data <- read.table(matrix_DEG_mRNA_cancer_path, header = T,
                       sep = "\t", check.names = F,
                       row.names = 1, quote = "", colClasses = classes)

    return(data)

  } else {

    return(apply_DEG_analysis_pvalue_only(
      data_matrix_path = data_mRNAseq_matrix,
      prc_IQR = 0.64,
      thr_pval = 0.01,
      thr_fc = 1,
      paired = T,
      generate_plot = F
    ))
  }
}
