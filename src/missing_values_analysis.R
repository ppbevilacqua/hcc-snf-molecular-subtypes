# Function to analyze and visualize missing values in clinical data
analyze_missing_values <- function(clinical_data, 
                                   threshold = 0.5,
                                   save_plots = FALSE,
                                   output_dir = "plots") {

  # Calculate missing value statistics
  missing_stats <- clinical_data %>%
    summarise_all(~sum(is.na(.))) %>%
    gather(key = "Feature", value = "Missing_Count") %>%
    mutate(
      Missing_Percentage = round((Missing_Count / nrow(clinical_data)) * 100, 2),
      Total_Patients = nrow(clinical_data)
    ) %>%
    arrange(desc(Missing_Percentage))
  
  # Create directory for saving plots if needed
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Plot 1: Bar plot of missing percentages
  p1 <- ggplot(missing_stats, aes(x = reorder(Feature, Missing_Percentage), 
                                  y = Missing_Percentage)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = threshold*100, color = "red", linetype = "dashed", 
               size = 1, alpha = 0.8) +
    coord_flip() +
    labs(title = "Missing Values by Clinical Feature",
         subtitle = paste("Red line indicates", threshold*100, "% threshold"),
         x = "Clinical Features",
         y = "Missing Percentage (%)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 14, face = "bold"),
          panel.background = element_rect(fill = 'white'),
          plot.background = element_rect(fill = 'white'))
  
  
  if (save_plots) {
    ggsave(filename = file.path(output_dir, "missing_values_barplot.png"), 
           plot = p1, width = 12, height = 8, dpi = 300)
  }
  
  # Plot 2: Heatmap of missing values pattern
  
  # Sample patients if dataset is too large for visualization
  sampled_data <- clinical_data
  
  # Create missing values heatmap
  missing_pattern <- sampled_data %>%
    mutate(Patient_ID = row_number()) %>%
    gather(key = "Feature", value = "Value", -Patient_ID) %>%
    mutate(Missing = is.na(Value))
  
  p2 <- ggplot(missing_pattern, aes(x = Feature, y = Patient_ID, fill = Missing)) +
    geom_tile() +
    scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "red"),
                      labels = c("Present", "Missing"),
                      name = "Data Status") +
    labs(title = "Missing Values Pattern Heatmap",
         subtitle = "Each row represents a patient, each column a clinical feature",
         x = "Clinical Features",
         y = "Patients") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 14, face = "bold"),
          panel.background = element_rect(fill = 'white'),
          plot.background = element_rect(fill = 'white'))
  
  if (save_plots) {
    ggsave(filename = file.path(output_dir, "missing_values_heatmap.png"), 
           plot = p2, width = 12, height = 8, dpi = 300)
  }

}