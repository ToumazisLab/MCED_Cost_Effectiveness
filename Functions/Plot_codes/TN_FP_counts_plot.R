FP_counts_plot <- function(final_result_list, screening_cancers, 
                           cancer_list) {
  # Define ordered_cancer_list
  ordered_cancer_list <- c(
    screening_cancers,
    setdiff(cancer_list, screening_cancers)
  )
  
  # Assign colors using RColorBrewer's Set3 palette
  num_cancers <- length(ordered_cancer_list)
  
  if (num_cancers <= 12) {
    base_palette <- brewer.pal(12, "Set3")[1:num_cancers]
  } else {
    base_palette <- colorRampPalette(brewer.pal(12, "Set3"))(num_cancers)
  }
  
  # Assign names to the palette for easy reference
  names(base_palette) <- ordered_cancer_list
  
  # Function to process MCED data
  process_mced_data <- function(final_result_list, ordered_cancer_list, screening_cancers) {
    mced_data <- lapply(ordered_cancer_list, function(cancer) {
      mced_entry <- final_result_list[["MCED"]][[cancer]][["mced_TN_FP"]]
      
      # Base subcategories
      FP_fu <- sum(mced_entry[["FP_fu_compliant"]])
      nc_FP <- sum(mced_entry[["nc_FP"]])
      nc_TN <- sum(mced_entry[["nc_TN"]])
      FP_i <- sum(mced_entry[["FP_i"]])
      
      # Initialize a list with base subcategories
      data <- list(
        FP_fu = FP_fu,
        FP_i = FP_i,
        nc_FP = nc_FP,
        nc_Invasive = nc_TN
      )
      
      # Add conditional subcategories for screening cancers
      if (cancer %in% screening_cancers) {
        FP_FP <- sum(mced_entry[["FP_FP"]])
        FP_FP_i <- sum(mced_entry[["FP_FP_i"]])
        TN_FP <- sum(mced_entry[["TN_FP"]])
        TN_FP_i <- sum(mced_entry[["TN_FP_i"]])
        
        data$FP_FP <- FP_FP
        data$FP_FP_i <- FP_FP_i
        data$TN_FP <- TN_FP
        data$TN_FP_i <- TN_FP_i
      }
      
      # Convert to a data frame
      data_frame <- as.data.frame(as.list(data))
      data_frame$cancer <- cancer
      return(data_frame)
    })
    
    # Combine all cancers into one data frame
    mced_df <- bind_rows(mced_data)
    
    # Replace NA with 0 for non-applicable subcategories
    mced_df[is.na(mced_df)] <- 0
    
    # Summarize by cancer and subcategory
    mced_long <- mced_df %>%
      pivot_longer(
        cols = -cancer,
        names_to = "subcategory",
        values_to = "value"
      ) %>%
      group_by(cancer, subcategory) %>%
      summarise(value = sum(value), .groups = 'drop')
    
    return(mced_long)
  }
  
  # Function to process Usual data
  process_usual_data <- function(final_result_list, ordered_cancer_list) {
    usual_data <- lapply(ordered_cancer_list, function(cancer) {
      usual_entry <- final_result_list[["usual"]][[cancer]][["FP_TN_Screening_Summary"]]
      
      # Extract subcategories
      FP <- sum(usual_entry[["FP"]])
      FP_i <- sum(usual_entry[["Invasive"]])

      # Create a data frame
      data_frame <- data.frame(
        FP = FP,
        FP_i = FP_i,
        cancer = cancer,
        stringsAsFactors = FALSE
      )
      
      return(data_frame)
    })
    
    # Combine all cancers into one data frame
    usual_df <- bind_rows(usual_data)
    
    # Replace NA with 0
    usual_df[is.na(usual_df)] <- 0
    
    # Summarize by cancer and subcategory
    usual_long <- usual_df %>%
      pivot_longer(
        cols = -cancer,
        names_to = "subcategory",
        values_to = "value"
      ) %>%
      group_by(cancer, subcategory) %>%
      summarise(value = sum(value), .groups = 'drop')
    
    return(usual_long)
  }
  
  # Process data for both arms
  mced_long <- process_mced_data(final_result_list, ordered_cancer_list, screening_cancers)
  usual_long <- process_usual_data(final_result_list, ordered_cancer_list)
  
  # Ensure the cancer factor is ordered correctly
  mced_long$cancer <- factor(mced_long$cancer, levels = ordered_cancer_list)
  usual_long$cancer <- factor(usual_long$cancer, levels = ordered_cancer_list)
  
  # Create MCED palette by repeating base_palette for subcategories
  # Define a palette for MCED subcategories
  mced_subcategories <- unique(mced_long$subcategory)
  num_subcategories_mced <- length(mced_subcategories)
  
  # Assign distinct colors for MCED subcategories
  mced_palette <- brewer.pal(max(3, num_subcategories_mced), "Paired")[1:num_subcategories_mced]
  names(mced_palette) <- mced_subcategories
  
  mced_max <- max(mced_long$value, na.rm = TRUE)
  usual_max <- max(usual_long$value, na.rm = TRUE)
  
  # Define a padding factor (e.g., 10%)
  padding_factor <- 0.15
  
  # Calculate the upper limit with padding
  mced_y_limit <- mced_max * (1 + padding_factor)
  usual_y_limit <- usual_max * (1 + padding_factor)
  
  
  # Create MCED plot
  mced_plot <- ggplot(mced_long, aes(x = cancer, y = value, fill = subcategory)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = mced_palette) +
    labs(
      title = "MCED Arm",
      x = "Cancer Type",
      y = "Count",
      fill = "Subcategory"
    ) +
    coord_cartesian(ylim = c(0, mced_y_limit)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Define a palette for Usual subcategories
  usual_subcategories <- unique(usual_long$subcategory)
  num_subcategories_usual <- length(usual_subcategories)
  
  # Assign distinct colors for Usual subcategories
  usual_palette <- brewer.pal(max(3, num_subcategories_usual), "Set2")[1:num_subcategories_usual]
  names(usual_palette) <- usual_subcategories
  
  # Create Usual plot
  usual_plot <- ggplot(usual_long, aes(x = cancer, y = value, fill = subcategory)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = usual_palette) +
    labs(
      title = "Usual Arm",
      x = "Cancer Type",
      y = "Count",
      fill = "Subcategory"
    ) +
    coord_cartesian(ylim = c(0, usual_y_limit)) + 
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )

  
  plot <- grid.arrange(mced_plot, usual_plot, ncol = 2)
  
  # Optionally, return the plots as a list
  return(list(mced_long = mced_long, usual_long = usual_long, 
              plot = plot, MCED_Plot = mced_plot, Usual_Plot = usual_plot))
}

