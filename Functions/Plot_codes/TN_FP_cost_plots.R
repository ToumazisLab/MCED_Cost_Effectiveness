create_fp_tn_plots <- function(final_result_list, 
                               screening_cancers, 
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
  
  # Define cost categories
  cost_categories <- c("Imaging", "Diagnostic", "Invasive")
  arms <- c("MCED", "Usual Care")
  
  # Initialize empty data frames for FP and TN
  fp_data <- data.frame(
    Category = character(),
    Cancer = character(),
    Cost = numeric(),
    Cost_Type = character(),
    stringsAsFactors = FALSE
  )
  
  tn_data <- data.frame(
    Category = character(),
    Cancer = character(),
    Cost = numeric(),
    Cost_Type = character(),
    stringsAsFactors = FALSE
  )
  
  # Iterate over each cancer to populate FP and TN data
  for (cancer in ordered_cancer_list) {
    for (category in cost_categories) {
      # MCED Arm FP Costs
      mced_FP <- sum(final_result_list[["MCED"]][[cancer]][["mced_tn_fp_cost"]][[paste0("mced_FP_", category)]], na.rm = TRUE)
      # MCED Arm nc Costs
      nc_FP_key <- paste0("nc_FP_", tolower(category), "_cost")
      nc_FP <- ifelse(!is.null(final_result_list[["MCED"]][[cancer]][["mced_tn_fp_cost"]][[nc_FP_key]]),
                      sum(final_result_list[["MCED"]][[cancer]][["mced_tn_fp_cost"]][[nc_FP_key]], na.rm = TRUE),
                      0)
      
      # MCED Arm TN Costs
      mced_TN_key <- paste0("mced_TN_", category)
      mced_TN <- sum(final_result_list[["MCED"]][[cancer]][["mced_tn_fp_cost"]][[mced_TN_key]], na.rm = TRUE)
      # MCED Arm nc TN Costs
      nc_TN_key <- paste0("nc_TN_", tolower(category), "_cost")
      nc_TN <- ifelse(!is.null(final_result_list[["MCED"]][[cancer]][["mced_tn_fp_cost"]][[nc_TN_key]]),
                      sum(final_result_list[["MCED"]][[cancer]][["mced_tn_fp_cost"]][[nc_TN_key]], na.rm = TRUE),
                      0)
      
      # Usual Care Arm FP Costs
      if (category == "Invasive") {
        usual_FP <- sum(final_result_list[["usual"]][[cancer]][["FP_TN_Screening_Summary"]][["FP_Invasive_cost"]], na.rm = TRUE)
      } else {
        usual_FP_key <- paste0("FP_", tolower(category), "_cost")
        usual_FP <- sum(final_result_list[["usual"]][[cancer]][["FP_TN_Screening_Summary"]][[usual_FP_key]], na.rm = TRUE)
      }
      
      # Usual Care Arm TN Costs
      usual_TN_key <- paste0("TN_", tolower(category), "_cost")
      usual_TN <- sum(final_result_list[["usual"]][[cancer]][["FP_TN_Screening_Summary"]][[usual_TN_key]], na.rm = TRUE)
      
      # Populate FP Data Frame for MCED Arm
      fp_data <- rbind(fp_data, data.frame(
        Category = paste(category, "MCED"),
        Cancer = cancer,
        Cost = mced_FP,
        Cost_Type = "mced",
        stringsAsFactors = FALSE
      ))
      
      # Populate FP Data Frame for MCED nc Costs if applicable
      if (nc_FP > 0) {
        fp_data <- rbind(fp_data, data.frame(
          Category = paste(category, "MCED"),
          Cancer = cancer,
          Cost = nc_FP,
          Cost_Type = "nc",
          stringsAsFactors = FALSE
        ))
      }
      
      # Populate TN Data Frame for MCED Arm
      tn_data <- rbind(tn_data, data.frame(
        Category = paste(category, "MCED"),
        Cancer = cancer,
        Cost = mced_TN,
        Cost_Type = "mced",
        stringsAsFactors = FALSE
      ))
      
      # Populate TN Data Frame for MCED nc Costs if applicable
      if (nc_TN > 0) {
        tn_data <- rbind(tn_data, data.frame(
          Category = paste(category, "MCED"),
          Cancer = cancer,
          Cost = nc_TN,
          Cost_Type = "nc",
          stringsAsFactors = FALSE
        ))
      }
      
      # Populate FP Data Frame for Usual Care Arm
      fp_data <- rbind(fp_data, data.frame(
        Category = paste(category, "Usual"),
        Cancer = cancer,
        Cost = usual_FP,
        Cost_Type = "Usual",
        stringsAsFactors = FALSE
      ))
      
      # Populate TN Data Frame for Usual Care Arm
      tn_data <- rbind(tn_data, data.frame(
        Category = paste(category, "Usual"),
        Cancer = cancer,
        Cost = usual_TN,
        Cost_Type = "Usual",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Set Cost_Type as factor with defined levels
  fp_data$Cost_Type <- factor(fp_data$Cost_Type, levels = c("mced", "nc", "Usual"))
  tn_data$Cost_Type <- factor(tn_data$Cost_Type, levels = c("mced", "nc", "Usual"))
  
  # Define alpha levels for Cost_Type
  alpha_levels <- c("mced" = 1, "nc" = 0.5, "Usual" = 1)
  
  # Ensure Category is a factor with a defined order
  categories_order <- c()
  for (category in cost_categories) {
    categories_order <- c(categories_order, paste(category, "MCED"), paste(category, "Usual"))
  }
  
  fp_data$Category <- factor(fp_data$Category, levels = categories_order)
  tn_data$Category <- factor(tn_data$Category, levels = categories_order)
  
  # Ensure Cancer is a factor with ordered levels
  fp_data$Cancer <- factor(fp_data$Cancer, levels = ordered_cancer_list)
  tn_data$Cancer <- factor(tn_data$Cancer, levels = ordered_cancer_list)
  
  # Calculate the maximum Cost across both datasets
  max_fp_cost <- max(fp_data$Cost, na.rm = TRUE)
  max_tn_cost <- max(tn_data$Cost, na.rm = TRUE)
  
  # Define a padding factor (e.g., 10%)
  padding_factor <- 0.10
  
  # Calculate the upper limit with padding
  y_fp_upper_limit <- max_fp_cost * (1 + padding_factor)
  y_tn_upper_limit <- max_tn_cost * (1 + padding_factor)
  
  
  # Create FP Plot
  fp_plot <- ggplot(fp_data, aes(x = Category, y = Cost, fill = Cancer, alpha = Cost_Type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = base_palette) +
    scale_alpha_manual(values = alpha_levels, guide = "none") +
    labs(
      title = "False Positives (FP) Costs by Category and Cancer",
      x = "FP Category",
      y = "Total FP Cost",
      fill = "Cancer"
    ) +
    coord_cartesian(ylim = c(0, y_fp_upper_limit)) +
    theme_minimal() +
    theme(
      text = element_text(size = 14),
      legend.position = "right",
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1), # X-axis text size
      axis.text.y = element_text(size = 14), # Y-axis text size
      axis.title.x = element_text(size = 14), # X-axis title size
      axis.title.y = element_text(size = 14) # Y-axis title size
    )
  
  # Create TN Plot
  tn_plot <- ggplot(tn_data, aes(x = Category, y = Cost, fill = Cancer, alpha = Cost_Type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = base_palette) +
    scale_alpha_manual(values = alpha_levels, guide = "none") +
    labs(
      title = "True Negatives (TN) Costs by Category and Cancer",
      x = "TN Category",
      y = "Total TN Cost",
      fill = "Cancer"
    ) +
    coord_cartesian(ylim = c(0, y_tn_upper_limit)) +
    theme_minimal() +
    theme(
      text = element_text(size = 14),
      legend.position = "right",
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1), # X-axis text size
      axis.text.y = element_text(size = 14), # Y-axis text size
      axis.title.x = element_text(size = 14), # X-axis title size
      axis.title.y = element_text(size = 14) # Y-axis title size
    )
  
  # Arrange the two plots vertically
  combined_plot <- grid.arrange(fp_plot, tn_plot, ncol = 1)
  
  # Return the combined plot
  return(list(combined_plot = combined_plot, fp_data = fp_data, tn_data = tn_data))
}
