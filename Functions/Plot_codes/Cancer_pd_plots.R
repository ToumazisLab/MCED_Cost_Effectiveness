# Define the plot_cancer_costs function
plot_cancer_costs <- function(final_result_list, cancer_list, screening_cancers, output_file = NULL) {
  
  ordered_cancer_list <- c(
    screening_cancers,
    setdiff(cancer_list, screening_cancers)
  )
  
  # Helper function to lighten a color by blending it with white
  lighten_color <- function(color, factor = 0.5) {
    # Ensure the color is valid
    if (is.na(col2rgb(color, TRUE)[1,1])) {
      stop(paste("Invalid color:", color))
    }
    rgb_vals <- col2rgb(color)
    # Blend with white
    new_rgb <- rgb_vals + (255 - rgb_vals) * factor
    # Convert back to hexadecimal color
    rgb(new_rgb[1,]/255, new_rgb[2,]/255, new_rgb[3,]/255)
  }
  
  # Validate inputs
  if (!is.list(final_result_list)) {
    stop("final_result_list must be a list.")
  }
  
  if (!is.character(ordered_cancer_list) || length(ordered_cancer_list) == 0) {
    stop("ordered_cancer_list must be a non-empty character vector.")
  }
  
  if (!is.character(screening_cancers)) {
    stop("screening_cancers must be a character vector.")
  }
  
  # Define arms and stages
  arms <- c("MCED", "usual")
  stages <- c("I", "II", "III", "IV")
  
  # Initialize an empty list to store rows (for performance)
  cost_list <- list()
  row_index <- 1
  
  # Iterate through each arm, stage, and cancer to extract and sum costs
  for (arm in arms) {
    for (stage in stages) {
      for (cancer in ordered_cancer_list) {
        if (arm == "MCED") {
          # Define suffixes for MCED
          suffixes <- c("TP", "FN", "NA")
          for (suffix in suffixes) {
            # Construct the path to the cost value
            # Using tryCatch to handle missing data gracefully
            cost_value <- tryCatch({
              # Navigate through final_result_list
              cost_vector <- final_result_list[["MCED"]][[cancer]][["MCED_Cancer_cases"]][[paste0(stage, "_post_cost_", suffix)]]
              # Sum the cost_vector if it's numeric
              if (is.numeric(cost_vector)) {
                sum(cost_vector, na.rm = TRUE)
              } else {
                NA
              }
            }, error = function(e) {
              NA
            })
            
            # Store the row in the list
            cost_list[[row_index]] <- list(
              Arm = arm,
              Stage = stage,
              Cancer = cancer,
              Suffix = suffix,
              Cost = ifelse(is.null(cost_value) || is.na(cost_value), 0, cost_value)
            )
            row_index <- row_index + 1
          }
        } else if (arm == "usual") {
          # For Usual Care, the suffix is "Total"
          suffix <- "Total"
          cost_value <- tryCatch({
            cost_vector <- final_result_list[["usual"]][[cancer]][["Cancer_cases_screening"]][[paste0(stage, "_post_cost")]]
            if (is.numeric(cost_vector)) {
              sum(cost_vector, na.rm = TRUE)
            } else {
              NA
            }
          }, error = function(e) {
            NA
          })
          
          # Store the row in the list
          cost_list[[row_index]] <- list(
            Arm = arm,
            Stage = stage,
            Cancer = cancer,
            Suffix = suffix,
            Cost = ifelse(is.null(cost_value) || is.na(cost_value), 0, cost_value)
          )
          row_index <- row_index + 1
        }
      }
    }
  }
  
  # Convert the list to a dataframe
  cost_data <- bind_rows(cost_list)
  
  # Replace any remaining NA costs with 0
  cost_data$Cost[is.na(cost_data$Cost)] <- 0
  
  # Create a new column for Cancer and Suffix separately
  cost_data <- cost_data %>%
    mutate(
      Cancer_Suffix = Suffix
    ) %>% 
    mutate(Cancer = factor(Cancer, levels = ordered_cancer_list))
  
  # Assign base colors using RColorBrewer's palette
  num_cancers <- length(ordered_cancer_list)
  
  # Determine the palette to use based on the number of cancers
  if (num_cancers <= 12) {
    base_palette <- brewer.pal(12, "Set3")[1:num_cancers]
  } else {
    # If more than 12 cancers, generate a larger palette
    base_palette <- colorRampPalette(brewer.pal(12, "Set3"))(num_cancers)
  }
  
  # Assign base colors to cancers
  names(base_palette) <- ordered_cancer_list
  
  # Define alpha levels for suffixes
  # TP: 0% transparency (fully opaque) => alpha = 1
  # FN: 25% transparency          => alpha = 0.75
  # NA: 50% transparency          => alpha = 0.5
  # Total: 0% transparency        => alpha = 1.0 (for Usual Care)
  suffix_alpha <- c("TP" = 1.0, "FN" = 0.75, "NA" = 0.5, "Total" = 1.0)
  
  # Summarize costs per Cancer, Arm, Stage, Suffix
  cost_summary <- cost_data %>%
    group_by(Cancer, Arm, Stage, Suffix) %>%
    summarize(Total_Cost = sum(Cost), .groups = 'drop')
  
  max_cost <- max(cost_summary$Total_Cost, na.rm = TRUE)
  
  padding_factor <- 0.15  # 10% padding
  y_upper_limit <- max_cost * (1 + padding_factor)
  
  # Handle cases where max_cost is 0 to avoid y_upper_limit being 0
  if (y_upper_limit == 0) {
    y_upper_limit <- 1  # Set a minimal upper limit
  }
  
  # Create label strings per Cancer
  label_strings <- cost_summary %>%
    group_by(Cancer, Arm) %>%
    arrange(Stage, Suffix) %>%
    mutate(
      Stage_Suffix_Cost = paste0(Stage, ": ", Suffix, " - ", sprintf("%.2e", Total_Cost))
    ) %>%
    summarize(
      Arm_Info = paste(Stage_Suffix_Cost, collapse = "; "),
      .groups = 'drop'
    ) %>%
    group_by(Cancer) %>%
    summarize(
      Label = paste0(Arm, ": ", Arm_Info, collapse = "; "),
      .groups = 'drop'
    ) %>%
    mutate(
      Label = str_replace_all(Label, "; ", "\n")
    )
  
  # Create a named vector for labels
  label_mapping <- setNames(
    paste0(ordered_cancer_list, ": ", label_strings$Label),  # Append cancer name to label
    ordered_cancer_list
  )
  
  # Create the stacked bar plot with fill = Cancer and alpha = Suffix
  plot <- ggplot(cost_summary, aes(x = Stage, y = Total_Cost, fill = Cancer, alpha = Suffix)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Arm, ncol = 1, scales = "free_y") +
    scale_fill_manual(
      values = base_palette,
      name = "Cancer",
      breaks = ordered_cancer_list,      # Ensure legend matches reordered list
      labels = label_mapping       # Ensure legend matches reordered list
    ) +
    scale_alpha_manual(
      values = suffix_alpha,
      name = "Category",
      breaks = names(suffix_alpha),
      labels = names(suffix_alpha)
    ) +
    labs(
      title = "Stacked Cancer Costs by Stage and Arm",
      x = "Stage",
      y = "Cost"
    ) +
    # Apply dynamic y-axis limits
    coord_cartesian(ylim = c(0, y_upper_limit)) +
    theme_minimal(base_size = 12) +
    theme(
      text = element_text(size = 15),
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 15)
    )
  
  
  # If output_file is specified, save the plot
  if (!is.null(output_file)) {
    ggsave(filename = output_file, plot = plot, width = 30, height = 15, dpi = 300)
    message(paste("Plot successfully saved to", output_file))
  }
  
  return(plot)
}

plot_cancer_utility <- function(final_result_list, cancer_list, screening_cancers, output_file = NULL) {
  
  ordered_cancer_list <- c(
    screening_cancers,
    setdiff(cancer_list, screening_cancers)
  )
  
  # Helper function to lighten a color by blending it with white
  lighten_color <- function(color, factor = 0.5) {
    # Ensure the color is valid
    if (is.na(col2rgb(color, TRUE)[1,1])) {
      stop(paste("Invalid color:", color))
    }
    rgb_vals <- col2rgb(color)
    # Blend with white
    new_rgb <- rgb_vals + (255 - rgb_vals) * factor
    # Convert back to hexadecimal color
    rgb(new_rgb[1,]/255, new_rgb[2,]/255, new_rgb[3,]/255)
  }
  
  # Validate inputs
  if (!is.list(final_result_list)) {
    stop("final_result_list must be a list.")
  }
  
  if (!is.character(ordered_cancer_list) || length(ordered_cancer_list) == 0) {
    stop("ordered_cancer_list must be a non-empty character vector.")
  }
  
  if (!is.character(screening_cancers)) {
    stop("screening_cancers must be a character vector.")
  }
  
  # Define arms and stages
  arms <- c("MCED", "usual")
  stages <- c("I", "II", "III", "IV")
  
  # Initialize an empty list to store rows (for performance)
  cost_list <- list()
  row_index <- 1
  
  # Iterate through each arm, stage, and cancer to extract and sum costs
  for (arm in arms) {
    for (stage in stages) {
      for (cancer in ordered_cancer_list) {
        if (arm == "MCED") {
          # Define suffixes for MCED
          suffixes <- c("TP", "FN", "NA")
          for (suffix in suffixes) {
            # Construct the path to the cost value
            # Using tryCatch to handle missing data gracefully
            cost_value <- tryCatch({
              # Navigate through final_result_list
              cost_vector <- final_result_list[["MCED"]][[cancer]][["MCED_Cancer_cases"]][[paste0(stage, "_post_utility_", suffix)]]
              # Sum the cost_vector if it's numeric
              if (is.numeric(cost_vector)) {
                sum(cost_vector, na.rm = TRUE)
              } else {
                NA
              }
            }, error = function(e) {
              NA
            })
            
            # Store the row in the list
            cost_list[[row_index]] <- list(
              Arm = arm,
              Stage = stage,
              Cancer = cancer,
              Suffix = suffix,
              Cost = ifelse(is.null(cost_value) || is.na(cost_value), 0, cost_value)
            )
            row_index <- row_index + 1
          }
        } else if (arm == "usual") {
          # For Usual Care, the suffix is "Total"
          suffix <- "Total"
          cost_value <- tryCatch({
            cost_vector <- final_result_list[["usual"]][[cancer]][["Cancer_cases_screening"]][[paste0(stage, "_post_utility")]]
            if (is.numeric(cost_vector)) {
              sum(cost_vector, na.rm = TRUE)
            } else {
              NA
            }
          }, error = function(e) {
            NA
          })
          
          # Store the row in the list
          cost_list[[row_index]] <- list(
            Arm = arm,
            Stage = stage,
            Cancer = cancer,
            Suffix = suffix,
            Cost = ifelse(is.null(cost_value) || is.na(cost_value), 0, cost_value)
          )
          row_index <- row_index + 1
        }
      }
    }
  }
  
  # Convert the list to a dataframe
  cost_data <- bind_rows(cost_list)
  
  # Replace any remaining NA costs with 0
  cost_data$Cost[is.na(cost_data$Cost)] <- 0
  
  # Create a new column for Cancer and Suffix separately
  cost_data <- cost_data %>%
    mutate(
      Cancer_Suffix = Suffix
    ) %>% 
    mutate(Cancer = factor(Cancer, levels = ordered_cancer_list))
  
  # Assign base colors using RColorBrewer's palette
  num_cancers <- length(ordered_cancer_list)
  
  # Determine the palette to use based on the number of cancers
  if (num_cancers <= 12) {
    base_palette <- brewer.pal(12, "Set3")[1:num_cancers]
  } else {
    # If more than 12 cancers, generate a larger palette
    base_palette <- colorRampPalette(brewer.pal(12, "Set3"))(num_cancers)
  }
  
  # Assign base colors to cancers
  names(base_palette) <- ordered_cancer_list
  
  # Define alpha levels for suffixes
  # TP: 0% transparency (fully opaque) => alpha = 1
  # FN: 25% transparency          => alpha = 0.75
  # NA: 50% transparency          => alpha = 0.5
  # Total: 0% transparency        => alpha = 1.0 (for Usual Care)
  suffix_alpha <- c("TP" = 1.0, "FN" = 0.75, "NA" = 0.5, "Total" = 1.0)
  
  # Summarize costs per Cancer, Arm, Stage, Suffix
  cost_summary <- cost_data %>%
    group_by(Cancer, Arm, Stage, Suffix) %>%
    summarize(Total_Cost = sum(Cost), .groups = 'drop')
  
  max_cost <- max(cost_summary$Total_Cost, na.rm = TRUE)
  
  padding_factor <- 0.15  # 10% padding
  y_upper_limit <- max_cost * (1 + padding_factor)
  
  # Handle cases where max_cost is 0 to avoid y_upper_limit being 0
  if (y_upper_limit == 0) {
    y_upper_limit <- 1  # Set a minimal upper limit
  }
  
  # Create label strings per Cancer
  label_strings <- cost_summary %>%
    group_by(Cancer, Arm) %>%
    arrange(Stage, Suffix) %>%
    mutate(
      Stage_Suffix_Cost = paste0(Stage, ": ", Suffix, " - ", sprintf("%.2e", Total_Cost))
    ) %>%
    summarize(
      Arm_Info = paste(Stage_Suffix_Cost, collapse = "; "),
      .groups = 'drop'
    ) %>%
    group_by(Cancer) %>%
    summarize(
      Label = paste0(Arm, ": ", Arm_Info, collapse = "; "),
      .groups = 'drop'
    ) %>%
    mutate(
      Label = str_replace_all(Label, "; ", "\n")
    )
  
  # Create a named vector for labels
  label_mapping <- setNames(
    paste0(ordered_cancer_list, ": ", label_strings$Label),  # Append cancer name to label
    ordered_cancer_list
  )
  
  # Create the stacked bar plot with fill = Cancer and alpha = Suffix
  plot <- ggplot(cost_summary, aes(x = Stage, y = Total_Cost, fill = Cancer, alpha = Suffix)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Arm, ncol = 1, scales = "fixed") +  # Use fixed scales
    scale_fill_manual(
      values = base_palette,
      name = "Cancer",
      breaks = ordered_cancer_list,      # Ensure legend matches reordered list
      labels = label_mapping             # Ensure legend matches reordered list
    ) +
    scale_alpha_manual(
      values = suffix_alpha,
      name = "Category",
      breaks = names(suffix_alpha),
      labels = names(suffix_alpha)
    ) +
    labs(
      title = "Stacked Cancer Post Dx Utility by Stage and Arm",
      x = "Stage",
      y = "Utility"
    ) +
    # Apply dynamic y-axis limits
    coord_cartesian(ylim = c(0, y_upper_limit)) +
    theme_minimal(base_size = 12) +
    theme(
      text = element_text(size = 15),
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 15)
    )
  
  
  # If output_file is specified, save the plot
  if (!is.null(output_file)) {
    ggsave(filename = output_file, plot = plot, width = 30, height = 15, dpi = 300)
    message(paste("Plot successfully saved to", output_file))
  }
  
  return(plot)
}
