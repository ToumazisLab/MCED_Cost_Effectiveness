library(dplyr)
library(ggplot2)
library(patchwork)
library(ggtext)


# Function to create a single stacked bar plot
update_results_csv <- function(final_result_list, 
                               new_params, 
                               cancer_list, 
                               file_path, 
                               run_id) {
  
  library(dplyr)
  library(tidyr)
  
  # Step 1: Extract incremental_values
  incremental_values <- final_result_list[["Overall"]][["final_icer_result"]] %>%
    filter(Variable %in% c("Incremental_Cost", "Incremental_Utility", "ICER")) %>%
    select(Variable, Value) %>%
    distinct() %>%  # Ensure uniqueness
    pivot_wider(names_from = Variable, values_from = Value)
  
  # Check if incremental_values has exactly one row
  if(nrow(incremental_values) != 1){
    stop("Error: 'incremental_values' should have exactly one row after pivot_wider.")
  }
  
  # Step 2: Create results_df
  results_df <- final_result_list[["Overall"]][["final_icer_result"]] %>%
    select(-Type, -State) %>%  # Drop `Type` and `State` columns
    arrange(Variable, Arm) %>%  # Arrange by Variable first, then by Arm
    distinct(Arm, Variable, .keep_all = TRUE) %>%  # Ensure unique combinations
    unite("Arm_Variable", Arm, Variable, sep = "_", remove = TRUE)
  
  # Step 3: Pivot to wide format
  results_wide <- results_df %>%
    pivot_wider(names_from = Arm_Variable, values_from = Value)
  
  # Check if results_wide has exactly one row
  if(nrow(results_wide) != 1){
    stop("Error: 'results_wide' should have exactly one row after pivot_wider.")
  }
  
  # Step 4: Convert new_params to a single-row DataFrame
  params_df <- as.data.frame(new_params)
  
  # If new_params is a list of parameters, ensure it's a single row
  if(nrow(params_df) != 1){
    params_df <- params_df %>%
      summarise(across(everything(), ~ first(.)))  # Modify aggregation as needed
  }
  
  # Step 5: Add incremental_values to params_df
  params_df <- params_df %>%
    mutate(
      Incremental_cost = incremental_values$Incremental_Cost,
      Incremental_utility = incremental_values$Incremental_Utility,
      ICER = incremental_values$ICER
    )
  
  # Step 6: Create the `cancer_list` column
  cancer_list_sorted <- paste(sort(cancer_list), collapse = " | ")
  
  # Step 7: Combine results and parameters ensuring single row
  final_row <- results_wide %>%
    mutate(cancer_list = cancer_list_sorted, run_id = run_id) %>%
    bind_cols(params_df)
  
  # Ensure final_row has only one row
  if(nrow(final_row) != 1){
    stop("Error: 'final_row' should have exactly one row.")
  }
  
  # Reorder columns to place `cancer_list` and `run_id` first
  final_row <- final_row %>%
    select(cancer_list, run_id, everything())
  
  # Step 8: Check if the file exists and append or create
  if (file.exists(file_path)) {
    # Load existing CSV
    existing_data <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Standardize column types
    for (col in names(existing_data)) {
      if (col %in% names(final_row)) {
        target_class <- class(existing_data[[col]])
        # Convert final_row's column to match existing_data's column class
        if(target_class == "numeric"){
          final_row[[col]] <- as.numeric(final_row[[col]])
        } else if(target_class == "integer"){
          final_row[[col]] <- as.integer(final_row[[col]])
        } else if(target_class == "character"){
          final_row[[col]] <- as.character(final_row[[col]])
        } else if(target_class == "factor"){
          final_row[[col]] <- as.factor(final_row[[col]])
        }
        # Add more type conversions if necessary
      }
    }
    
    # Append the new row
    updated_data <- bind_rows(existing_data, final_row)
  } else {
    # If the file doesn't exist, create a new file
    updated_data <- final_row
  }
  
  # Step 9: Save to CSV
  write.csv(updated_data, file_path, row.names = FALSE)
  
  # Return the updated DataFrame for inspection
  return(updated_data)
}


make_single_bar_plot <- function(df, arm_label, type_label) {
  
  # Get y axis limit
  sum_values <- tapply(df$Value, list(df$Arm, df$Type), sum, na.rm = TRUE)
  
  # Finding the max value for each Type
  max_values <- apply(sum_values, 2, max, na.rm = TRUE)

  cost_ylim <- as.numeric(max_values[1]) * (1.15)
  util_ylim <- as.numeric(max_values[2]) * (1.15)
  # Subset the data for the chosen Arm & Type
  sub_df <- df %>%
    filter(Arm == arm_label, Type == type_label)
  
  # Apply the appropriate order for Variable
  if (arm_label == "mced" && type_label == "cost") {
    sub_df <- sub_df %>%
      mutate(Variable = factor(Variable, levels = mced_cost_order))
    
    y_limit <- cost_ylim
    
  } else if (arm_label == "mced" && type_label == "utility") {
    sub_df <- sub_df %>%
      mutate(Variable = factor(Variable, levels = mced_utility_order))
    
    y_limit <- util_ylim
  } else if (arm_label == "usual" && type_label == "utility") {
    sub_df <- sub_df %>%
      mutate(Variable = factor(Variable, levels = usual_utility_order))
    
    y_limit <- util_ylim
    
  } else if (arm_label == "usual" && type_label == "cost") {
    sub_df <- sub_df %>%
      mutate(Variable = factor(Variable, levels = usual_cost_order))
    
    y_limit <- cost_ylim
    
  } else {
    # Fallback for unexpected Arm/Type combinations
    sub_df <- sub_df %>%
      mutate(Variable = factor(Variable, levels = unique(Variable)))
  }
  
  # Calculate formatted labels with values in scientific notation
  variable_labels <- sub_df %>%
    group_by(Variable) %>%
    summarise(TotalValue = sum(Value, na.rm = TRUE)) %>%
    mutate(
      Label = paste0(Variable, " (", formatC(TotalValue, format = "e", digits = 2), ")")
    ) %>%
    select(Variable, Label) %>%
    deframe()
  
  # Calculate total value for the entire bar
  total_value <- sum(sub_df$Value, na.rm = TRUE)
  formatted_total <- formatC(total_value, format = "e", digits = 3)
  
  # Create the single stacked bar plot
  p <- ggplot(sub_df, aes(x = "", y = Value, fill = Variable)) +
    geom_col(position = "stack") +
    # Add the total value as annotation on top of the bar
    annotate(
      "text",
      x = 1,  # Position at the center of the bar
      y = total_value, 
      label = formatted_total,
      vjust = -0.5, # Position above the bar
      size = 8      # Text size
    ) +
    # Apply custom color mapping
    scale_fill_manual(
      values = var_colors, 
      labels = variable_labels,
      drop = FALSE
    ) +
    labs(
      title = paste0(toupper(substring(arm_label, 1,1)), substring(arm_label, 2), " - ", 
                     toupper(substring(type_label, 1,1)), substring(type_label, 2)),
      x     = NULL,
      y     = "Value",
      fill  = "Variable"
    ) +
    # Apply dynamic y-axis limits
    coord_cartesian(ylim = c(0, y_limit)) +
    # Minimal theme and remove x-axis text
    theme_minimal() +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title   = element_text(size = 20, face = "bold", hjust = 0.5),    # Centered and bold title
      axis.title.y = element_text(size = 12),    # Adjusted y-axis label size
      axis.text.y  = element_text(size = 10),    # Adjusted y-axis text size
      legend.text  = element_text(size = 12),     # Adjusted legend text size
      legend.title = element_text(size = 12),     # Adjusted legend title size
      legend.key.size = unit(0.5, "cm"),         # Smaller color boxes
      legend.position = "right",                  # Positioned legend to the right
      legend.box = "vertical"                     # Arrange legend items vertically
    )
  
  # Return the plot
  return(p)
}


generate_unique_filename <- function(base_path, plot_type, extension) {
  # Ensure the base path ends with a slash
  if (!grepl("/$", base_path)) {
    base_path <- paste0(base_path, "/")
  }
  
  # Construct the base filename using plot_type
  base_name <- paste0(base_path, plot_type)
  
  # Initialize a counter
  counter <- 0
  
  # Loop to find a unique filename
  repeat {
    # Generate the candidate filename
    candidate <- if (counter == 0) {
      paste0(base_name, ".", extension)
    } else {
      paste0(base_name, "_", counter, ".", extension)
    }
    
    # Check if the file exists
    if (!file.exists(candidate)) {
      return(list(filename = candidate, run_id = counter)) # Return the filename and run ID
    }
    
    # Increment the counter
    counter <- counter + 1
  }
}