time_shift_matrix <- function(data, cancer_type, scenario) {
  # Filter data for the specified cancer type
  filtered_data <- data[data$Cancer == cancer_type, ]
  # Order the data by stage
  filtered_data <- filtered_data[order(filtered_data$Stage), ]
  # Extract the mean dwell times for the specified scenario
  mean_dwell_times <- filtered_data[[scenario]]
  # Number of stages
  n_stages <- length(mean_dwell_times)
  # Initialize the matrix
  matrix <- matrix(NA, nrow = n_stages, ncol = n_stages)
  rownames(matrix) <- filtered_data$Stage
  colnames(matrix) <- filtered_data$Stage
  
  # Populate the matrix
  for (i in 1:n_stages) {
    for (j in i:n_stages) {
      matrix[i, j] <- sum(mean_dwell_times[i:j])
    }
  }
  
  # Setting lower triangle of the matrix to 0
  matrix[lower.tri(matrix)] <- 0
  
  return(matrix)
}

generate_ts_matrix <- function(cancer_list, dwell_data = data_list$dwell_time) {
  unique_scenarios <- c("Scenario_1", "Scenario_2", "Scenario_3", "Scenario_4")
  
  # Initialize a list for time shift matrices
  time_shift_matrices <- setNames(vector("list", length(cancer_list)), cancer_list)
  
  # Populate the list with matrices for each scenario and cancer type
  for (cancer_type in cancer_list) {
    time_shift_matrices[[cancer_type]] <- list()
    
    for (scenario in unique_scenarios) {
      matrix <- time_shift_matrix(dwell_data, cancer_type, scenario)
      time_shift_matrices[[cancer_type]][[scenario]] <- matrix
    }
    
    # Naming the scenario lists for each cancer type
    names(time_shift_matrices[[cancer_type]]) <- paste("time_shift_matrix", seq_along(unique_scenarios), sep = "_")
  }
  return(time_shift_matrices)
}