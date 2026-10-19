# Function to run the Markov model for each cancer and store the results in a list
create_markov_model_list <- function(cancer_list, tp_matrices, v_names_states, num_people) {
  m_M_SoC_list <- list()
  
  for (cancer in cancer_list) {
    # Create the 3D array a_P_SoC for the current cancer
    a_P_SoC <- array(0, dim = c(6, 6, length(40:100)))
    for (age in 40:100) {
      matrix_name <- paste(age, cancer, sep = ".")
      a_P_SoC[, , age - 39] <- tp_matrices[[matrix_name]]
    }
    
    # Run the Markov model for the current cancer
    m_M_SoC <- run_markov_model_SoC(a_P_SoC, v_names_states, num_people)
    
    # Store the resulting m_M_SoC matrix in the list
    m_M_SoC_list[[cancer]] <- m_M_SoC
  }
  
  return(m_M_SoC_list)
}


# Function to create the final cumulative matrix
create_final_cumulative_matrix <- function(cancer_list, m_M_SoC_list) {
  # Initialize the final cumulative matrix
  final_cumulative_matrix <- data.frame(
    H = numeric(nrow(m_M_SoC_list[[1]])),
    D = numeric(nrow(m_M_SoC_list[[1]]))
  )
  
  # Add columns for each cancer and stage in the final matrix
  for (cancer in cancer_list) {
    for (stage in c("I", "II", "III", "IV")) {
      col_name <- paste(stage, cancer, sep = "_")
      final_cumulative_matrix[[col_name]] <- numeric(nrow(m_M_SoC_list[[1]]))
    }
  }
  
  # Sum the H and D columns across all m_M_SoC matrices
  for (cancer in cancer_list) {
    m_M_SoC <- m_M_SoC_list[[cancer]]
    final_cumulative_matrix$H <- final_cumulative_matrix$H + m_M_SoC$H
    final_cumulative_matrix$D <- final_cumulative_matrix$D + m_M_SoC$D
    
    # Add stage-specific columns to the final cumulative matrix
    for (stage in c("I", "II", "III", "IV")) {
      col_name <- paste(stage, cancer, sep = "_")
      final_cumulative_matrix[[col_name]] <- m_M_SoC[[stage]]
    }
  }
  
  # Average the H and D columns across all cancers
  final_cumulative_matrix$H <- final_cumulative_matrix$H / length(cancer_list)
  final_cumulative_matrix$D <- final_cumulative_matrix$D / length(cancer_list)
  
  return(final_cumulative_matrix)
}


# Function to normalize the final cumulative matrix
normalize_final_cumulative_matrix <- function(final_cumulative_matrix, num_people) {
  # Identify the columns to be normalized (H and D)
  columns_to_normalize <- c("H", "D")
  
  # Calculate the row sums (excluding the 'Age' column if present)
  row_sums <- rowSums(final_cumulative_matrix)
  
  # Normalize only the H and D columns to ensure each row sums to num_people
  normalized_matrix <- final_cumulative_matrix
  for (i in 1:nrow(final_cumulative_matrix)) {
    if (row_sums[i] > 0) {
      # Calculate the total of the other columns
      other_columns_sum <- sum(final_cumulative_matrix[i, setdiff(colnames(final_cumulative_matrix), columns_to_normalize)])
      
      # Adjust the H and D columns proportionally to make the row sum to num_people
      normalized_matrix[i, columns_to_normalize] <- final_cumulative_matrix[i, columns_to_normalize] * ((num_people - other_columns_sum) / sum(final_cumulative_matrix[i, columns_to_normalize]))
    }
  }
  
  return(normalized_matrix)
}


# Wrapper function to run all steps
run_cancer_model_pipeline <- function(cancer_list, incidence_final, all_cause, female_prop, male_prop, female_cancers, v_names_states, num_people) {
  # Step 1: Generate transition probability matrices
  tp_matrices <- generate_tp_matrices(cancer_list, incidence_final, all_cause, female_prop, male_prop, female_cancers)
  
  # Step 2: Create the list of Markov model results for each cancer
  m_M_SoC_list <- create_markov_model_list(cancer_list, tp_matrices, v_names_states, num_people)
  
  # Step 3: Create the final cumulative matrix
  cumulative_matrix <- create_final_cumulative_matrix(cancer_list, m_M_SoC_list)
  
  final_cumulative_matrix <- normalize_final_cumulative_matrix(cumulative_matrix, num_people)
  
  return(final_cumulative_matrix)
}
