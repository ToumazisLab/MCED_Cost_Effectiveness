calculate_deaths_by_cancer_stage <- function(matrix, transposed_survival) {
  # Get the dimensions of the incident matrix
  num_cycles <- nrow(matrix)
  num_cancers <- ncol(matrix)
  
  # Initialize a matrix to store the number of deaths by year and cancer stage
  deaths_matrix <- matrix(0, nrow = num_cycles * 1.5, ncol = num_cancers)
  colnames(deaths_matrix) <- colnames(matrix)
  
  # Filter transposed_survival to include only relevant columns
  common_columns <- intersect(colnames(matrix), colnames(transposed_survival))
  transposed_survival_filtered <- transposed_survival[, common_columns, drop = FALSE]
  
  # Loop through the incident matrix and calculate the deaths
  for (i in 1:num_cycles) {
    for (j in 1:num_cancers) {
      cancer_stage <- colnames(matrix)[j]
      cases <- matrix[i, cancer_stage]
      survival_years <- transposed_survival_filtered[i, cancer_stage]
      
      # Calculate the cycle in which these cases will die
      if (!is.na(survival_years) && survival_years > 0 && cases > 0) {
        death_cycle <- as.numeric(i + floor(survival_years))
        
        # Ensure the death cycle is within the matrix bounds
        if (death_cycle <= nrow(deaths_matrix)) {
          deaths_matrix[death_cycle, cancer_stage] <- deaths_matrix[death_cycle, cancer_stage] + cases
        }
      }
    }
  }
  
  return(deaths_matrix)
}
