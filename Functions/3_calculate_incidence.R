calculate_incidents <- function(matrix) {
  # Initialize a matrix to store incident cases
  incident_matrix <- matrix(0, nrow = nrow(matrix) - 1, ncol = ncol(matrix))
  colnames(incident_matrix) <- colnames(matrix)
  
  # Handle H, D, and Age columns
  incident_matrix[, "H"] <- matrix[-1, "H"]
  
  # Calculate incident cases for each state (column) except H, D, and Age
  for (col in setdiff(colnames(matrix), c("H"))) {
    for (row in 2:nrow(matrix)) {
      incident_cases = matrix[row, col] - matrix[row - 1, col]
      incident_matrix[row - 1, col] <- max(0, incident_cases) # Ensure no negative values
    }
  }
  
  # Trim the incident matrix to remove columns H, D, Age
  incident_matrix_trimmed <- incident_matrix[, !colnames(incident_matrix) %in% c("H", "D")]
  
  return(list(incident_matrix = incident_matrix, incident_matrix_trimmed = incident_matrix_trimmed))
}