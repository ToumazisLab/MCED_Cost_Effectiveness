initialize_cancer_metrics <- function(cancer_list) {
  # Define the arms
  arm_names <- c("MCED", "usual", "Overall")
  
  # Initialize the main list to store each arm
  metrics <- list()
  
  # Loop through each arm
  for (arm in arm_names) {
    # For each arm, create a list of cancers, each initialized to an empty list
    metrics[[arm]] <- setNames(vector("list", length(cancer_list)), cancer_list)
  }
  
  # Return the structured list
  return(metrics)
}
