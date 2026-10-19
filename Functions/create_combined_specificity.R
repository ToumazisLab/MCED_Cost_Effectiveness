create_combined_specificity <- function(specificity, cancer_list, min_specificity, max_specificity) {
  # Step 1: Find the missing cancer types
  missing_cancers <- setdiff(cancer_list, specificity$cancer_type)
  
  # Step 2: Get the mean and standard deviation of specificity values from the existing dataframe
  min_specificity <- min(specificity$specificity, na.rm = TRUE)
  max_specificity <- max(specificity$specificity, na.rm = TRUE)
  
  # Step 3: Create a dataframe with the missing cancer types and assign random specificity values
  set.seed(123)  # Optional: Set a seed for reproducibility
  missing_specificity <- data.frame(
    cancer_type = missing_cancers,
    specificity = runif(length(missing_cancers), min = min_specificity, max = max_specificity)
  )
  
  # Step 4: Combine the original specificity dataframe with the newly created one
  combined_specificity <- rbind(specificity, missing_specificity)
  
  # Step 5: Rank the specificity values (lower values get a lower rank)
  combined_specificity <- combined_specificity %>%
    arrange(specificity) %>%
    mutate(rank = rank(specificity, ties.method = "min"))
  
  # Step 6: Return the combined dataframe
  return(combined_specificity)
}