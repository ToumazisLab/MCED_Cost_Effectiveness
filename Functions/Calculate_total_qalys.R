calculate_total_qalys_old <- function(input_matrix, transposed_survival, utility_long, expanded_utility) {
  results <- data.frame(Cancer_Stage = character(),
                        Age = integer(),
                        QALYs = numeric(),
                        stringsAsFactors = FALSE)
  
  for (cancer_stage in colnames(input_matrix)) {
    for (age_index in 1:nrow(input_matrix)) {
      cases <- input_matrix[age_index, cancer_stage]
      survival_years <- as.numeric(transposed_survival[age_index, cancer_stage])
      
      total_survival <- cases * survival_years
      # Calculate full survival years and the fractional part of the last year
      full_years <- floor(survival_years)
      partial_year_fraction <- survival_years - full_years
      
      # Generate age-adjusted utility vector
      max_age_for_utility_vector <- age_index + full_years - 1 + ifelse(partial_year_fraction > 0, 1, 0)
      if(max_age_for_utility_vector > nrow(expanded_utility)) {
        max_age_for_utility_vector <- nrow(expanded_utility)
      }
      age_adjusted_utility_vector <- expanded_utility$Discounted_avg[age_index:min(max_age_for_utility_vector, nrow(expanded_utility))]
      
      # Generate utility multipliers vector
      utility_multipliers_for_stage <- utility_long[utility_long$Cancer_Stage == cancer_stage, "utility", drop = TRUE]
      utility_multipliers <- rep(1, length(age_adjusted_utility_vector)) # Default to 1 for years beyond year 5
      for (year in 1:min(5, full_years + ifelse(partial_year_fraction > 0, 1, 0))) {
        utility_multipliers[year] <- utility_multipliers_for_stage[year]
      }
      
      # Ensure both vectors are of equal length
      if (length(age_adjusted_utility_vector) < length(utility_multipliers)) {
        age_adjusted_utility_vector <- c(age_adjusted_utility_vector, rep(age_adjusted_utility_vector[length(age_adjusted_utility_vector)], length(utility_multipliers) - length(age_adjusted_utility_vector)))
      }
      
      # Adjust the utility and multiplier for the partial year
      if (partial_year_fraction > 0 && length(age_adjusted_utility_vector) > length(utility_multipliers_for_stage)) {
        utility_multipliers[length(utility_multipliers)] <- 1
      }
      if (partial_year_fraction > 0) {
        age_adjusted_utility_vector[length(age_adjusted_utility_vector)] <- age_adjusted_utility_vector[length(age_adjusted_utility_vector)] * partial_year_fraction
        utility_multipliers[length(utility_multipliers)] <- utility_multipliers[length(utility_multipliers)] # * partial_year_fraction
      }
      
      # Calculate QALYs
      qalys_for_age_and_cancer_stage <- sum(age_adjusted_utility_vector * utility_multipliers) * cases
      
      # Append to results
      results <- rbind(results, data.frame(Cancer_Stage = cancer_stage,
                                           Age = age_index + 39,
                                           QALYs = qalys_for_age_and_cancer_stage,
                                           stringsAsFactors = FALSE))
    }
  }
  
  return(results)
}

calculate_total_qalys <- function(input_matrix, transposed_survival, utility_long, expanded_utility) {
  n_rows <- nrow(input_matrix)
  n_cols <- ncol(input_matrix)
  cancer_stages <- colnames(input_matrix)
  
  # Preallocate vectors instead of rbind()
  total_entries <- n_rows * n_cols
  res_stage <- character(total_entries)
  res_age   <- integer(total_entries)
  res_qaly  <- numeric(total_entries)
  counter <- 1
  
  # Preprocess utility multipliers: list of numeric vectors
  utility_lookup <- split(utility_long$utility, utility_long$Cancer_Stage)
  
  for (j in seq_along(cancer_stages)) {
    cancer_stage <- cancer_stages[j]
    utility_multipliers_for_stage <- utility_lookup[[cancer_stage]]
    
    for (i in seq_len(n_rows)) {
      cases <- input_matrix[i, j]
      survival_years <- as.numeric(transposed_survival[i, j])
      
      full_years <- floor(survival_years)
      partial_year <- survival_years - full_years
      
      max_years <- full_years + if (partial_year > 0) 1 else 0
      max_age_index <- min(i + max_years - 1, nrow(expanded_utility))
      
      # Age-adjusted utility vector
      age_util <- expanded_utility$Discounted_avg[i:max_age_index]
      
      # Utility multipliers (up to 5 years, rest default to 1)
      utility_multipliers <- rep(1, length(age_util))
      if (length(utility_multipliers_for_stage) > 0) {
        utility_multipliers[seq_len(min(5, length(utility_multipliers)))] <-
          utility_multipliers_for_stage[seq_len(min(5, length(utility_multipliers_for_stage), length(utility_multipliers)))]
      }
      
      # Handle partial year adjustment
      if (partial_year > 0 && length(age_util) > 0) {
        age_util[length(age_util)] <- age_util[length(age_util)] * partial_year
      }
      
      # QALYs for this cell
      qalys <- sum(age_util * utility_multipliers) * cases
      
      # Store results
      res_stage[counter] <- cancer_stage
      res_age[counter]   <- i + 39
      res_qaly[counter]  <- qalys
      counter <- counter + 1
    }
  }
  
  # Return as data.frame
  return(data.frame(Cancer_Stage = res_stage, Age = res_age, QALYs = res_qaly, stringsAsFactors = FALSE))
}
