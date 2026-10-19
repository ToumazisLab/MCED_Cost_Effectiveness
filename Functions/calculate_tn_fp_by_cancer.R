calculate_tn_fp_by_cancer <- function(cancer_types, input_matrix, screening_data, screening_cancers, female_cancers, male_cancers, arm = "usual", tn_skip_years = tn_skip_years) {
  
  # Initialize data frames for TN, FP, and Unscreened with rows as cycles and columns as cancer types
  tn_df <- data.frame(matrix(0, nrow = nrow(input_matrix), ncol = length(cancer_types)))
  colnames(tn_df) <- cancer_types
  fp_df <- data.frame(matrix(0, nrow = nrow(input_matrix), ncol = length(cancer_types)))
  colnames(fp_df) <- cancer_types
  unscreened_df <- data.frame(matrix(0, nrow = nrow(input_matrix), ncol = length(cancer_types)))
  colnames(unscreened_df) <- cancer_types
  
  # List to store detailed breakdowns for each cancer type
  detailed_breakdown <- data.frame()
  
  # Filter screening_data to include only cancers in screening_cancers
  filtered_screening_data <- screening_data[screening_data$Cancer %in% screening_cancers, ]
  
  # Sort screening_cancers by specificity (descending) and then by compliance frequency (ascending)
  filtered_screening_data <- filtered_screening_data[order(-filtered_screening_data$Specificity, filtered_screening_data$Compliance_frequency), ]
  screening_cancers_sorted <- filtered_screening_data$Cancer
  
  # Iterate through each cancer type in cancer_types
  for (cancer in cancer_types) {
    if (cancer %in% colnames(input_matrix)) {
      H_vector <- input_matrix[, cancer]  # Initialize H_vector once for each cancer
      H_vector_original <- H_vector
      
      # Prioritize screening for the specific cancer if it is also a screening cancer
      if (cancer %in% screening_cancers_sorted) {
        screening_cancers_sorted <- c(cancer, setdiff(screening_cancers_sorted, cancer))
      }
      
      # Iterate over each cancer in the sorted list of screening_cancers
      for (screening_cancer in screening_cancers_sorted) {
        print(screening_cancer)
        # Get screening info for the current screening cancer
        screening_info <- filtered_screening_data[filtered_screening_data$Cancer == screening_cancer, ]
        
        # Prioritize when cancer == screening_cancer
        if (cancer == screening_cancer) {
          eligibility <- 1
          compliance_freq <- 1  # MCED FPs are fully compliant for their own cancer
        } else if (cancer %in% female_cancers && screening_cancer %in% female_cancers) {
          eligibility <- 1
          compliance_freq <- screening_info$Compliance_frequency
        } else if (cancer %in% male_cancers && screening_cancer %in% female_cancers) {
          eligibility <- 0
          compliance_freq <- 0
        } else if (cancer %in% female_cancers && screening_cancer %in% male_cancers) {
          eligibility <- 0
          compliance_freq <- 0
        } else {
          eligibility <- screening_info$Eligibility
          compliance_freq <- screening_info$Compliance_frequency
        }
        
        if (arm == "mced_tn") {
          tn_skip <- tn_skip_years + 1
          compliance_freq <- compliance_freq / tn_skip
        }
        
        # Apply screening rates to calculate TN and FP for each cycle
        Screened <- H_vector * eligibility * compliance_freq
        Unscreened <- H_vector - Screened
        TN_vector <- Screened * screening_info$Specificity
        FP_vector <- Screened * (1 - screening_info$Specificity)
        
        # Update the H_vector to exclude the screened individuals
        # H_vector <- Unscreened
        
        
        # Assign calculated TN and FP to the respective data frames
        tn_df[, screening_cancer] <- tn_df[, screening_cancer] + TN_vector
        fp_df[, screening_cancer] <- fp_df[, screening_cancer] + FP_vector
        
        # Store detailed breakdown with cycle information
        detailed_breakdown <- rbind(
          detailed_breakdown,
          data.frame(
            Cycle = 1:nrow(input_matrix),
            MCED_Cancer = cancer,
            Screening_Cancer = screening_cancer,
            Outcome = "FP",
            Count = FP_vector
          ),
          data.frame(
            Cycle = 1:nrow(input_matrix),
            MCED_Cancer = cancer,
            Screening_Cancer = screening_cancer,
            Outcome = "TN",
            Count = TN_vector
          ),
          data.frame(
            Cycle = 1:nrow(input_matrix),
            MCED_Cancer = cancer,
            Screening_Cancer = screening_cancer,
            Outcome = "Unscreened",
            Count = Unscreened
          )
        )
      }
      
      # After processing all screening_cancers, store the final Unscreened value in unscreened_df
      unscreened_df[, cancer] <- H_vector
    }
  }
  
  # Return a list containing the detailed breakdown data frame and the aggregated TN, FP, and Unscreened data frames in wide format
  list(
    Detailed_Breakdown = detailed_breakdown,
    Aggregated = list(TN_df = tn_df, FP_df = fp_df, Unscreened_df = unscreened_df)
  )
}
