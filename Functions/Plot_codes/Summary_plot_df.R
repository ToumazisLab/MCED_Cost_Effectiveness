sum_mced_usual_costs <- function(final_result_list) {
  # Initialize a list to store the results
  summary_list <- list()
  
  # MCED Cancer Costs
  m_cancer_cost <- final_result_list[["Overall"]][["MCED_Cancer_costs_summary"]]
  if (!is.null(m_cancer_cost)) {
    for (stage in c("I", "II", "III", "IV")) {
      for (suffix in c("TP", "FN", "NA")) {
        stage_cost_col <- paste0(stage, "_post_cost_", suffix)
        stage_utility_col <- paste0(stage, "_post_utility_", suffix)
        
        summary_list[[paste0("mced_", stage, "_pd_cost_", suffix)]] <- sum(m_cancer_cost[[stage_cost_col]], na.rm = TRUE)
        summary_list[[paste0("mced_", stage, "_pd_utility_", suffix)]] <- sum(m_cancer_cost[[stage_utility_col]], na.rm = TRUE)
      }
    }
    summary_list$mced_cancer_diagnosed_cost <- sum(m_cancer_cost$total_m_dx_per_cycle, na.rm = TRUE)
    summary_list$mced_cancer_screen_diagnosed_cost <- sum(m_cancer_cost$total_s_dx_per_cycle, na.rm = TRUE)
    summary_list$mced_cancer_incidentally_diagnosed_cost <- sum(m_cancer_cost$total_i_dx_per_cycle, na.rm = TRUE)
  }
  
  # MCED Cancer Test Costs
  m_cancer_test_cost <- final_result_list[["Overall"]][["MCED_Cancer_cases_screening_summary"]]
  if (!is.null(m_cancer_test_cost)) {
    summary_list$mced_cancer_test_cost <- sum(m_cancer_test_cost$mced_test_cost, na.rm = TRUE)
  }
  
  # MCED TN/FP Costs
  m_tn_fp_cost <- final_result_list[["Overall"]][["mced_TN_FP_cost"]]
  if (!is.null(m_tn_fp_cost)) {
    summary_list$mced_FP_Imaging_cost <- sum(m_tn_fp_cost$mced_FP_Imaging, na.rm = TRUE)
    summary_list$mced_FP_Diagnostic_cost <- sum(m_tn_fp_cost$mced_FP_Diagnostic, na.rm = TRUE)
    summary_list$mced_FP_Invasive_cost <- sum(m_tn_fp_cost$mced_FP_Invasive, na.rm = TRUE)
    summary_list$mced_TN_Imaging_cost <- sum(m_tn_fp_cost$mced_TN_Imaging, na.rm = TRUE)
    summary_list$mced_TN_Diagnostic_cost <- sum(m_tn_fp_cost$mced_TN_Diagnostic, na.rm = TRUE)
    summary_list$mced_TN_Invasive_cost <- sum(m_tn_fp_cost$mced_TN_Invasive, na.rm = TRUE)
    summary_list$mced_H_nc_cost <- sum(m_tn_fp_cost$nc_FP_TN_cost_sum, na.rm = TRUE)
  }
  
  # MCED TN/FP Utility
  m_tn_fp_utility <- final_result_list[["Overall"]][["mced_TN_FP"]]
  if (!is.null(m_tn_fp_utility)) {
    summary_list$mced_tn_fp_mced_test_cost <- sum(m_tn_fp_utility$mced_test_cost, na.rm = TRUE)
    summary_list$mced_screening_disutility <- sum(m_tn_fp_utility$screening_disutility, na.rm = TRUE)
    summary_list$mced_diagnostic_disutility <- sum(m_tn_fp_utility$diagnosis_disutility, na.rm = TRUE)
    summary_list$mced_invasive_disutility <- sum(m_tn_fp_utility$invasive_disutility, na.rm = TRUE)
    summary_list$mced_m_invasive_disutility <- sum(m_tn_fp_utility$m_invasive_disutility, na.rm = TRUE)
    summary_list$mced_adjusted_utility <- sum(m_tn_fp_utility$Adjusted_population_utility, na.rm = TRUE)
  }
  
  # Usual Care Cancer Costs
  u_cancer_cost <- final_result_list[["Overall"]][["Cancer_costs_summary"]]
  if (!is.null(u_cancer_cost)) {
    summary_list$usual_cancer_screen_diagnosed_cost <- sum(u_cancer_cost$total_s_dx_per_cycle, na.rm = TRUE)
    summary_list$usual_cancer_incidentally_diagnosed_cost <- sum(u_cancer_cost$total_i_dx_per_cycle, na.rm = TRUE)
  }
  
  # Usual Care Cancer Summary
  u_cancer_summary <- final_result_list[["Overall"]][["Cancer_cases_screening_summary"]]
  if (!is.null(u_cancer_summary)) {
    for (stage in c("I", "II", "III", "IV")) {
      stage_cost_col <- paste0(stage, "_post_cost")
      stage_utility_col <- paste0(stage, "_post_utility")
      
      summary_list[[paste0("usual_", stage, "_pd_cost")]] <- sum(u_cancer_summary[[stage_cost_col]], na.rm = TRUE)
      summary_list[[paste0("usual_", stage, "_pd_utility")]] <- sum(u_cancer_summary[[stage_utility_col]], na.rm = TRUE)
    }
  }
  
  # Usual Care TN/FP Summary
  u_tn_fp_summary <- final_result_list[["Overall"]][["FP_TN_Screening_Summary"]]
  if (!is.null(u_tn_fp_summary)) {
    summary_list$usual_FP_Imaging_cost <- sum(u_tn_fp_summary$FP_imaging_cost, na.rm = TRUE)
    summary_list$usual_TN_Imaging_cost <- sum(u_tn_fp_summary$TN_imaging_cost, na.rm = TRUE)
    summary_list$usual_FP_Diagnostic_cost <- sum(u_tn_fp_summary$FP_diagnosis_cost, na.rm = TRUE)
    summary_list$usual_TN_Diagnostic_cost <- sum(u_tn_fp_summary$TN_diagnosis_cost, na.rm = TRUE)
    summary_list$usual_FP_Invasive_cost <- sum(u_tn_fp_summary$Invasive_cost, na.rm = TRUE)
    summary_list$usual_total_diagnostic_disutility <- sum(u_tn_fp_summary$total_fp_disutility, na.rm = TRUE)
    summary_list$usual_total_invasive_disutility <- sum(u_tn_fp_summary$total_invasive_disutility, na.rm = TRUE)
    summary_list$usual_total_tn_disutility <- sum(u_tn_fp_summary$total_tn_disutility, na.rm = TRUE)
    summary_list$usual_adjusted_utility <- sum(u_tn_fp_summary$adjusted_utility, na.rm = TRUE)
  }
  
  # Convert summary_list to a data frame
  summary_df <- as.data.frame(t(as.data.frame(summary_list)))
  colnames(summary_df) <- "Value"
  summary_df$Variable <- rownames(summary_df)
  
  # Add 'Arm' column to differentiate between MCED and Usual Care
  summary_df$Arm <- ifelse(grepl("^mced_", summary_df$Variable), "mced", "usual")
  summary_df$Type <- ifelse(grepl("_cost", summary_df$Variable), "cost", "utility")
  
  # Add 'State' column based on Variable name
  summary_df$State <- ifelse(
    grepl("FP|TN|adjusted|tn_fp|H_nc", summary_df$Variable, ignore.case = FALSE), "Healthy",
    ifelse(grepl("FN|TP|NA|pd_cost|pd_utility|cancer", summary_df$Variable, ignore.case = TRUE), "Cancer", NA)
  )
  
  # Standardize naming conventions
  summary_df$Variable <- gsub("^mced_", "", summary_df$Variable)
  summary_df$Variable <- gsub("^usual_", "", summary_df$Variable)
  summary_df$Variable <- gsub("FP_imaging_cost", "FP_Imaging_Cost", summary_df$Variable)
  summary_df$Variable <- gsub("mced_FP_Imaging", "FP_Imaging_Cost", summary_df$Variable)
  # Repeat for other variables with similar logic
  
  # Restructure data for plotting
  plot_df <- summary_df %>%
    dplyr::group_by(Arm, Type, Variable, State) %>%
    dplyr::summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop")
  
  # Add the modified data frame back to the final_result_list
  final_result_list[["Overall"]][["plot_df"]] <- plot_df
  
  return(final_result_list)
}


calculate_icer <- function(final_result_list) {
  # Extract the plot_df created in the previous function
  plot_df <- final_result_list[["Overall"]][["plot_df"]]
  
  if (is.null(plot_df)) {
    stop("The plot_df is not available in the final_result_list. Ensure the data is prepared correctly.")
  }
  
  # Exclude rows with "disutility" in the Variable name
  plot_df_filtered <- plot_df[!grepl("disutility", plot_df$Variable, ignore.case = TRUE), ]
  
  # Sum costs and utilities by arm using the "Type" column
  summary_by_arm <- plot_df_filtered %>%
    dplyr::group_by(Arm, Type) %>%
    dplyr::summarise(
      Total_Value = sum(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = Type,
      values_from = Total_Value,
      values_fill = list(Total_Value = 0)
    )
  
  # Ensure both arms are present and costs/utilities are available
  if (!all(c("mced", "usual") %in% summary_by_arm$Arm)) {
    stop("Both MCED and Usual Care arms must be present to calculate ICER.")
  }
  if (!all(c("cost", "utility") %in% colnames(summary_by_arm))) {
    stop("The 'Type' column must contain both 'cost' and 'utility' labels.")
  }
  
  # Extract values for each arm
  mced <- summary_by_arm[summary_by_arm$Arm == "mced", ]
  usual <- summary_by_arm[summary_by_arm$Arm == "usual", ]
  
  # Calculate incremental costs and utilities
  incremental_cost <- mced$cost - usual$cost
  incremental_utility <- mced$utility - usual$utility
  
  # Calculate ICER
  icer <- if (incremental_utility == 0) {
    NA  # Avoid division by zero
  } else {
    incremental_cost / incremental_utility
  }
  
  # Create a DataFrame for the ICER results
  icer_results <- data.frame(
    Arm = "summary", # These are summary results
    Type = c("cost", "utility", "ICER"),
    Variable = c("Incremental_Cost", "Incremental_Utility", "ICER"),
    State = NA, # No state for summary rows
    Value = c(incremental_cost, incremental_utility, icer)
  )
  
  # Add arm-specific summaries for MCED and Usual Care
  arm_results <- summary_by_arm %>%
    tidyr::pivot_longer(cols = c("cost", "utility"), names_to = "Type", values_to = "Value") %>%
    dplyr::mutate(
      Variable = paste0(Arm, "_", Type), # e.g., "mced_cost", "usual_utility"
      State = NA
    ) %>%
    dplyr::select(Arm, Type, Variable, State, Value)
  
  # Combine the original plot_df with the new ICER and summary results
  final_icer_results <- dplyr::bind_rows(plot_df, arm_results, icer_results)
  
  # Add the final_icer_results to the final_result_list
  final_result_list[["Overall"]][["final_icer_result"]] <- final_icer_results
  
  return(final_result_list)
}

