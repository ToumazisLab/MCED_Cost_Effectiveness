calculate_overdiagnosis <- function(m_count_TP,
                                    incident_matrix,
                                    sensitivity_data,
                                    overdiagnosis_factor,
                                    mced_compliance_rate,
                                    cancer_list,
                                    all_cancer_list,
                                    mced_start_cycle,
                                    mced_end_cycle,
                                    final_result_list,
                                    df,
                                    two_step,
                                    overall_sens) {
  #----------------------------------------------------------------------
  # 1) Pre-calculate the proportions by age overall (Cancer mix)
  #    and by age-stage (Stage distribution for each Cancer).
  #----------------------------------------------------------------------
  df_overall <- calc_proportions_by_age_overall(df, all_cancer_list)
  #   Columns: Cancer, Age, proportion (fraction of all cancers at that age)
  
  df_stage   <- calc_proportions_by_age_by_stage(df, all_cancer_list)
  #   Columns: Cancer, Age, Stage, Count, total_for_cancer_age, prop_by_stage
  
  #----------------------------------------------------------------------
  # 2) Create a matrix to store the overdiagnosed cases (same size as m_count_TP).
  #----------------------------------------------------------------------
  overdiagnosed_matrix <- matrix(0, nrow = nrow(m_count_TP), ncol = ncol(m_count_TP))
  colnames(overdiagnosed_matrix) <- colnames(m_count_TP)
  
  #----------------------------------------------------------------------
  # 3) For each cycle in [mced_start_cycle : mced_end_cycle], 
  #    compute overdiagnosed_cases_this_cycle, then allocate to
  #    each Cancer & Stage according to the proportions & sensitivity.
  #----------------------------------------------------------------------
  for (cycle in mced_start_cycle:mced_end_cycle) {
    # Total potential overdiagnosed cases this cycle
    overdiagnosed_cases_this_cycle <- incident_matrix[cycle, "D"] * 
      mced_compliance_rate * overdiagnosis_factor
    
    if (overdiagnosed_cases_this_cycle <= 0) {
      # If there's nothing to overdiagnose this cycle, skip
      next
    }
    
    # Map cycle to age (assuming cycle = 1 => Age=50, etc.)
    current_age <- cycle + 39
    
    # Get the subset from df_overall for this age
    df_overall_age <- df_overall[df_overall$Age == current_age, ]
    # This tells us, at this age, what fraction each cancer is of the total.
    
    #--------------------------------------------------------------------
    #  For each cancer in df_overall_age, get the fraction_of_cancer
    #  Then for each stage of that cancer (from df_stage), 
    #  multiply fraction_of_stage, then multiply by sensitivity.
    #--------------------------------------------------------------------
    for (idx_cancer in seq_len(nrow(df_overall_age))) {
      this_cancer <- df_overall_age$Cancer[idx_cancer]
      frac_cancer <- df_overall_age$proportion[idx_cancer]
      
      # If fraction is zero or NA, skip
      if (is.na(frac_cancer) || frac_cancer <= 0) next
      
      # Grab the stage distribution for this cancer at this age
      df_stage_subset <- df_stage[df_stage$Cancer == this_cancer &
                                    df_stage$Age == current_age, ]
      
      # For each stage row in this subset
      for (idx_stage in seq_len(nrow(df_stage_subset))) {
        this_stage <- df_stage_subset$Stage[idx_stage]
        frac_stage <- df_stage_subset$prop_by_stage[idx_stage]
        
        # If stage fraction is zero or NA, skip
        if (is.na(frac_stage) || frac_stage <= 0) next
        
        # Potential overdiagnosed for this cancer & stage (before sensitivity)
        potential_overdiagnosed <- overdiagnosed_cases_this_cycle * frac_cancer * frac_stage
        
        # Find the stage-specific sensitivity from sensitivity_data
        # We'll filter by the cancer + stage
        sens_row <- sensitivity_data[sensitivity_data$Cancer == this_cancer &
                                       sensitivity_data$Stage == this_stage, ]
        if (nrow(sens_row) == 1) {
          # There is exactly one row for this Cancer + Stage
          stage_sensitivity <- sens_row$Sensitivity
        } else {
          # If not found or multiple rows, assume 0 or handle differently
          stage_sensitivity <- 0
        }
        
        if (two_step) {
          actual_overdiagnosed <- potential_overdiagnosed * overall_sens * stage_sensitivity
        } else {
          actual_overdiagnosed <- potential_overdiagnosed * stage_sensitivity
        }
        # Actual overdiagnosed cases, after applying sensitivity
        
        
        # We need to place this in the appropriate column in overdiagnosed_matrix
        stage_column <- paste0(this_stage, "_", this_cancer)  # e.g. "I_Lung"
        
        # Only update if that column actually exists in m_count_TP
        if (stage_column %in% colnames(overdiagnosed_matrix)) {
          overdiagnosed_matrix[cycle, stage_column] <- overdiagnosed_matrix[cycle, stage_column] +
            actual_overdiagnosed
          # Also update m_count_TP
          m_count_TP[cycle, stage_column] <- m_count_TP[cycle, stage_column] +
            actual_overdiagnosed
        }
      } # end of stages
    } # end of cancers
  } # end of cycle
  
  #----------------------------------------------------------------------
  # 4) Update final_result_list for only the cancers in `cancer_list`.
  #----------------------------------------------------------------------
  #    We'll assume you have the same structure as before,
  #    with final_result_list[["MCED"]][[cancer]]$MCED_Cancer_cases, etc.
  #----------------------------------------------------------------------
  for (cancer in cancer_list) {
    if (cancer %in% names(final_result_list[["MCED"]])) {
      cancer_data <- final_result_list[["MCED"]][[cancer]]$MCED_Cancer_cases
      stages <- c("I", "II", "III", "IV")
      
      for (stage in stages) {
        column_name <- paste0(stage, "_", cancer)
        overdiagnosed_col <- paste0("overdiagnosed_", stage)
        tp_overdiagnosed_col <- paste0("TP_o_", stage)
        
        if (column_name %in% colnames(overdiagnosed_matrix)) {
          # Set or update overdiagnosed counts for that stage
          cancer_data[[overdiagnosed_col]] <- overdiagnosed_matrix[, column_name]
        } else {
          # If the column doesn't exist in overdiagnosed_matrix, fill with 0
          cancer_data[[overdiagnosed_col]] <- 0
        }
        
        # Now the total TP that includes overdiagnosis
        # (assuming you previously store the normal TP counts in "TP_<stage>")
        if (!paste0("TP_", stage) %in% colnames(cancer_data)) {
          cancer_data[[paste0("TP_", stage)]] <- 0
        }
        
        cancer_data[[tp_overdiagnosed_col]] <- cancer_data[[paste0("TP_", stage)]] +
          cancer_data[[overdiagnosed_col]]
      }
      
      # Summaries
      cancer_data$Total_overdiagnosed <- rowSums(cancer_data[, paste0("overdiagnosed_", stages)], na.rm = TRUE)
      cancer_data$Total_TP_o         <- rowSums(cancer_data[, paste0("TP_o_", stages)], na.rm = TRUE)
      
      # Write back
      final_result_list[["MCED"]][[cancer]]$MCED_Cancer_cases <- cancer_data
    }
  }
  
  #----------------------------------------------------------------------
  # 5) Update the overall summary for cancers in `cancer_list`
  #    (If that is still desired for the summary. This is the same logic as before.)
  #----------------------------------------------------------------------
  summary_df <- final_result_list[["Overall"]][["MCED_Cancer_cases_screening_summary"]]
  stages <- c("I", "II", "III", "IV")
  
  for (stage in stages) {
    overdiagnosed_col <- paste0("overdiagnosed_", stage)
    tp_overdiagnosed_col <- paste0("TP_o_", stage)
    # get all columns that match this stage (e.g. "I_...")
    cancer_columns_stage <- grep(paste0("^", stage, "_"), colnames(overdiagnosed_matrix), value = TRUE)
    
    # sum across the relevant columns *only for the included cancers in cancer_list* 
    # If you strictly want only cancer_list, you can also subset those columns
    # by matching cancer_list. Below is the simpler approach summing all.
    summary_df[[overdiagnosed_col]] <- rowSums(overdiagnosed_matrix[, cancer_columns_stage, drop = FALSE])
    
    # If your summary_df has a "TP_<stage>", we just add it:
    summary_df[[tp_overdiagnosed_col]] <- summary_df[[paste0("TP_", stage)]] + summary_df[[overdiagnosed_col]]
  }
  
  summary_df$Total_overdiagnosed <- rowSums(overdiagnosed_matrix)
  summary_df$Total_TP_o         <- summary_df$TP + summary_df$Total_overdiagnosed
  
  final_result_list[["Overall"]][["MCED_Cancer_cases_screening_summary"]] <- summary_df
  
  #----------------------------------------------------------------------
  # 6) Return outputs
  #----------------------------------------------------------------------
  return(list(
    final_result_list   = final_result_list,
    overdiagnosed_matrix = overdiagnosed_matrix,
    m_count_TP          = m_count_TP
  ))
}



calc_proportions_by_age_overall <- function(df, all_cancer_list) {
  df %>%
    # 1. Keep only the cancers we're interested in
    filter(Cancer %in% all_cancer_list) %>%
    
    # 2. Sum counts across all Stage values for each Cancer-Age combination
    group_by(Cancer, Age) %>%
    summarise(total_for_cancer_age = sum(Rate, na.rm = TRUE), .groups = "drop") %>%
    
    # 3. For each Age, compute total counts (across all cancers)
    group_by(Age) %>%
    mutate(total_for_age = sum(total_for_cancer_age, na.rm = TRUE)) %>%
    
    # 4. Compute the proportion
    mutate(proportion = total_for_cancer_age / total_for_age) %>%
    ungroup() %>%
    
    # 5. Select final columns (no Stage in output)
    select(Cancer, Age, proportion)
}

calc_proportions_by_age_by_stage <- function(df, all_cancer_list) {
  df %>%
    # 1. Keep only the cancers we're interested in
    dplyr::filter(Cancer %in% all_cancer_list) %>%
    
    # 2. For each Cancer-Age group, compute the total count (across all stages)
    dplyr::group_by(Cancer, Age) %>%
    dplyr::mutate(total_for_cancer_age = sum(Rate, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    
    # 3. Calculate the proportion of Count among total_for_cancer_age
    dplyr::mutate(prop_by_stage = Rate / total_for_cancer_age)
}


