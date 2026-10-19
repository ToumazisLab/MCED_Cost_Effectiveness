shift_forward <- function(fn_matrix, dwell_time_df, cancer_types, scenario,
                          fn_shift_cancers, fn_lag_years, final_result_list,
                          mced_start_cycle, mced_end_cycle, screening_data,
                          screening_cancers) {
  # Initialize matrices
  shifted_fn_matrix <- matrix(0, nrow = nrow(fn_matrix), ncol = ncol(fn_matrix))
  colnames(shifted_fn_matrix) <- colnames(fn_matrix)
  
  updated_fn_matrix <- fn_matrix
  
  # New matrices for s_e_s, s_e_i, s_ne_i after final distribution
  fn_s_e_s_matrix <- matrix(0, nrow = nrow(fn_matrix), ncol = ncol(fn_matrix))
  fn_s_e_i_matrix <- matrix(0, nrow = nrow(fn_matrix), ncol = ncol(fn_matrix))
  fn_s_ne_i_matrix <- matrix(0, nrow = nrow(fn_matrix), ncol = ncol(fn_matrix))
  
  colnames(fn_s_e_s_matrix) <- colnames(fn_matrix)
  colnames(fn_s_e_i_matrix) <- colnames(fn_matrix)
  colnames(fn_s_ne_i_matrix) <- colnames(fn_matrix)
  
  stages_order <- c("I", "II", "III", "IV")
  
  for (cancer in cancer_types) {
    # Note: We do NOT skip if !(cancer %in% fn_shift_cancers). 
    # Instead, we handle the classification and record s_ne_i if needed.
    
    all_stages_for_this_cancer <- dwell_time_df$Stage[dwell_time_df$Cancer == cancer]
    stages <- intersect(stages_order, all_stages_for_this_cancer)  # keeps correct order
    
    for (stage_index in seq_along(stages)) {
      current_stage <- stages[stage_index]
      current_stage_col <- paste0(current_stage, "_", cancer)
      
      for (cycle in mced_start_cycle:mced_end_cycle) {
        total_cases <- fn_matrix[cycle, current_stage_col]
        
        # Classify the cases
        result <- classify_cases(total_cases, cancer, cycle, screening_data, screening_cancers)
        
        # Distribute s_ne_i and s_e_s immediately, they do not shift
        fn_s_ne_i_matrix[cycle, current_stage_col] <- fn_s_ne_i_matrix[cycle, current_stage_col] + result$s_ne_i
        fn_s_e_s_matrix[cycle, current_stage_col]  <- fn_s_e_s_matrix[cycle, current_stage_col] + result$s_e_s
        
        # s_e_i are shift candidates only if cancer is in fn_shift_cancers
        # and we have fn_lag_years > 0. Otherwise, no shift occurs and they stay put.
        
        shift_cases <- result$s_e_i
        remain_cases <- total_cases - shift_cases
        
        # Update remain_cases in updated_fn_matrix
        updated_fn_matrix[cycle, current_stage_col] <- remain_cases
        
        # Check if this cancer should actually shift cases
        if ((cancer %in% fn_shift_cancers) && (shift_cases > 0) && (fn_lag_years > 0)) {
          out <- forward_fn_shift(shift_cases, cycle, current_stage, stages, stage_index, 
                                  dwell_time_df, scenario, cancer, fn_lag_years, 
                                  dim(fn_matrix), updated_fn_matrix, shifted_fn_matrix)
          updated_fn_matrix <- out$updated_fn_matrix
          shifted_fn_matrix <- out$shifted_fn_matrix
          
          if (out$actually_shifted) {
            # s_e_i moved to a new position
            final_stage_col <- paste0(out$final_stage, "_", cancer)
            fn_s_e_i_matrix[out$shift_cycle, final_stage_col] <- fn_s_e_i_matrix[out$shift_cycle, final_stage_col] + shift_cases
          } else {
            # No actual shift in stage or cycle, s_e_i stays at original position
            fn_s_e_i_matrix[cycle, current_stage_col] <- fn_s_e_i_matrix[cycle, current_stage_col] + shift_cases
            
            # Also add these cases back since no shift occurred
            updated_fn_matrix[cycle, current_stage_col] <- updated_fn_matrix[cycle, current_stage_col] + shift_cases
          }
        } else {
          # No shifting occurs
          fn_s_e_i_matrix[cycle, current_stage_col] <- fn_s_e_i_matrix[cycle, current_stage_col] + shift_cases
          
          # If we are not shifting at all, add them back to updated_fn_matrix
          updated_fn_matrix[cycle, current_stage_col] <- updated_fn_matrix[cycle, current_stage_col] + shift_cases
        }
      }
    }
  }
  
  # Create unshifted_fn_matrix by subtracting shifted cases from original
  unshifted_fn_matrix <- fn_matrix - shifted_fn_matrix
  
  # Update final_result_list for each cancer type
  for (cancer in cancer_types) {
    if (cancer %in% names(final_result_list[["MCED"]])) {
      cancer_data <- final_result_list[["MCED"]][[cancer]]$MCED_Cancer_cases
      
      # For each stage, create the *_FN_shifted_only and *_FN_shifted columns
      for (stage in stages_order) {
        col_name <- paste0(stage, "_", cancer)
        fn_shifted_only_col <- paste0(stage, "_FN_shifted_only")
        fn_shifted_col <- paste0(stage, "_FN_shifted")
        
        # Stage_FN_shifted_only from shifted_fn_matrix
        cancer_data[[fn_shifted_only_col]] <- shifted_fn_matrix[, col_name]
        
        # Stage_FN_shifted from updated_fn_matrix (final distribution)
        cancer_data[[fn_shifted_col]] <- updated_fn_matrix[, col_name]
      }
      
      # Total_FN_missed = sum of all stage_FN_shifted_only for this cancer
      stage_shifted_only_cols <- paste0(stages_order, "_FN_shifted_only")
      cancer_data$Total_FN_missed <- rowSums(cancer_data[, stage_shifted_only_cols, drop = FALSE], na.rm = TRUE)
      
      # Total_FN_shifted = sum of all stage_FN_shifted for this cancer
      stage_shifted_cols <- paste0(stages_order, "_FN_shifted")
      cancer_data$Total_FN_shifted <- rowSums(cancer_data[, stage_shifted_cols, drop = FALSE], na.rm = TRUE)
      
      final_result_list[["MCED"]][[cancer]]$MCED_Cancer_cases <- cancer_data
    }
  }
  
  # Update the overall summary
  summary_df <- final_result_list[["Overall"]][["MCED_Cancer_cases_screening_summary"]]
  
  for (stage in stages_order) {
    fn_shifted_only_col <- paste0("Total_FN_shifted_only_", stage)
    fn_shifted_col <- paste0("Total_FN_shifted_", stage)
    
    # Identify all columns for this stage across all cancers
    cancer_columns <- grep(paste0("^", stage, "_"), colnames(shifted_fn_matrix), value = TRUE)
    
    # Sum across all cancers for shifted_only
    summary_df[[fn_shifted_only_col]] <- rowSums(shifted_fn_matrix[, cancer_columns, drop = FALSE], na.rm = TRUE)
    
    # Sum across all cancers for shifted
    summary_df[[fn_shifted_col]] <- rowSums(updated_fn_matrix[, cancer_columns, drop = FALSE], na.rm = TRUE)
  }
  
  # Total_FN_missed = sum of all shifted_only (all stages, all cancers)
  total_shifted_only_cols <- paste0("Total_FN_shifted_only_", stages_order)
  summary_df$Total_FN_missed <- rowSums(summary_df[, total_shifted_only_cols, drop = FALSE], na.rm = TRUE)
  
  # Total_FN_shifted = sum of all shifted (updated) (all stages, all cancers)
  total_shifted_cols <- paste0("Total_FN_shifted_", stages_order)
  summary_df$Total_FN_shifted <- rowSums(summary_df[, total_shifted_cols, drop = FALSE], na.rm = TRUE)
  
  final_result_list[["Overall"]][["MCED_Cancer_cases_screening_summary"]] <- summary_df
  
  updated_fn_matrix <- fn_s_e_s_matrix + fn_s_e_i_matrix + fn_s_ne_i_matrix
  return(list(
    final_result_list = final_result_list,
    shifted_fn_matrix = shifted_fn_matrix,
    unshifted_fn_matrix = unshifted_fn_matrix,
    updated_fn_matrix = updated_fn_matrix,
    fn_s_e_s_matrix = fn_s_e_s_matrix,
    fn_s_e_i_matrix = fn_s_e_i_matrix,
    fn_s_ne_i_matrix = fn_s_ne_i_matrix
  ))
}
