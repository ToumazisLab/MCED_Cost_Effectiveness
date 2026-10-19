generate_m_count_cancers <- function(final_result_list, cancer_list, incident_matrix_trimmed) {
  stages <- c("I", "II", "III", "IV")
  
  # Initialize a list to store columns
  col_data_list <- list()
  
  for (cancer_type in cancer_list) {
    # Extract the data frame for the given cancer
    cancer_data <- final_result_list[["MCED"]][[cancer_type]][["mced_cancer_cases_original"]]
    
    # For each stage, sum the required columns
    for (st in stages) {
      non_compliant_col <- paste0("mced_non_compliant_", st)
      tp_shifted_col <- paste0(st, "_TP_shifted")
      fn_shifted_col <- paste0(st, "_FN_shifted")
      
      # Check existence of columns
      required_cols <- c(non_compliant_col, tp_shifted_col, fn_shifted_col)
      missing_cols <- required_cols[!required_cols %in% colnames(cancer_data)]
      if (length(missing_cols) > 0) {
        stop(paste("Missing required columns for", st, cancer_type, ":", paste(missing_cols, collapse = ", ")))
      }
      
      combined_values <- cancer_data[[non_compliant_col]] + 
        cancer_data[[tp_shifted_col]] + 
        cancer_data[[fn_shifted_col]]
      
      final_col_name <- paste(st, cancer_type, sep = "_")
      col_data_list[[final_col_name]] <- combined_values
    }
  }
  
  # Combine into a matrix
  final_matrix <- do.call(cbind, col_data_list)
  
  # Reorder columns to match incident_matrix_trimmed
  desired_col_order <- colnames(incident_matrix_trimmed)
  
  # Check if all desired columns are present
  missing_in_final <- setdiff(desired_col_order, colnames(final_matrix))
  if (length(missing_in_final) > 0) {
    warning("Some columns in incident_matrix_trimmed are not found in final_matrix: ",
            paste(missing_in_final, collapse = ", "))
  }
  
  # Match order and subset
  common_cols <- intersect(desired_col_order, colnames(final_matrix))
  final_matrix <- final_matrix[, common_cols, drop = FALSE]
  
  # Save under final_result_list[["Overall"]] as MCED_count_cancers_final
  final_result_list[["Overall"]][["MCED_count_cancers_final"]] <- final_matrix
  
  return(final_result_list)
}
