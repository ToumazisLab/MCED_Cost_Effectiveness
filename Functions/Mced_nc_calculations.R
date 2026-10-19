apply_mced_nc_calculations <- function(final_result_list, screening_data, discount_factor, mced_compliance_rate) {
  
  discount_fn <- function(cycle) {
    1 / ((1 + discount_factor)^(cycle - 1))
  }
  
  usual_cancers <- names(final_result_list[["usual"]])
  n_cycles <- 61
  
  # Initialize overall aggregates for all cancers
  overall_nc_FP <- numeric(n_cycles)
  overall_nc_TN <- numeric(n_cycles)
  overall_nc_compliant <- numeric(n_cycles)
  overall_nc_Invasive <- numeric(n_cycles)
  
  overall_nc_FP_imaging_cost <- numeric(n_cycles)
  overall_nc_FP_diagnostic_cost <- numeric(n_cycles)
  overall_nc_FP_invasive_cost <- numeric(n_cycles)
  overall_nc_TN_imaging_cost <- numeric(n_cycles)
  
  for (cancer in usual_cancers) {
    screening_df <- final_result_list[["usual"]][[cancer]][["FP_TN_Screening_Summary"]]
    if (is.null(screening_df)) next
    
    # Ensure columns exist
    if (!"FP" %in% colnames(screening_df)) screening_df$FP <- 0
    if (!"TN" %in% colnames(screening_df)) screening_df$TN <- 0
    if (!"Invasive" %in% colnames(screening_df)) screening_df$Invasive <- 0
    
    FP <- screening_df$FP
    TN <- screening_df$TN
    Invasive <- screening_df$Invasive
    
    scr_row <- screening_data[screening_data$Cancer == cancer, , drop = FALSE]
    if (nrow(scr_row) == 0) {
      Screening_cost <- 0
      FP_cost <- 0
      invasive_cost <- 0
    } else {
      Screening_cost <- scr_row$Screening_cost
      FP_cost <- scr_row$FP_cost
      invasive_cost <- scr_row$invasive_cost
    }
    
    cycles <- screening_df$Cycle
    discount_vec <- sapply(cycles, discount_fn)
    
    # Compute total costs
    FP_imaging_cost <- FP * Screening_cost * discount_vec
    TN_imaging_cost <- TN * Screening_cost * discount_vec
    FP_diagnostic_cost <- FP * FP_cost * discount_vec
    FP_Invasive_cost <- Invasive * invasive_cost * discount_vec
    
    # Non-compliant fractions
    nc_FP <- FP * (1 - mced_compliance_rate)
    nc_TN <- TN * (1 - mced_compliance_rate)
    nc_Invasive <- Invasive * (1 - mced_compliance_rate)
    
    nc_FP_imaging_cost <- FP_imaging_cost * (1 - mced_compliance_rate)
    nc_FP_diagnostic_cost <- FP_diagnostic_cost * (1 - mced_compliance_rate)
    nc_FP_invasive_cost <- FP_Invasive_cost * (1 - mced_compliance_rate)
    nc_TN_imaging_cost <- TN_imaging_cost * (1 - mced_compliance_rate)
    
    # Summations for per-cancer NC costs
    nc_FP_cost_sum <- nc_FP_imaging_cost + nc_FP_diagnostic_cost + nc_FP_invasive_cost
    nc_TN_cost_sum <- nc_TN_imaging_cost
    nc_FP_TN_cost_sum <- nc_FP_cost_sum + nc_TN_cost_sum
    
    # Update per-cancer MCED data
    mced_tn_fp_df <- final_result_list[["MCED"]][[cancer]][["mced_TN_FP"]]
    # Ensure mced_tn_fp_df has correct length
    if (nrow(mced_tn_fp_df) != length(nc_FP)) {
      stop("Dimension mismatch between mced_TN_FP and screening_df for cancer: ", cancer)
    }
    
    mced_tn_fp_df$nc_FP <- nc_FP
    mced_tn_fp_df$nc_TN <- nc_TN
    mced_tn_fp_df$nc_compliant <- nc_FP + nc_TN
    mced_tn_fp_df$nc_Invasive <- nc_Invasive
    
    final_result_list[["MCED"]][[cancer]][["mced_TN_FP"]] <- mced_tn_fp_df
    
    # Update per-cancer MCED costs
    mced_tn_fp_cost_df <- final_result_list[["MCED"]][[cancer]][["mced_tn_fp_cost"]]
    # Ensure mced_tn_fp_cost_df has correct length
    if (nrow(mced_tn_fp_cost_df) != length(nc_FP)) {
      stop("Dimension mismatch between mced_tn_fp_cost and screening_df for cancer: ", cancer)
    }
    
    mced_tn_fp_cost_df$nc_FP_imaging_cost <- nc_FP_imaging_cost
    mced_tn_fp_cost_df$nc_FP_diagnostic_cost <- nc_FP_diagnostic_cost
    mced_tn_fp_cost_df$nc_FP_invasive_cost <- nc_FP_invasive_cost
    mced_tn_fp_cost_df$nc_TN_imaging_cost <- nc_TN_imaging_cost
    
    mced_tn_fp_cost_df$nc_FP_cost_sum <- nc_FP_cost_sum
    mced_tn_fp_cost_df$nc_TN_cost_sum <- nc_TN_cost_sum
    mced_tn_fp_cost_df$nc_FP_TN_cost_sum <- nc_FP_TN_cost_sum
    
    final_result_list[["MCED"]][[cancer]][["mced_tn_fp_cost"]] <- mced_tn_fp_cost_df
    
    # Update overall sums across all cancers
    overall_nc_FP <- overall_nc_FP + nc_FP
    overall_nc_TN <- overall_nc_TN + nc_TN
    overall_nc_Invasive <- overall_nc_Invasive + nc_Invasive
    
    overall_nc_FP_imaging_cost <- overall_nc_FP_imaging_cost + nc_FP_imaging_cost
    overall_nc_FP_diagnostic_cost <- overall_nc_FP_diagnostic_cost + nc_FP_diagnostic_cost
    overall_nc_FP_invasive_cost <- overall_nc_FP_invasive_cost + nc_FP_invasive_cost
    overall_nc_TN_imaging_cost <- overall_nc_TN_imaging_cost + nc_TN_imaging_cost
  }
  
  # After processing all cancers, compute overall sums
  overall_nc_FP_cost_sum <- overall_nc_FP_imaging_cost + overall_nc_FP_diagnostic_cost + overall_nc_FP_invasive_cost
  overall_nc_TN_cost_sum <- overall_nc_TN_imaging_cost
  overall_nc_FP_TN_cost_sum <- overall_nc_FP_cost_sum + overall_nc_TN_cost_sum
  
  # Update overall MCED data
  overall_mced_TN_FP <- final_result_list[["Overall"]][["mced_TN_FP"]]
  overall_mced_TN_FP$nc_FP <- overall_nc_FP
  overall_mced_TN_FP$nc_TN <- overall_nc_TN
  overall_mced_TN_FP$nc_compliant <- overall_nc_FP + overall_nc_TN
  overall_mced_TN_FP$nc_Invasive <- overall_nc_Invasive
  final_result_list[["Overall"]][["mced_TN_FP"]] <- overall_mced_TN_FP
  
  # Update overall MCED costs
  overall_mced_TN_FP_cost <- final_result_list[["Overall"]][["mced_TN_FP_cost"]]
  overall_mced_TN_FP_cost$nc_FP_imaging_cost <- overall_nc_FP_imaging_cost
  overall_mced_TN_FP_cost$nc_FP_diagnostic_cost <- overall_nc_FP_diagnostic_cost
  overall_mced_TN_FP_cost$nc_FP_invasive_cost <- overall_nc_FP_invasive_cost
  overall_mced_TN_FP_cost$nc_TN_imaging_cost <- overall_nc_TN_imaging_cost
  
  overall_mced_TN_FP_cost$nc_FP_cost_sum <- overall_nc_FP_cost_sum
  overall_mced_TN_FP_cost$nc_TN_cost_sum <- overall_nc_TN_cost_sum
  overall_mced_TN_FP_cost$nc_FP_TN_cost_sum <- overall_nc_FP_TN_cost_sum
  
  final_result_list[["Overall"]][["mced_TN_FP_cost"]] <- overall_mced_TN_FP_cost
  
  return(final_result_list)
}
