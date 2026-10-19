apply_usual_costs <- function(final_result_list, screening_data, discount_factor) {
  
  discount_fn <- function(cycle) {
    1 / ((1 + discount_factor)^(cycle - 1))
  }
  
  usual_cancers <- names(final_result_list[["usual"]])
  n_cycles <- nrow(final_result_list[["Overall"]][["FP_TN_Screening_Summary"]])
  
  # Initialize overall sums for usual care costs
  FP_imaging_cost_sum <- numeric(n_cycles)
  TN_imaging_cost_sum <- numeric(n_cycles)
  FP_diagnostic_cost_sum <- numeric(n_cycles)
  TN_diagnostic_cost_sum <- numeric(n_cycles)
  FP_Invasive_cost_sum <- numeric(n_cycles)
  TN_Invasive_cost_sum <- numeric(n_cycles)
  Imaging_cost_sum <- numeric(n_cycles)
  Diagnostic_cost_sum <- numeric(n_cycles)
  Invasive_cost_sum <- numeric(n_cycles)
  FP_cost_sum <- numeric(n_cycles)
  TN_cost_sum <- numeric(n_cycles)
  FP_TN_cost_sum <- numeric(n_cycles)
  
  # Process each cancer for usual care
  for (cancer in usual_cancers) {
    screening_df <- final_result_list[["usual"]][[cancer]][["FP_TN_Screening_Summary"]]
    if (is.null(screening_df)) next
    
    scr_row <- screening_data[screening_data$Cancer == cancer, , drop = FALSE]
    if (nrow(scr_row) == 0) {
      stop("UC TN/FP: No Screening data found")
    } else {
      Screening_cost <- scr_row$Screening_cost
      FP_cost <- scr_row$FP_cost
      invasive_cost <- scr_row$invasive_cost
    }
    
    # Get cycles and discount
    cycles <- screening_df$Cycle
    discount_vec <- sapply(cycles, discount_fn)
    TN_cost <- 0
    TN_invasive_cost <- 0
    
    # Compute usual care costs by cancer
    FP_imaging_cost <- screening_df$FP * Screening_cost * discount_vec
    TN_imaging_cost <- screening_df$TN * Screening_cost * discount_vec
    FP_diagnostic_cost <- screening_df$FP * FP_cost * discount_vec
    TN_diagnostic_cost <- screening_df$TN * TN_cost * discount_vec
    FP_Invasive_cost <- screening_df$Invasive * invasive_cost * discount_vec
    TN_Invasive_cost <- TN_invasive_cost * discount_vec
    
    Imaging_cost <- FP_imaging_cost + TN_imaging_cost
    Diagnostic_cost <- FP_diagnostic_cost + TN_diagnostic_cost
    Invasive_cost <- FP_Invasive_cost + TN_Invasive_cost
    
    FP_cost <- FP_imaging_cost + FP_diagnostic_cost + FP_Invasive_cost
    TN_cost <- TN_imaging_cost + TN_diagnostic_cost + TN_Invasive_cost
    FP_TN_cost <- FP_cost + TN_cost
    
    # Update the screening_df with these cost columns
    screening_df$FP_imaging_cost <- FP_imaging_cost
    screening_df$TN_imaging_cost <- TN_imaging_cost
    screening_df$FP_diagnostic_cost <- FP_diagnostic_cost
    screening_df$TN_diagnostic_cost <- TN_diagnostic_cost
    screening_df$FP_Invasive_cost <- FP_Invasive_cost
    screening_df$TN_Invasive_cost <- TN_Invasive_cost
    screening_df$Imaging_cost <- Imaging_cost
    screening_df$Diagnostic_cost <- Diagnostic_cost
    screening_df$Invasive_cost <- Invasive_cost
    screening_df$FP_cost <- FP_cost
    screening_df$TN_cost <- TN_cost
    screening_df$FP_TN_cost <- FP_TN_cost
    
    # Save updated screening_df back to final_result_list
    final_result_list[["usual"]][[cancer]][["FP_TN_Screening_Summary"]] <- screening_df
    
    # Update overall sums (usual care)
    FP_imaging_cost_sum <- FP_imaging_cost_sum + FP_imaging_cost
    TN_imaging_cost_sum <- TN_imaging_cost_sum + TN_imaging_cost
    FP_diagnostic_cost_sum <- FP_diagnostic_cost_sum + FP_diagnostic_cost
    TN_diagnostic_cost_sum <- TN_diagnostic_cost_sum + TN_diagnostic_cost
    FP_Invasive_cost_sum <- FP_Invasive_cost_sum + FP_Invasive_cost
    TN_Invasive_cost_sum <- TN_Invasive_cost_sum + TN_Invasive_cost
    Imaging_cost_sum <- Imaging_cost_sum + Imaging_cost
    Diagnostic_cost_sum <- Diagnostic_cost_sum + Diagnostic_cost
    Invasive_cost_sum <- Invasive_cost_sum + Invasive_cost
    FP_cost_sum <- FP_cost_sum + FP_cost
    TN_cost_sum <- TN_cost_sum + TN_cost
    FP_TN_cost_sum <- FP_TN_cost_sum + FP_TN_cost
  }
  
  # Update overall Screening_Summary with total usual care costs
  overall_summary_df <- final_result_list[["Overall"]][["FP_TN_Screening_Summary"]]
  overall_summary_df$FP_imaging_cost <- FP_imaging_cost_sum
  overall_summary_df$TN_imaging_cost <- TN_imaging_cost_sum
  overall_summary_df$FP_diagnosis_cost <- FP_diagnostic_cost_sum
  overall_summary_df$TN_diagnosis_cost <- TN_diagnostic_cost_sum
  overall_summary_df$Imaging_cost <- Imaging_cost_sum
  overall_summary_df$Diagnostic_cost <- Diagnostic_cost_sum
  overall_summary_df$Invasive_cost <- Invasive_cost_sum
  overall_summary_df$FP_cost <- FP_cost_sum
  overall_summary_df$TN_cost <- TN_cost_sum
  overall_summary_df$FP_TN_cost <- FP_TN_cost_sum 
  final_result_list[["Overall"]][["FP_TN_Screening_Summary"]] <- overall_summary_df
  
  return(final_result_list)
}


apply_mced_test_cost <- function(final_result_list,
                                 mced_test_cost,           # unit cost for step‑1 test
                                 mced_secondary_test_cost, # unit cost for step‑2 test
                                 discount_factor) {
  ## helper
  discount_fn <- function(cycle) 1 / ((1 + discount_factor)^(cycle - 1))
  
  # ────────────────────────────────────────────────────────────────
  # 1) MCED_Cancer_cases_screening_summary  (unchanged)
  # ────────────────────────────────────────────────────────────────
  mced_sum <- final_result_list[["Overall"]][["MCED_Cancer_cases_screening_summary"]]
  if (!"Cycle" %in% names(mced_sum))
    mced_sum$Cycle <- seq_len(nrow(mced_sum))
  for (col in c("screen_step1", "screen_step2", "Total_overdiagnosed"))
    if (!col %in% names(mced_sum)) mced_sum[[col]] <- 0
  
  disc_vec <- discount_fn(mced_sum$Cycle)
  
  mced_sum$mced_test_cost_step1 <- (mced_sum$screen_step1 + mced_sum$Total_overdiagnosed) * mced_test_cost           * disc_vec
  mced_sum$mced_test_cost_step2 <- (mced_sum$screen_step2 + mced_sum$Total_overdiagnosed) * mced_test_cost * disc_vec
  mced_sum$mced_test_cost <- mced_sum$mced_test_cost_step1 + mced_sum$mced_test_cost_step2
  
  final_result_list[["Overall"]][["MCED_Cancer_cases_screening_summary"]] <- mced_sum
  
  
  # ────────────────────────────────────────────────────────────────
  # 2) mced_TN_FP  – add cost columns here as well
  # ────────────────────────────────────────────────────────────────
  mced_tn_fp <- final_result_list[["Overall"]][["mced_TN_FP"]]
  if (!"Cycle" %in% names(mced_tn_fp))
    mced_tn_fp$Cycle <- seq_len(nrow(mced_tn_fp))
  
  ## a) guarantee the step‑1 columns (0 if single‑step run)
  for (col in c("step1_FP_f","step1_FP_m","step1_TN_f","step1_TN_m"))
    if (!col %in% names(mced_tn_fp)) mced_tn_fp[[col]] <- 0
  for (col in c("FP_f","FP_m","True_TN"))
    if (!col %in% names(mced_tn_fp)) stop("mced_TN_FP is missing expected column: ", col)
  
  ## b) derive screen counts
  # step‑1 screens = everyone who took the first test
  mced_tn_fp$screen_step1 <- 0
  
  # step‑2 screens = step‑2 FP  +  step‑2 TN
  step2_FP <- with(mced_tn_fp, fp_counts)
  step2_TN <- with(mced_tn_fp, True_TN)
  mced_tn_fp$screen_step2 <- step2_FP + step2_TN
  
  ## c) discounted costs
  disc_vec_tn <- discount_fn(mced_tn_fp$Cycle)
  
  
  mced_tn_fp$mced_test_cost_step1 <- mced_tn_fp$screen_step1 * mced_test_cost           * disc_vec_tn
  mced_tn_fp$mced_test_cost_step2 <- mced_tn_fp$screen_step2 * mced_test_cost * disc_vec_tn
  mced_tn_fp$mced_test_cost <- mced_tn_fp$mced_test_cost_step1 + mced_tn_fp$mced_test_cost_step2
  
  final_result_list[["Overall"]][["mced_TN_FP"]] <- mced_tn_fp
  
  return(final_result_list)
}


