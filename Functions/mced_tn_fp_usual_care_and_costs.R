discount <- function(cycle, discount_factor) {
  1 / ((1 + discount_factor)^(cycle - 1))
}


apply_mced_usual_care_and_costs <- function(final_result_list,
                                            screening_data,
                                            discount_factor,
                                            screening_cancers,
                                            tn_skip_factor,
                                            invasive_multiplier) {
  
  harmonized_cols <- c(
    "FP_fu_compliant", "FP_i", "FP_ni", "FP_FP", "FP_TN",
    "FP_FP_i", "FP_FP_ni",
    "TN", "TN_c", "TN_nc",
    "TN_c_TN", "TN_c_FP", "TN_c_FP_i", "TN_c_FP_ni",
    "m_invasive", "u_screened", "u_unscreened",
    "u_FP", "u_TN", "u_invasive"
  )
  
  
  overall_df <- final_result_list[["Overall"]][["mced_TN_FP"]]
  n_cycles   <- nrow(overall_df)
  discount_vec <- sapply(seq_len(n_cycles), discount, discount_factor = discount_factor)
  final_result_list[["Overall"]][["mced_TN_FP"]] <- overall_df
  
  mced_cancers <- names(final_result_list[["MCED"]])
  for (cancer in mced_cancers) {
    df <- final_result_list[["MCED"]][[cancer]][["mced_TN_FP"]]
    if (is.null(df)) next
    n_cycles <- nrow(df)
    ages     <- (1:n_cycles) + 39
    
    scr_row <- screening_data[screening_data$Cancer == cancer, , drop = FALSE]
    discount_vec <- sapply(seq_len(n_cycles), discount, discount_factor = discount_factor)
    
    if (!(cancer %in% screening_cancers)) {
      # ensure these fields exist
      if (!"FP_fu_compliant" %in% colnames(df)) {
        df$FP_fu_compliant <- 0
        message("  -> FP_fu_compliant was missing; set to 0")
      }
      
      incidental_cost <- scr_row$incidental_cost
      invasive_cost   <- scr_row$invasive_cost
      
      # recalc costs
      df$FP_i <- df$FP_fu_compliant * scr_row$Invasive_procedure
      df$FP_ni <- df$FP_fu_compliant * (1- scr_row$Invasive_procedure)
      df$m_invasive <- df$FP_i
      
      Imaging_cost    <- df$FP_fu_compliant * incidental_cost * discount_vec
      Diagnostic_cost <- rep(0, n_cycles)
      Invasive_cost   <- df$FP_i * invasive_cost * discount_vec
      
      
      # assemble and store
      mced_TN_cost <- rep(0, n_cycles)  # always zero here
      final_result_list[["MCED"]][[cancer]][["mced_tn_fp_cost"]] <-
        data.frame(
          Imaging_cost    = Imaging_cost,
          Diagnostic_cost = Diagnostic_cost,
          Invasive_cost   = Invasive_cost,
          mced_TN_cost    = mced_TN_cost,
          mced_FP_Imaging = Imaging_cost,
          mced_FP_Diagnostic = Diagnostic_cost,
          mced_FP_Invasive = Invasive_cost,
          stringsAsFactors = FALSE
        )
      
      for (col in harmonized_cols) {
        if (!col %in% colnames(df)) {
          df[[col]] <- 0
        }
      }
      
      final_result_list[["MCED"]][[cancer]][["mced_TN_FP"]] <- df
      next
    }
    
    # Screening cancer logic:
    Specificity <- scr_row$Specificity
    Eligibility <- scr_row$Eligibility
    Invasive_procedure <- scr_row$Invasive_procedure
    Screening_cost <- scr_row$Screening_cost
    FP_cost <- scr_row$FP_cost
    invasive_cost <- scr_row$invasive_cost
    incidental_cost <- scr_row$incidental_cost
    Start_age <- scr_row$Start_age
    End_age <- scr_row$end_age
    if (tn_skip_factor == 0) {
      Compliance_freq <- 0
    } else if (tn_skip_factor != 1) {
      Compliance_freq <- scr_row$Compliance_frequency / tn_skip_factor
    } else {
      Compliance_freq <- scr_row$Compliance_frequency
    }
    
    if (!"FP_fu_compliant" %in% colnames(df)) {
      stop("There is an issue with FP followups for MCED")
    }
    if (!"TN" %in% colnames(df)) {
      df$TN <- 0
    }
    

    
    if (cancer == "prostate") {
      df$FP_i <-  df$FP_fu_compliant * Invasive_procedure
      df$FP_ni <- df$FP_fu_compliant * (1 - Invasive_procedure)
      
      df$FP_FP <- 0
      df$FP_TN <- 0
      
      df$FP_FP_i <- 0
      df$FP_FP_ni <- 0
      
      df$m_invasive <- 0

    } else {
      
      df$FP_i <- 0 # Only for non-screening cancers
      df$FP_ni <- 0
      
      # Split FP_fu_compliant into usual care FP and TN
      
      df$FP_FP <- df$FP_fu_compliant * (1 - Specificity)
      df$FP_TN <- df$FP_fu_compliant * Specificity
      
      FP_FP_invasive_factor <- Invasive_procedure * (1 + invasive_multiplier)
      # FP_FP_i, FP_FP_ni
      df$FP_FP_i <- df$FP_FP * FP_FP_invasive_factor
      df$FP_FP_ni <- df$FP_FP * (1 - FP_FP_invasive_factor)
      
      df$m_invasive <- df$FP_i
    }
    
    
    TN_c <- numeric(n_cycles)
    TN_nc <- numeric(n_cycles)
    TN_c_TN <- numeric(n_cycles)
    TN_c_FP <- numeric(n_cycles)
    TN_c_FP_i <- numeric(n_cycles)
    TN_c_FP_ni <- numeric(n_cycles)
    
    for (i in seq_len(n_cycles)) {
      age_i <- ages[i]
      TN_i <- df$TN[i]
      if (age_i >= Start_age && age_i <= End_age) {
        tested_prop <- Eligibility * Compliance_freq
        TN_c[i] <- TN_i * tested_prop
        TN_nc[i] <- TN_i * (1 - tested_prop)
        
        TN_c_TN[i] <- TN_c[i] * Specificity
        TN_c_FP[i] <- TN_c[i] * (1 - Specificity)
        
        TN_c_FP_i[i] <- TN_c_FP[i] * Invasive_procedure
        TN_c_FP_ni[i] <- TN_c_FP[i] * (1 - Invasive_procedure)
      } else {
        TN_c[i] <- 0
        TN_nc[i] <- TN_i
        TN_c_TN[i] <- 0
        TN_c_FP[i] <- 0
        TN_c_FP_i[i] <- 0
        TN_c_FP_ni[i] <- 0
      }
    }
    
    df$TN_c <- TN_c
    df$TN_nc <- TN_nc
    df$TN_c_TN <- TN_c_TN
    df$TN_c_FP <- TN_c_FP
    df$TN_c_FP_i <- TN_c_FP_i
    df$TN_c_FP_ni <- TN_c_FP_ni
    
    
    df$u_screened <- df$FP_fu_compliant + df$TN_c
    df$u_unscreened <- df$TN_nc + (df$FP_fu - df$FP_fu_compliant)
    df$u_FP <- df$FP_FP + df$TN_c_FP
    df$u_TN <- df$FP_TN + df$TN_c_TN
    df$u_invasive <- df$FP_FP_i + df$TN_c_FP_i + df$FP_i
    
    
    # Costs
    Imaging_cost <- (df$FP_FP_i + df$FP_FP_ni + df$FP_TN + df$TN_c_TN + df$TN_c_FP_i + df$TN_c_FP_ni)*Screening_cost * discount_vec
    Diagnostic_cost <- (df$FP_FP_i + df$FP_FP_ni + df$TN_c_FP_i + df$TN_c_FP_ni + df$FP_i + df$FP_ni)*FP_cost * discount_vec
    Invasive_cost <- (df$FP_FP_i + df$TN_c_FP_i + df$FP_i)*invasive_cost * discount_vec
    Total_FP_TN_cost <- Imaging_cost + Diagnostic_cost + Invasive_cost
    
    # MCED-FP costs
    mced_FP_Imaging <- (df$FP_fu_compliant)*Screening_cost * discount_vec
    mced_FP_Diagnostic <- (df$FP_FP_i + df$FP_FP_ni + df$FP_i + df$FP_ni)*FP_cost * discount_vec
    mced_FP_Invasive <- (df$FP_FP_i + df$FP_i)*invasive_cost * discount_vec
    mced_FP_cost <- mced_FP_Imaging + mced_FP_Diagnostic + mced_FP_Invasive
    
    # MCED-TN costs
    mced_TN_Imaging <- (df$TN_c)*Screening_cost * discount_vec
    mced_TN_Diagnostic <- (df$TN_c_FP)*FP_cost * discount_vec
    mced_TN_Invasive <- (df$TN_c_FP_i)*invasive_cost * discount_vec
    mced_TN_cost <- mced_TN_Imaging + mced_TN_Diagnostic + mced_TN_Invasive  # Fixed typo here
    
    final_result_list[["MCED"]][[cancer]][["mced_TN_FP"]] <- df
    
    final_result_list[["MCED"]][[cancer]][["mced_tn_fp_cost"]] <- data.frame(
      Imaging_cost = Imaging_cost,
      Diagnostic_cost = Diagnostic_cost,
      Invasive_cost = Invasive_cost,
      Total_FP_TN_cost = Total_FP_TN_cost,
      
      mced_FP_Imaging = mced_FP_Imaging,
      mced_FP_Diagnostic = mced_FP_Diagnostic,
      mced_FP_Invasive = mced_FP_Invasive,
      mced_FP_cost = mced_FP_cost,
      
      mced_TN_Imaging = mced_TN_Imaging,
      mced_TN_Diagnostic = mced_TN_Diagnostic,
      mced_TN_Invasive = mced_TN_Invasive,
      mced_TN_cost = mced_TN_cost,
      stringsAsFactors = FALSE
    )
  }
  
  return(final_result_list)
}


summarize_overall_results <- function(final_result_list) {
  
  overall_df <- final_result_list[["Overall"]][["mced_TN_FP"]]
  n_cycles <- nrow(overall_df)
  
  # Initialize aggregates
  count_cols <- c(
    "m_invasive", "u_screened", "u_unscreened", "u_FP", "u_TN", "u_invasive", "FP_i", "FP_ni",
    "FP_FP", "FP_TN", "FP_FP_i", "FP_FP_ni",
    "TN_c", "TN_nc", "TN_c_TN", "TN_c_FP", "TN_c_FP_i", "TN_c_FP_ni"
  )
  
  # Initialize sum holders
  count_sums <- matrix(0, nrow = n_cycles, ncol = length(count_cols))
  colnames(count_sums) <- count_cols
  
  cost_cols <- c("Imaging_cost", "Diagnostic_cost", "Invasive_cost",
                 "mced_FP_Imaging", "mced_FP_Diagnostic", "mced_FP_Invasive",
                 "mced_TN_Imaging", "mced_TN_Diagnostic", "mced_TN_Invasive",
                 "mced_FP_cost", "mced_TN_cost", "Total_FP_TN_cost")
  
  cost_sums <- matrix(0, nrow = n_cycles, ncol = length(cost_cols))
  colnames(cost_sums) <- cost_cols
  
  mced_cancers <- names(final_result_list[["MCED"]])
  
  get_col_or_zero <- function(df, colname, n) {
    if (colname %in% names(df)) df[[colname]] else rep(0, n)
  }
  
  for (cancer in mced_cancers) {
    cancer_df <- final_result_list[["MCED"]][[cancer]][["mced_TN_FP"]]
    if (is.null(cancer_df)) next
    
    n_cycles <- nrow(cancer_df)
    for (col in count_cols) {
      count_sums[, col] <- count_sums[, col] + get_col_or_zero(cancer_df, col, n_cycles)
    }
    
    cancer_cost_df <- final_result_list[["MCED"]][[cancer]][["mced_tn_fp_cost"]]
    
    if (is.null(cancer_cost_df)) next
    for (cc in cost_cols) {
      if (cc %in% colnames(cancer_cost_df)) {
        cost_sums[, cc] <- cost_sums[, cc] + cancer_cost_df[[cc]]
      } else {
        cost_sums[, cc] <- cost_sums[, cc] + 0
      }
    }
  }
  
  # Update overall_df with aggregated counts
  for (col in count_cols) {
    overall_df[[col]] <- count_sums[, col]
  }
  
  final_result_list[["Overall"]][["mced_TN_FP"]] <- overall_df
  
  # Create a data frame for overall costs
  overall_cost_df <- as.data.frame(cost_sums, stringsAsFactors = FALSE)
  final_result_list[["Overall"]][["mced_TN_FP_cost"]] <- overall_cost_df
  
  return(final_result_list)
}




