calculate_mced_tp_fn_compliance <- function(
    incident_matrix,
    sensitivity_data,
    final_result_list,
    screening_data,            
    case_type = "TP",          
    reduced_comply = 1,        
    cycles = mced_start_cycle:mced_end_cycle,
    mced_compliance_rate = mced_compliance_rate,
    two_step,
    overall_sens) {
  # Standardize column names in the incident matrix if necessary
  colnames(incident_matrix) <- gsub("\\.{2,}", ".", colnames(incident_matrix))
  
  # Define stages and initialize result matrices
  stages <- c("I", "II", "III", "IV")
  m_count_TP <- matrix(0, nrow = nrow(incident_matrix), ncol = ncol(incident_matrix),
                       dimnames = list(NULL, colnames(incident_matrix)))
  m_count_FN <- matrix(0, nrow = nrow(incident_matrix), ncol = ncol(incident_matrix),
                       dimnames = list(NULL, colnames(incident_matrix)))
  m_count_mced_NA <- matrix(0, nrow = nrow(incident_matrix), ncol = ncol(incident_matrix),
                            dimnames = list(NULL, colnames(incident_matrix)))
  
  # Identify cancer types and loop through each
  cancer_types <- unique(sensitivity_data$Cancer)
  for (cancer_type in cancer_types) {
    
    # --- 1) Get standard screening info for this cancer_type ---
    screening_info <- get_screening_info(
      cancer_type = cancer_type,
      screening_data = screening_data,
      case_type = "TP",          # e.g., "FN"
      reduced_comply = reduced_comply # if applying that logic
    )
    screening_start_cycle <- screening_info$start_age - 39
    screening_end_cycle   <- screening_info$end_age - 39
    screening_compliance <- screening_info$compliance
    screening_compliance <- as.numeric(screening_compliance)
    

    
    # Build the stage column names for this cancer
    stage_columns <- paste(stages, cancer_type, sep = "_")
    
    # Skip if these columns do not exist in incident_matrix
    if (!all(stage_columns %in% colnames(incident_matrix))) next
    
    # Create a data frame for this cancer
    cancer_df <- data.frame(
      Cycle = seq_len(nrow(incident_matrix)),
      I = incident_matrix[, stage_columns[stages == "I"]],
      II = incident_matrix[, stage_columns[stages == "II"]],
      III = incident_matrix[, stage_columns[stages == "III"]],
      IV = incident_matrix[, stage_columns[stages == "IV"]]
    )
    
    # --- 2) Determine final compliance for each row/cycle ---
    # We assume "Cycle" is effectively the person's age in years.
    # If you have an offset from cycle to actual age, adjust the logic accordingly.
    
    # By default, if the cycle is in the MCED range, use mced_compliance_rate;
    # otherwise, no MCED screening (0).
    cancer_df$final_compliance <- ifelse(
      cancer_df$Cycle %in% cycles, 
      mced_compliance_rate, 
      0
    )
    
    # Now override with the standard screening compliance IF:
    # (1) we are within [screening_start_age, screening_end_age], AND
    # (2) we are in the MCED cycles.
    # AND if standard compliance > MCED compliance.
    in_screening_age <- (cancer_df$Cycle >= screening_start_cycle & 
                           cancer_df$Cycle <= screening_end_cycle)
    in_mced_cycle    <- (cancer_df$Cycle %in% cycles)
    
    # Identify the overlap
    overlap_index <- which(in_screening_age & in_mced_cycle)
    
    # replicate to match the vector length pmax expects
    if (length(screening_compliance) == 1) {
      screening_compliance <- rep(screening_compliance,
                                  length(overlap_index))
    }
    
    
    # --- 3) Retrieve stage-specific sensitivity values ---
    sensitivity_for_cancer <- sensitivity_data[sensitivity_data$Cancer == cancer_type, ]
    sensitivity_map <- setNames(sensitivity_for_cancer$Sensitivity, sensitivity_for_cancer$Stage)
    
    # --- 4) Compute compliant vs non‑compliant counts by stage ---
    for (stg in stages) {
      base_cases         <- cancer_df[[stg]]
      compl              <- cancer_df$final_compliance
      
      compliant_cases    <- base_cases * compl        # N_c
      non_compliant_cases<- base_cases - compliant_cases
      
      screen_step1   <- screen_step2 <- screen_total <- overall_FN <- NULL
      
      if (two_step) {
        # ---- STEP A : overall MCED rule‑in ----
        overall_pos   <- compliant_cases * overall_sens
        overall_FN    <- compliant_cases - overall_pos
        
        # ---- STEP B : stage‑specific sensitivity ----
        S_stage       <- sensitivity_map[stg]
        TP_stage      <- overall_pos * S_stage
        stage_FN      <- overall_pos - TP_stage
        
        FN_stage      <- overall_FN + stage_FN
        
        ## -------- NEW: screen counters ----------
        screen_step1  <- compliant_cases          # first blood draw
        screen_step2  <- overall_pos              # second work‑up
        screen_total  <- screen_step1 + screen_step2
      } else {
        overall_FN    <- overall_FN <- rep(0, length(compliant_cases))
        S_stage       <- sensitivity_map[stg]
        TP_stage      <- compliant_cases * S_stage
        FN_stage      <- compliant_cases - TP_stage
        
        ## -------- NEW: screen counters ----------
        screen_step1  <- compliant_cases          # only one test
        screen_step2  <- 0
        screen_total  <- screen_step1
      }
      
      # store into dataframe (these column names already exist later in your code)
      cancer_df[[paste0("mced_compliant_", stg)]] <- compliant_cases
      cancer_df[[paste0("mced_non_compliant_", stg)]] <- non_compliant_cases
      cancer_df[[paste0("overall_FN_", stg)]] <- overall_FN
      cancer_df[[paste0("TP_", stg)]] <- TP_stage
      cancer_df[[paste0("FN_", stg)]] <- FN_stage
      cancer_df[[paste0("screen_step1_", stg)]]  <- screen_step1
      cancer_df[[paste0("screen_step2_", stg)]]  <- screen_step2
      cancer_df[[paste0("screen_total_", stg)]]  <- screen_total
    }
    
    
    # --- 6) Compute totals across stages within this cancer-specific dataframe ---
    total_columns <- c("mced_compliant", "mced_non_compliant", "overall_FN",  "TP", "FN")
    for (total_col in total_columns) {
      stage_specific_cols <- paste0(total_col, "_", stages)
      cancer_df[[total_col]] <- rowSums(
        cancer_df[, stage_specific_cols, drop = FALSE], 
        na.rm = TRUE
      )
    }
    
    screen_cols <- c("screen_step1", "screen_step2", "screen_total")
    for (scr in screen_cols) {
      stage_cols <- paste0(scr, "_", stages)
      cancer_df[[scr]] <- rowSums(cancer_df[, stage_cols, drop = FALSE], na.rm = TRUE)
    }
    # --- 7) Fill the output matrices with stage-specific values (TP, FN, NA) ---
    # TP
    m_count_TP[, stage_columns[stages == "I"]]   <- cancer_df$TP_I
    m_count_TP[, stage_columns[stages == "II"]]  <- cancer_df$TP_II
    m_count_TP[, stage_columns[stages == "III"]] <- cancer_df$TP_III
    m_count_TP[, stage_columns[stages == "IV"]]  <- cancer_df$TP_IV
    
    # FN
    m_count_FN[, stage_columns[stages == "I"]]   <- cancer_df$FN_I
    m_count_FN[, stage_columns[stages == "II"]]  <- cancer_df$FN_II
    m_count_FN[, stage_columns[stages == "III"]] <- cancer_df$FN_III
    m_count_FN[, stage_columns[stages == "IV"]]  <- cancer_df$FN_IV
    
    # MCED_NA (noncompliant)
    m_count_mced_NA[, stage_columns[stages == "I"]]   <- cancer_df$mced_non_compliant_I
    m_count_mced_NA[, stage_columns[stages == "II"]]  <- cancer_df$mced_non_compliant_II
    m_count_mced_NA[, stage_columns[stages == "III"]] <- cancer_df$mced_non_compliant_III
    m_count_mced_NA[, stage_columns[stages == "IV"]]  <- cancer_df$mced_non_compliant_IV
    
    # --- 8) Store this cancer's data frame in the final_result_list ---
    final_result_list[["MCED"]][[cancer_type]] <- list(
      MCED_Cancer_cases = cancer_df
    )
  }
  
  # --- 9) Create a summary data frame for Overall ---
  summary_df <- data.frame(Cycle = seq_len(nrow(incident_matrix)))
  
  for (stg in stages) {
    # Add stage-specific totals for compliance, non-compliance, TP, and FN
    summary_df[[paste0(stg)]] <- Reduce("+", lapply(final_result_list[["MCED"]], function(x) x$MCED_Cancer_cases[[stg]]))
    summary_df[[paste0("mced_compliant_", stg)]] <- Reduce("+", lapply(final_result_list[["MCED"]], function(x) x$MCED_Cancer_cases[[paste0("mced_compliant_", stg)]]))
    summary_df[[paste0("mced_non_compliant_", stg)]] <- Reduce("+", lapply(final_result_list[["MCED"]], function(x) x$MCED_Cancer_cases[[paste0("mced_non_compliant_", stg)]]))
    summary_df[[paste0("overall_FN_", stg)]] <- Reduce("+", lapply(final_result_list[["MCED"]], function(x) x$MCED_Cancer_cases[[paste0("overall_FN_", stg)]]))
    summary_df[[paste0("TP_", stg)]] <- Reduce("+", lapply(final_result_list[["MCED"]], function(x) x$MCED_Cancer_cases[[paste0("TP_", stg)]]))
    summary_df[[paste0("FN_", stg)]] <- Reduce("+", lapply(final_result_list[["MCED"]], function(x) x$MCED_Cancer_cases[[paste0("TP_", stg)]]))
  }
  
  for (stg in stages) {
    summary_df[[paste0("screen_step1_", stg)]] <- Reduce("+", lapply(final_result_list[["MCED"]], function(x) x$MCED_Cancer_cases[[paste0("screen_step1_", stg)]]))
    summary_df[[paste0("screen_step2_", stg)]] <- Reduce("+", lapply(final_result_list[["MCED"]], function(x) x$MCED_Cancer_cases[[paste0("screen_step2_", stg)]]))
    summary_df[[paste0("screen_total_", stg)]] <- Reduce("+", lapply(final_result_list[["MCED"]], function(x) x$MCED_Cancer_cases[[paste0("screen_total_", stg)]]))
  }
  
  # and the grand totals
  for (scr in c("screen_step1", "screen_step2", "screen_total")) {
    stage_cols <- paste0(scr, "_", stages)
    summary_df[[scr]] <- rowSums(summary_df[, stage_cols, drop = FALSE], na.rm = TRUE)
  }
  
  # Totals across stages
  total_columns <- c("mced_compliant", "mced_non_compliant", "overall_FN", "TP", "FN")
  for (total_col in total_columns) {
    stage_specific_cols <- paste0(total_col, "_", stages)
    summary_df[[total_col]] <- rowSums(
      summary_df[, stage_specific_cols, drop = FALSE], 
      na.rm = TRUE
    )
  }
  
  
  # Update the Overall section of final_result_list
  final_result_list[["Overall"]][["MCED_Cancer_cases_screening_summary"]] <- summary_df
  
  # Return results
  return(list(
    final_result_list = final_result_list,
    m_count_TP = m_count_TP,
    m_count_FN = m_count_FN,
    m_count_mced_NA = m_count_mced_NA
  ))
}
