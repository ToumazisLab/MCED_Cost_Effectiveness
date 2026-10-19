shift_backward <- function(
    cancer_types,
    df,
    ss_matrix_list,
    opt_val = 1,
    ts_matrix_list,
    scenario,
    mced_start_cycle,
    mced_end_cycle,
    final_result_list
) {
  final_df <- NULL
  shifted_cases_matrix <- NULL
  stages <- c("I", "II", "III", "IV")
  
  # Sankey (unchanged structure: cancer_type, cycle, source_stage, target_stage, n_cases)
  sankey_list <- vector("list", length = 0)
  sankey_counter <- 1
  
  # NEW: transition records carrying timing for the floored death calc
  tp_trans_list <- vector("list", length = 0)
  tp_trans_counter <- 1
  
  if (opt_val %in% c(1, 2)) {
    scenario <- "time_shift_matrix_1"
  } else if (opt_val %in% c(3, 4)) {
    scenario <- "time_shift_matrix_2"
  } else if (opt_val %in% c(5, 6)) {
    scenario <- "time_shift_matrix_3"
  } else if (opt_val %in% c(7, 8)) {
    scenario <- "time_shift_matrix_4"
  } else {
    scenario <- "unknown_scenario"
  }
  
  for (cancer_type in cancer_types) {
    
    cancer_df <- df[, grep(paste0("_", cancer_type, "$"), colnames(df)), drop = FALSE]
    
    # NEW: force I, II, III, IV column order so positional indexing is safe
    wanted <- paste0(stages, "_", cancer_type)
    wanted <- wanted[wanted %in% colnames(cancer_df)]
    cancer_df <- cancer_df[, wanted, drop = FALSE]
    
    if (!cancer_type %in% names(ss_matrix_list[[opt_val]]) ||
        !cancer_type %in% names(ts_matrix_list)) {
      stop("Cancer type matrix not found in list!")
    }
    
    ss_matrix <- ss_matrix_list[[opt_val]][[cancer_type]]
    ts_matrix <- ts_matrix_list[[cancer_type]][[scenario]]
    
    temp_df <- matrix(0, nrow = nrow(cancer_df), ncol = ncol(cancer_df),
                      dimnames = list(NULL, colnames(cancer_df)))
    shifted_cases_df <- matrix(0, nrow = nrow(cancer_df), ncol = ncol(cancer_df),
                               dimnames = list(NULL, colnames(cancer_df)))
    
    for (cycle in mced_start_cycle:mced_end_cycle) {
      original_cases <- cancer_df[cycle, ]
      
      if (all(original_cases == 0)) {
        temp_df[cycle, ] <- original_cases
        next
      }
      
      adjusted_cases <- original_cases  # placeholder for compliance/sensitivity logic
      
      shifted_cases <- numeric(length(stages))
      
      for (i in seq_along(stages)) {
        from_stage_count <- adjusted_cases[i]
        
        for (j in 1:i) {
          n_ij <- from_stage_count * ss_matrix[j, i]
          shifted_cases[j] <- shifted_cases[j] + n_ij
          
          if (n_ij > 0) {
            # --- existing sankey record (unchanged) ---
            sankey_list[[sankey_counter]] <- list(
              cancer_type  = cancer_type,
              cycle        = cycle,
              source_stage = stages[i],
              target_stage = stages[j],
              n_cases      = n_ij
            )
            sankey_counter <- sankey_counter + 1
            
            # --- NEW: transition record with destination timing ---
            # Time shift is a function of the DESTINATION stage j (matches the
            # pooled placement below: temp_df[new_cycle, j] gets all of stage j).
            L_j  <- ts_matrix[j, j]
            nc_j <- max(mced_start_cycle, ceiling(cycle - L_j))
            tp_trans_list[[tp_trans_counter]] <- list(
              cancer_type  = cancer_type,
              origin_cycle = cycle,
              source_stage = stages[i],   # usual-arm stage (for baseline survival)
              dest_stage   = stages[j],   # shifted-down stage (for new survival)
              new_cycle    = nc_j,
              n_cases      = n_ij
            )
            tp_trans_counter <- tp_trans_counter + 1
          }
        }
      }
      
      shifted_cases_df[cycle, ] <- shifted_cases
      
      # Time-shift placement (UNCHANGED)
      for (stage_idx in seq_along(stages)) {
        if (shifted_cases[stage_idx] > 0) {
          time_shift <- ts_matrix[stage_idx, stage_idx]
          new_cycle <- max(mced_start_cycle, ceiling(cycle - time_shift))
          temp_df[new_cycle, stage_idx] <- temp_df[new_cycle, stage_idx] + shifted_cases[stage_idx]
        }
      }
      
      # Non-shifted portion stays put (currently zero because adjusted == original)
      non_shifted_cases <- original_cases - adjusted_cases
      temp_df[cycle, ] <- temp_df[cycle, ] + non_shifted_cases
      
      # NEW: record non-shifted cases so their (baseline) deaths are also captured
      # if/when you implement partial compliance. No-op while non_shifted == 0.
      for (s_idx in seq_along(stages)) {
        if (non_shifted_cases[s_idx] > 0) {
          tp_trans_list[[tp_trans_counter]] <- list(
            cancer_type  = cancer_type,
            origin_cycle = cycle,
            source_stage = stages[s_idx],
            dest_stage   = stages[s_idx],
            new_cycle    = cycle,
            n_cases      = non_shifted_cases[s_idx]
          )
          tp_trans_counter <- tp_trans_counter + 1
        }
      }
    }
    
    adjusted_cancer_df <- temp_df
    
    if (cancer_type %in% names(final_result_list[["MCED"]])) {
      cancer_data <- final_result_list[["MCED"]][[cancer_type]][["MCED_Cancer_cases"]]
      
      for (stage in stages) {
        column_name <- paste0(stage, "_", cancer_type)
        shifted_col <- paste0(stage, "_TP_shifted")
        cancer_data[[shifted_col]] <- shifted_cases_df[, column_name]
      }
      
      cancer_data$Total_TP_shifted <- cancer_data$`I_TP_shifted` +
        cancer_data$`II_TP_shifted` +
        cancer_data$`III_TP_shifted` +
        cancer_data$`IV_TP_shifted`
      
      final_result_list[["MCED"]][[cancer_type]]$MCED_Cancer_cases <- cancer_data
    }
    
    final_df <- cbind(final_df, temp_df)
    shifted_cases_matrix <- cbind(shifted_cases_matrix, shifted_cases_df)
  }
  
  sankey_records <- do.call(rbind.data.frame, c(sankey_list, stringsAsFactors = FALSE))
  
  # NEW: assemble transition table (empty-safe)
  tp_transitions <- if (length(tp_trans_list) > 0) {
    do.call(rbind.data.frame, c(tp_trans_list, stringsAsFactors = FALSE))
  } else {
    data.frame(cancer_type = character(0), origin_cycle = integer(0),
               source_stage = character(0), dest_stage = character(0),
               new_cycle = integer(0), n_cases = numeric(0))
  }
  
  summary_df <- final_result_list[["Overall"]][["MCED_Cancer_cases_screening_summary"]]
  for (stage in stages) {
    shifted_col <- paste0("Total_TP_shifted_", stage)
    cancer_columns <- grep(paste0("^", stage, "_"), colnames(shifted_cases_matrix), value = TRUE)
    summary_df[[shifted_col]] <- rowSums(shifted_cases_matrix[, cancer_columns, drop = FALSE])
  }
  summary_df$Total_TP_shifted <- rowSums(shifted_cases_matrix)
  final_result_list[["Overall"]][["MCED_Cancer_cases_screening_summary"]] <- summary_df
  
  return(list(
    final_result_list    = final_result_list,
    updated_matrix       = final_df,
    shifted_cases_matrix = shifted_cases_matrix,
    sankey_df            = sankey_records,
    tp_transitions       = tp_transitions   # <- NEW
  ))
}


# Discounted post-diagnosis cost per CASE. Mirrors decision_tree()'s inner loop:
#   full_years capped at 5; year k discounted at exponent (start_row-1)+(k-1);
#   partial year only when full_years < 5.
.cohort_post_cost <- function(start_row, stage_col, surv_years,
                              cost_list, discount_rate, n_years = 5) {
  if (is.na(surv_years) || surv_years <= 0) return(0)
  cost_vec <- cost_list[[stage_col]]
  if (is.null(cost_vec)) return(0)
 
  full_years   <- min(floor(surv_years), n_years)
  partial_year <- surv_years - full_years
  total <- 0
 
  if (full_years >= 1) {
    fc <- unlist(cost_vec[as.character(seq_len(full_years))])
    if (length(fc) > 0) {
      k <- seq_along(fc)
      total <- total + sum(fc / ((1 + discount_rate) ^ ((start_row - 1) + (k - 1))))
    }
  }
  if (partial_year > 0 && full_years < n_years) {
    pc <- cost_vec[[as.character(full_years + 1)]]
    if (!is.null(pc)) {
      disc_part <- (start_row - 1) + full_years
      total <- total + pc * partial_year / ((1 + discount_rate) ^ disc_part)
    }
  }
  total
}
 
# Discounted QALYs per CASE. Mirrors calculate_total_qalys()'s inner loop:
#   full_years = floor(surv) (uncapped); age-adjusted discounted utility from
#   expanded_utility starting at start_row; stage multipliers for years 1-5,
#   1.0 thereafter; last (partial) year scaled by the fractional part.
.cohort_qalys <- function(start_row, stage_col, surv_years,
                          expanded_utility, utility_lookup) {
  if (is.na(surv_years) || surv_years <= 0) return(0)
 
  full_years   <- floor(surv_years)
  partial_year <- surv_years - full_years
  max_years    <- full_years + if (partial_year > 0) 1 else 0
  if (max_years < 1) return(0)
 
  max_age_index <- min(start_row + max_years - 1, nrow(expanded_utility))
  if (start_row > nrow(expanded_utility)) return(0)
 
  age_util <- expanded_utility$Discounted_avg[start_row:max_age_index]
 
  utility_multipliers_for_stage <- utility_lookup[[stage_col]]
  utility_multipliers <- rep(1, length(age_util))
  if (length(utility_multipliers_for_stage) > 0) {
    k <- seq_len(min(5, length(utility_multipliers)))
    utility_multipliers[k] <-
      utility_multipliers_for_stage[seq_len(min(5,
        length(utility_multipliers_for_stage), length(utility_multipliers)))]
  }
 
  if (partial_year > 0 && length(age_util) > 0) {
    age_util[length(age_util)] <- age_util[length(age_util)] * partial_year
  }
 
  sum(age_util * utility_multipliers)
}
 
# --- effective survival from the detection cycle -----------------------------
.effective_survival <- function(method, new_surv, old_surv, eff_lead) {
  switch(method,
    "naive"       = new_surv,
    "floor_max"   = max(new_surv, eff_lead + old_surv),               # Approach 1
    "add_shift"   = eff_lead + max(new_surv, old_surv),               # Approach 2 (+safety floor)
    "max_new_old" = max(new_surv, old_surv)                           # crossover-guarded, NOT lead-time floored
  )
}
 
# --- driver ------------------------------------------------------------------
calculate_tp_post_diagnosis_floored <- function(tp_transitions,
                                                transposed_survival,
                                                cost_per_year,
                                                utility_long,
                                                expanded_utility,
                                                discount_rate,
                                                n_cycles,
                                                col_names,
                                                survival_method = c("floor_max", "add_shift", "naive", "max_new_old"),
                                                cost_horizon = c("effective")) {
  survival_method <- match.arg(survival_method)
  cost_horizon    <- match.arg(cost_horizon)
  
  n_cols      <- length(col_names)
  death_rows  <- ceiling(n_cycles * 1.5)
  
  costs   <- matrix(0, nrow = n_cycles,   ncol = n_cols, dimnames = list(NULL, col_names))
  ly      <- matrix(0, nrow = n_cycles,   ncol = n_cols, dimnames = list(NULL, col_names))
  qaly    <- matrix(0, nrow = n_cycles,   ncol = n_cols, dimnames = list(NULL, col_names))
  deaths  <- matrix(0, nrow = death_rows, ncol = n_cols, dimnames = list(NULL, col_names))
  
  # NEW: empty-safe template for the per-transition diagnostic table
  empty_transition_table <- function() data.frame(
    cancer_type        = character(0),
    original_stage     = character(0),
    shifted_stage      = character(0),
    original_cycle     = numeric(0),
    shifted_cycle      = numeric(0),
    original_survival  = numeric(0),
    new_survival       = numeric(0),
    lead_time          = numeric(0),
    effective_survival = numeric(0),
    n_cases            = numeric(0),
    stringsAsFactors   = FALSE
  )
  
  empty_out <- function() list(
    total_costs_incidence = as.data.frame(costs),
    survival_df           = as.data.frame(ly),
    qaly_wide             = as.data.frame(qaly),
    deaths                = deaths,
    transition_table      = empty_transition_table()   # <- NEW
  )
  if (is.null(tp_transitions) || nrow(tp_transitions) == 0) return(empty_out())
  
  # Lookups, built once (same construction as your originals)
  cost_list <- cost_per_year %>%
    dplyr::mutate(Year = as.integer(gsub("Year ", "", Year))) %>%
    dplyr::group_by(Cancer_Stage) %>%
    dplyr::summarize(cost_vec = list(setNames(Cost, Year)), .groups = "drop") %>%
    tibble::deframe()
  
  utility_lookup <- split(utility_long$utility, utility_long$Cancer_Stage)
  
  surv_cols <- colnames(transposed_survival)
  surv_rows <- nrow(transposed_survival)
  
  trans_list    <- vector("list", length = 0)   # <- NEW
  trans_counter <- 1                             # <- NEW
  
  for (r in seq_len(nrow(tp_transitions))) {
    cancer  <- tp_transitions$cancer_type[r]
    c0      <- as.numeric(tp_transitions$origin_cycle[r])
    nc      <- as.numeric(tp_transitions$new_cycle[r])
    n       <- as.numeric(tp_transitions$n_cases[r])
    src_col <- paste0(tp_transitions$source_stage[r], "_", cancer)
    dst_col <- paste0(tp_transitions$dest_stage[r],   "_", cancer)
    
    if (n <= 0) next
    if (!(dst_col %in% col_names)) next
    if (!(src_col %in% surv_cols) || !(dst_col %in% surv_cols)) next
    if (c0 < 1 || nc < 1 || c0 > surv_rows || nc > surv_rows) next
    
    old_surv <- as.numeric(transposed_survival[[src_col]][c0])
    new_surv <- as.numeric(transposed_survival[[dst_col]][nc])
    if (length(old_surv) != 1 || is.na(old_surv)) old_surv <- 0
    if (length(new_surv) != 1 || is.na(new_surv)) new_surv <- 0
    
    eff_lead <- c0 - nc
    eff_surv <- .effective_survival(survival_method, new_surv, old_surv, eff_lead)
    
    # cost horizon choice
    cost_surv <- if (cost_horizon == "clinical") new_surv else eff_surv
    
    # NEW: record the diagnostic row for this transition
    trans_list[[trans_counter]] <- data.frame(
      cancer_type          = cancer,
      original_stage       = tp_transitions$source_stage[r],
      new_stage            = tp_transitions$dest_stage[r],
      original_cycle       = c0,
      new_cycle            = nc,
      time_shift_length    = eff_lead,                # c0 - nc
      original_survival    = old_surv,
      new_survival         = new_surv,
      delta_survival       = new_surv - old_surv,     # new - original
      orig_surv_plus_shift = old_surv + eff_lead,     # the floor_max candidate
      effective_survival   = eff_surv,
      n_cases              = n,
      stringsAsFactors     = FALSE
    )
    trans_counter <- trans_counter + 1
    
    # --- accrue at the detection cycle (nc), destination stage -------------
    if (nc <= n_cycles) {
      ly[nc, dst_col]    <- ly[nc, dst_col]    + n * eff_surv
      qaly[nc, dst_col]  <- qaly[nc, dst_col]  + n * .cohort_qalys(nc, dst_col, eff_surv,
                                                                   expanded_utility, utility_lookup)
      costs[nc, dst_col] <- costs[nc, dst_col] + n * .cohort_post_cost(nc, dst_col, cost_surv,
                                                                       cost_list, discount_rate)
    }
    
    # --- death cycle (identical to calculate_tp_deaths_with_floor) ---------
    death_cycle <- nc + floor(eff_surv)
    if (death_cycle >= 1 && death_cycle <= death_rows) {
      deaths[death_cycle, dst_col] <- deaths[death_cycle, dst_col] + n
    }
  }
  
  # NEW: assemble the transition table (empty-safe)
  transition_table <- if (length(trans_list) > 0) {
    do.call(rbind, trans_list)          # all rows are 1-row data.frames now
  } else {
    empty_transition_table()
  }
  
  list(
    total_costs_incidence = as.data.frame(costs),  # -> m_post_diagnosis_tp
    survival_df           = as.data.frame(ly),     # -> m_survival_years_tp
    qaly_wide             = as.data.frame(qaly),    # -> m_utility_cancers_wide_tp
    deaths                = deaths,                 # -> m_deaths_tp
    transition_table      = transition_table        # <- NEW
  )
}

expand_sankey_with_FN_NA <- function(sankey_df, m_count_FN, m_count_mced_NA, cancer_list) {
  
  # Get the maximum cycle number from FN matrix
  total_cycles <- nrow(m_count_FN)
  
  # Dataframe to accumulate FN+NA identity flows
  fnna_list <- vector("list", length = 0)
  fnna_counter <- 1
  
  for (cancer_type in cancer_list) {
    stage_cols <- grep(paste0("_", cancer_type, "$"), colnames(m_count_FN))
    
    if (length(stage_cols) == 0) {
      stop(sprintf("No matching columns for cancer '%s'", cancer_type))
    }
    
    for (cycle in 1:total_cycles) {
      
      for (col_name in colnames(m_count_FN)[stage_cols]) {
        
        # Extract stage from colname
        stage_label <- sub(paste0("_", cancer_type), "", col_name)
        
        # Sum FN and NA counts for the current cycle-stage
        n_fn <- m_count_FN[cycle, col_name]
        n_na <- m_count_mced_NA[cycle, col_name]
        
        total_cases <- n_fn + n_na
        
        if (total_cases > 0) {
          fnna_list[[fnna_counter]] <- list(
            cancer_type  = cancer_type,
            cycle        = cycle,
            source_stage = stage_label,
            target_stage = stage_label,
            n_cases      = total_cases
          )
          fnna_counter <- fnna_counter + 1
        }      
      }
    }
  }
  
  fnna_df <- do.call(rbind.data.frame, c(fnna_list, stringsAsFactors = FALSE))
  
  # Now combine your original sankey_df with this expanded fnna_df
  combined_df <- rbind(sankey_df, fnna_df)
  
  # Ensure combined_df covers all cycles from 1:total_cycles
  return(combined_df)
}
