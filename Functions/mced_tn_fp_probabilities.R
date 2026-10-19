# Function to generate compliant_vector
generate_compliant_vector <- function(incident_matrix, mced_compliance_rate,
                                      mced_start_cycle, mced_end_cycle) {
  H_vector <- as.matrix(incident_matrix[, "H"])

  # Get the number of cycles
  n_cycles <- nrow(H_vector)

  # Create a compliant_vector matrix initialized to 0
  compliant_vector <- matrix(0, nrow = n_cycles, ncol = 1)

  # Set compliant_vector to H_vector * mced_compliance_rate only for the specified range
  for (i in seq_len(n_cycles)) {
    if (i >= mced_start_cycle && i <= mced_end_cycle) {
      compliant_vector[i] <- H_vector[i] * mced_compliance_rate
    } else {
      compliant_vector[i] <- 0
    }
  }

  return(list(H_vector = H_vector, compliant_vector = compliant_vector))
}


generate_scenarios_fast <- function(n) {
  scen <- do.call(
    cbind,
    lapply(seq_len(n), function(k) {
      rep(rep(0:1, each = 2^(n - k)), times = 2^(k - 1))
    })
  )
  
  # Ensure matrix structure even when n = 1 (scen becomes a vector otherwise)
  if (n == 1 && is.null(dim(scen))) {
    scen <- matrix(scen, ncol = 1)
  }
  
  return(scen)
}

build_mced_cache <- function(specificity_df, cancer_list, K) {
  sub_df   <- specificity_df[specificity_df$Cancer %in% cancer_list, ]
  sub_df   <- sub_df[order(sub_df$Specificity), ]
  cancers  <- sub_df$Cancer
  
  if (length(cancer_list) == 0) {
    return(list(
      scen = matrix(nrow = 1, ncol = 0),
      scen_prob = 1,
      TN_mat = matrix(nrow = 1, ncol = 0),
      FP_mat = matrix(nrow = 1, ncol = 0),
      FU_mat = matrix(nrow = 1, ncol = 0),
      cancer_order = character(0)
    ))
  }
  
  
  spec     <- sub_df$Specificity
  
  n        <- length(cancers)
  
  ## 2.1  scenario lattice ----------------------------------------------------
  scen            <- generate_scenarios_fast(n)         # 2^n × n  (0 = TN, 1 = FP)
  is_FP           <- as.matrix(scen == 1)
  is_TN           <- !is_FP
  
  log_fp <-  matrix(log1p(-spec), nrow = nrow(is_FP), ncol = n, byrow = TRUE)
  log_tn <-  matrix(log(spec),    nrow = nrow(is_FP), ncol = n, byrow = TRUE)
  

  logP   <-  is_FP * log_fp + is_TN * log_tn
  scen_prob <- exp(rowSums(logP))
  stopifnot(abs(sum(scen_prob) - 1) < 1e-12)   # sanity check
  # 
  # print(head(scen))
  # print(head(is_FP))
  # ## 2.2  scenario probabilities (vectorised) ---------------------------------
  # logP  <- is_FP  * log1p(-spec) +            #  log(1 − specificity)
  #   is_TN  * log(spec)
  # 
  # print(head(logP))
  # logP_mat <- as.matrix(logP)
  # 
  # # now rowSums() will work even if there's only one column
  # scen_prob <- exp(rowSums(logP_mat))       # length = 2^n
  # 
  # print(scen_prob[1:6])
  # 
  # print(sum(scen_prob))
  
  ## 2.3  follow-up mask (first K FP per row) ---------------------------------
  if (K > 0) {
    max_K <- min(K, ncol(is_FP)) 
    
    if (n < 2) {
      cumsum_FP     <- apply(is_FP, 1, cumsum)
    } else {
      cumsum_FP     <- t(apply(is_FP, 1, cumsum))
    }
    
    scen_followup <- is_FP & (cumsum_FP <= max_K)
  } else {
    scen_followup <- matrix(FALSE, nrow(is_FP), ncol(is_FP))
  }
  
  
  return(list(
    cancer_order  = cancers,
    scen          = scen,            # keep for diagnostics
    scen_prob     = scen_prob,
    TN_mat        = is_TN,           # logical   2^n × n
    FP_mat        = is_FP,
    FU_mat        = scen_followup    # logical   2^n × n
  ))
}

run_mced_cached <- function(cache, H) {
  ## scale the one-time scenario probabilities by current cohort size
  scen_counts <- cache$scen_prob * H          # length = 2^n
  
  
  ## Fast aggregation: crossprod(logical->numeric , counts)  → n-vector
  total_TN <- as.numeric(crossprod(cache$TN_mat, scen_counts))
  total_FP <- as.numeric(crossprod(cache$FP_mat, scen_counts))
  total_FU <- as.numeric(crossprod(cache$FU_mat, scen_counts))
  
  list(
    total_TN         = total_TN,
    total_FP         = total_FP,
    total_followups  = total_FU,
    scenario_outcomes   = cache$scen,         # identical to old output
    scenario_counts     = scen_counts,
    scenario_followups  = cache$FU_mat,
    cancer_order        = cache$cancer_order
  )
}


# Run All Cycles
run_all_cycles <- function(compliant_vector,
                           H_vector,
                           female_prop,
                           fp_track_num,
                           specificity_df,
                           cancer_list,
                           female_cancers,
                           male_cancers,
                           final_result_list,
                           two_step,
                           overall_specificity, 
                           incident_matrix_trimmed) {

  non_gender_cancers <- setdiff(cancer_list, c(female_cancers, male_cancers))
  f_cancer_list      <- c(intersect(female_cancers, cancer_list), non_gender_cancers)
  m_cancer_list      <- c(intersect(male_cancers,   cancer_list), non_gender_cancers)

  # full_order <- specificity_df |>
  #   dplyr::filter(Cancer %in% cancer_list) |>
  #   dplyr::arrange(Specificity) |>
  #   dplyr::pull(Cancer)
  
  # # Build prevalence-sorted specificity order
  # prevalence_df <- incident_matrix_trimmed |>
  #   as.data.frame() |>
  #   dplyr::summarize(across(everything(), \(x) sum(x, na.rm = TRUE))) |>
  #   tidyr::pivot_longer(everything(), names_to = "Stage_Cancer", values_to = "Total") |>
  #   dplyr::mutate(
  #     Cancer = sub("^[^_]+_", "", Stage_Cancer)
  #   ) |>
  #   dplyr::group_by(Cancer) |>
  #   dplyr::summarize(Prevalence = sum(Total, na.rm = TRUE), .groups = "drop")
  # 
  # print(prevalence_df)
  # spec_rank_df <- specificity_df |>
  #   dplyr::filter(Cancer %in% cancer_list) |>
  #   dplyr::left_join(prevalence_df, by = "Cancer") |>
  #   dplyr::arrange(Specificity, desc(Prevalence))
  # 
  # full_order <- spec_rank_df |>
  #   dplyr::pull(Cancer)
  # 
  # print(full_order)
  full_order_all <- c(
    "breast",
    "prostate",
    "lung and bronchus",
    "colon and rectum",
    "lymphoma",
    "kidney and renal pelvis",
    "uterus",
    "pancreas",
    "bladder",
    "head and neck",
    "liver and intrahepatic bile duct",
    "stomach",
    "ovarian",
    "esophagus",
    "cervix",
    "anus"
  )
  
  full_order <- full_order_all[full_order_all %in% cancer_list]

  f_cache <- build_mced_cache(specificity_df, f_cancer_list, fp_track_num)
  m_cache <- build_mced_cache(specificity_df, m_cancer_list, fp_track_num)
  
  ## per‑sex matrix column names
  base_cols <- c("H", "compliant_pop",
                 as.vector(t(outer(full_order, c("_TN", "_FP", "_FP_fu"), paste0))),
                 "step1_FP", "step1_TN")          # NEW columns
  f_colnames <- m_colnames <- base_cols

  n_cycles   <- nrow(compliant_vector)
  f_mced_tn_fp <- matrix(0, n_cycles, length(f_colnames), dimnames = list(NULL, f_colnames))
  m_mced_tn_fp <- matrix(0, n_cycles, length(m_colnames), dimnames = list(NULL, m_colnames))

  ## overall matrix columns
  n_cancers <- length(cancer_list)
  overall_cols <- c("Total_FP", "Total_FP_fu",
                    paste0(1:n_cancers, "_FP"), "True_TN",
                    "FP_f", "FP_fu_f", "FP_m", "FP_fu_m",
                    "step1_FP_f", "step1_FP_m",
                    "step1_TN_f", "step1_TN_m")      # NEW

  f_overall <- m_overall <- matrix(0, n_cycles, length(overall_cols),
                                   dimnames = list(NULL, overall_cols))
  ## ── 1.  Cycle loop ──────────────────────────────────────────────────────────
  for (i in seq_len(n_cycles)) {

    total_H       <- H_vector[i]
    compliant_pop <- compliant_vector[i, 1]

    ## initialise step‑1 counters to 0 for single‑step runs
    step1_FP_total <- step1_TN_total <- 0

    if (two_step) {
      step1_FP_total <- compliant_pop * (1 - overall_specificity)
      step1_TN_total <- compliant_pop - step1_FP_total
      step2_pop_f    <- step1_FP_total * female_prop
      step2_pop_m    <- step1_FP_total * (1 - female_prop)
    } else {
      step2_pop_f    <- compliant_pop * female_prop
      step2_pop_m    <- compliant_pop * (1 - female_prop)
    }

    f_res <- run_mced_cached(f_cache, step2_pop_f)
    m_res <- run_mced_cached(m_cache, step2_pop_m)
    
    ## run per‑cancer FP engine (unchanged)
    # f_res <- simulate_mced(specificity_df, f_cancer_list, H = step2_pop_f, K = fp_track_num)
    # m_res <- simulate_mced(specificity_df, m_cancer_list, H = step2_pop_m, K = fp_track_num)

    ## ---------- populate per‑sex matrices ------------------------------------
    for (sex in c("f","m")) {

      mat   <- if (sex=="f") f_mced_tn_fp else m_mced_tn_fp
      res   <- if (sex=="f") f_res        else m_res
      prop  <- if (sex=="f") female_prop  else (1 - female_prop)
      idx   <- if (sex=="f") "f" else "m"

      mat[i,"H"]             <- total_H
      mat[i,"compliant_pop"] <- compliant_pop

      map <- match(full_order, res$cancer_order)
      TN  <- ifelse(is.na(map), 0, res$total_TN[map])
      FP  <- ifelse(is.na(map), 0, res$total_FP[map])
      FU  <- ifelse(is.na(map), 0, res$total_followups[map])

      for (j in seq_along(full_order)) {
        c_ <- full_order[j]
        mat[i, paste0(c_,"_TN")]    <- TN[j]
        mat[i, paste0(c_,"_FP")]    <- FP[j]
        mat[i, paste0(c_,"_FP_fu")] <- FU[j]
      }

      mat[i,"step1_FP"] <- step1_FP_total * prop
      mat[i,"step1_TN"] <- step1_TN_total * prop

      ## write back
      if (sex=="f") f_mced_tn_fp <- mat else m_mced_tn_fp <- mat
    }
    
    safe_row_sum <- function(x) {
      if (is.null(x)) {
        return(NULL)  # Or rep(0, n_cycles) if you want explicit zeros
      } else if (is.vector(x)) {
        return(matrix(x, ncol = 1))  # coerce to matrix
      } else if (is.matrix(x) && ncol(x) == 1) {
        return(x[, 1])  # drop to vector
      } else {
        return(rowSums(x, na.rm = TRUE))
      }
    }
    

    ## ---------- scenario‑based summaries (unchanged) -------------------------
    f_fp_per_scen <- safe_row_sum(f_res$scenario_outcomes)
    m_fp_per_scen <- safe_row_sum(m_res$scenario_outcomes)

    f_True_TN     <- sum(f_res$scenario_counts[f_fp_per_scen == 0])
    m_True_TN     <- sum(m_res$scenario_counts[m_fp_per_scen == 0])

    f_dist <- sapply(1:n_cancers, \(k) sum(f_res$scenario_counts[f_fp_per_scen == k]))
    m_dist <- sapply(1:n_cancers, \(k) sum(m_res$scenario_counts[m_fp_per_scen == k]))

    ## ---------- totals for overall matrices ----------------------------------
    fp_cols_f  <- grep("_FP$", colnames(f_mced_tn_fp))
    fp_cols_m  <- grep("_FP$", colnames(m_mced_tn_fp))
    fu_cols_f  <- grep("_FP_fu$", colnames(f_mced_tn_fp))
    fu_cols_m  <- grep("_FP_fu$", colnames(m_mced_tn_fp))

    step2_FP_f <- sum(f_mced_tn_fp[i, fp_cols_f])
    step2_FP_m <- sum(m_mced_tn_fp[i, fp_cols_m])
    step2_FU_f <- sum(f_mced_tn_fp[i, fu_cols_f])
    step2_FU_m <- sum(m_mced_tn_fp[i, fu_cols_m])

    step1_FP_f <- step1_FP_total * female_prop
    step1_FP_m <- step1_FP_total * (1 - female_prop)
    step1_TN_f <- step1_TN_total * female_prop
    step1_TN_m <- step1_TN_total * (1 - female_prop)

    ## female row
    f_overall[i,] <- c(
      step1_FP_f + step2_FP_f,          # Total_FP
      step2_FU_f,                       # Total_FP_fu  (no FU on step‑1)
      f_dist,                           # 1_FP … n_FP
      f_True_TN,           # True_TN  (step‑1 TN + step‑2 TN)
      step2_FP_f, step2_FU_f,           # FP_f / FP_fu_f  (step‑2 only)
      0, 0,                             # FP_m / FP_fu_m
      step1_FP_f, 0,           # step1_FP_f / _m
      step1_TN_f, 0            # step1_TN_f / _m
    )

    ## male row
    m_overall[i,] <- c(
      step1_FP_m + step2_FP_m,
      step2_FU_m,
      m_dist,
      m_True_TN,
      0, 0,
      step2_FP_m, step2_FU_m,
      0, step1_FP_m,
      0, step1_TN_m
    )
  } # end cycle loop

  ## ── 2.  Aggregate & store (existing logic, untouched) ──────────────────────



  # Initialize aggregated results for MCED by cancer
  aggregated_mced_tn_fp <- array(0, dim = c(n_cycles, 11, length(cancer_list)),
                                 dimnames = list(NULL,
                                                 c("H", "compliant_pop", "TN", "FP", "FP_fu", "TN_f",
                                                   "FP_f", "FP_fu_f", "TN_m", "FP_m", "FP_fu_m"),
                                                 cancer_list))

  # For each cycle and each cancer, aggregate:
  for (c_idx in seq_along(full_order)) {
    cancer <- full_order[c_idx]
    if (cancer %in% cancer_list) {
      # Extract female and male columns for this cancer
      f_TN_col <- paste0(cancer, "_TN")
      f_FP_col <- paste0(cancer, "_FP")
      f_FP_fu_col <- paste0(cancer, "_FP_fu")

      m_TN_col <- f_TN_col
      m_FP_col <- f_FP_col
      m_FP_fu_col <- f_FP_fu_col


      aggregated_mced_tn_fp[, "H", cancer] <- f_mced_tn_fp[, "H"]
      aggregated_mced_tn_fp[, "compliant_pop", cancer] <- f_mced_tn_fp[, "compliant_pop"]  # same as m_mced_tn_fp
      aggregated_mced_tn_fp[, "TN", cancer] <- f_mced_tn_fp[, f_TN_col] + m_mced_tn_fp[, m_TN_col]
      aggregated_mced_tn_fp[, "FP", cancer] <- f_mced_tn_fp[, f_FP_col] + m_mced_tn_fp[, m_FP_col]
      aggregated_mced_tn_fp[, "FP_fu", cancer] <- f_mced_tn_fp[, f_FP_fu_col] + m_mced_tn_fp[, m_FP_fu_col]
      aggregated_mced_tn_fp[, "TN_f", cancer] <- f_mced_tn_fp[, f_TN_col]
      aggregated_mced_tn_fp[, "FP_f", cancer] <- f_mced_tn_fp[, f_FP_col]
      aggregated_mced_tn_fp[, "FP_fu_f", cancer] <- f_mced_tn_fp[, f_FP_fu_col]
      aggregated_mced_tn_fp[, "TN_m", cancer] <- m_mced_tn_fp[, m_TN_col]
      aggregated_mced_tn_fp[, "FP_m", cancer] <- m_mced_tn_fp[, m_FP_col]
      aggregated_mced_tn_fp[, "FP_fu_m", cancer] <- m_mced_tn_fp[, m_FP_fu_col]
    }
  }

  # Store aggregated results into final_result_list for each cancer
  for (cancer in cancer_list) {
    final_result_list[["MCED"]][[cancer]][["mced_TN_FP"]] <- as.data.frame(aggregated_mced_tn_fp[, , cancer])
  }

  # Aggregate overall results:
  aggregated_overall <- f_overall + m_overall

  compliant_pop_vec <- as.vector(compliant_vector[,1])
  female_compliant_pop_vec <- compliant_pop_vec * female_prop
  male_compliant_pop_vec   <- compliant_pop_vec * (1 - female_prop)

  # Determine if all cancers in cancer_list belong to one gender category
  if (all(cancer_list %in% male_cancers)) {
    compliant_pop_vec <- male_compliant_pop_vec  # Use male-only population
    female_prop <- 0
  } else if (all(cancer_list %in% female_cancers)) {
    compliant_pop_vec <- female_compliant_pop_vec  # Use female-only population
    female_prop <- 1
  } else {
    compliant_pop_vec <- compliant_pop_vec  # Default: use total compliant population
  }
  # Convert aggregated_overall to a data.frame so we can add columns easily
  aggregated_overall_df <- as.data.frame(aggregated_overall)
  aggregated_overall_df$H              <- as.vector(H_vector[,1])
  aggregated_overall_df$compliant_pop  <- compliant_pop_vec
  aggregated_overall_df$non_compliant_pop <- aggregated_overall_df$H - aggregated_overall_df$compliant_pop
  aggregated_overall_df$female_pop     <- aggregated_overall_df$compliant_pop * female_prop
  aggregated_overall_df$male_pop       <- aggregated_overall_df$compliant_pop * (1 - female_prop)
  aggregated_overall_df$FP_step1 <- aggregated_overall_df$step1_FP_f + aggregated_overall_df$step1_FP_m
  aggregated_overall_df$FP_step2 <- aggregated_overall_df$FP_f       + aggregated_overall_df$FP_m
  aggregated_overall_df$TN_step1 <- aggregated_overall_df$step1_TN_f + aggregated_overall_df$step1_TN_m
  aggregated_overall_df$TN_step2 <- aggregated_overall_df$True_TN
  


  # Add the new columns:
  aggregated_overall_df$fp_counts <- rowSums(aggregated_overall_df[1:nrow(aggregated_overall_df), 3:(2 + length(cancer_list)), drop = FALSE])
  aggregated_overall_df$H <- as.vector(H_vector[ ,1])
  aggregated_overall_df$compliant_pop <- compliant_pop_vec
  aggregated_overall_df$non_compliant_pop <- aggregated_overall_df$H - aggregated_overall_df$compliant_pop

  final_result_list[["Overall"]][["mced_TN_FP"]] <- aggregated_overall_df


  # Return results
  # Gender-specific matrices and aggregated final_result_list are returned
  return(list(
    f_mced_tn_fp = f_mced_tn_fp,
    m_mced_tn_fp = m_mced_tn_fp,
    f_overall = f_overall,
    m_overall = m_overall,
    final_result_list = final_result_list
  ))
}
