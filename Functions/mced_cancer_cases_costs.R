assign_costs <- function(final_result_list, 
                         screening_data,
                         screening_cancers,
                         arm = "usual",
                         discount_factor = d_c) {
  
  # Helper: discount factor by cycle
  discount_fn <- function(cycle) {
    (1 / ((1 + discount_factor)^(cycle - 1)))
  }
  
  # Get list of cancers processed
  cancers <- setdiff(names(final_result_list[[arm]]), "Cancer_cases_screening_summary")
  cancers <- cancers[cancers != "Overall"] # exclude overall if it exists
  
  # Function to compute costs for a given cancer_type
  compute_cancer_costs <- function(df, cancer_type) {
    # Lookup costs in screening_data
    cost_info <- screening_data[screening_data$Cancer == cancer_type, ]
    if (nrow(cost_info) == 0) {
      # If no cost info found, default to zeros
      Screening_cost <- 0
      FP_cost <- 0
      incidental_cost <- 0
      invasive_cost <- 0
    } else {
      Screening_cost <- cost_info$Screening_cost
      FP_cost <- cost_info$FP_cost
      incidental_cost <- cost_info$incidental_cost
      invasive_cost <- cost_info$invasive_cost
    }
    
    n <- nrow(df)
    stages <- c("I", "II", "III", "IV")
    buckets <- c("s_e_s", "s_e_i", "s_ne_i")
    cost_types <- c("imaging", "diagnostic", "invasive")
    
    # Create a template of zero columns
    cost_cols <- c()
    for (ct in cost_types) {
      for (b in buckets) {
        for (st in stages) {
          cost_cols <- c(cost_cols, paste(ct, b, st, sep = "_"))
        }
      }
    }
    
    cost_df <- data.frame(Cycle = df$Cycle)
    for (colname in cost_cols) cost_df[[colname]] <- 0
    
    # Determine if this is a screening cancer
    is_screening_cancer <- cancer_type %in% screening_cancers
    
    # Cost lookup rules:
    # s_e_s: imaging = Screening_cost, diagnostic = FP_cost, invasive = invasive_cost
    # s_e_i: imaging = incidental_cost, diagnostic = FP_cost, invasive = invasive_cost
    # s_ne_i (screening cancer): imaging = incidental_cost, diagnostic = FP_cost, invasive = invasive_cost
    # s_ne_i (non-screening cancer): imaging = incidental_cost, diagnostic = 0, invasive = invasive_cost
    
    bucket_costs <- list(
      s_e_s = list(imaging = Screening_cost, diagnostic = FP_cost, invasive = invasive_cost),
      s_e_i = list(imaging = incidental_cost, diagnostic = FP_cost, invasive = invasive_cost)
    )
    if (is_screening_cancer) {
      bucket_costs$s_ne_i <- list(imaging = incidental_cost, diagnostic = FP_cost, invasive = invasive_cost)
    } else {
      bucket_costs$s_ne_i <- list(imaging = incidental_cost, diagnostic = 0, invasive = invasive_cost)
    }
    
    # Apply costs per cycle and stage
    for (i in 1:n) {
      d_factor <- discount_fn(i)
      for (b in buckets) {
        for (st in stages) {
          count_col <- paste(b, st, sep = "_")
          count_val <- df[i, count_col]
          if (count_val == 0) next
          
          # Apply costs
          for (ct in cost_types) {
            cost_val <- bucket_costs[[b]][[ct]] * count_val * d_factor
            cost_df[i, paste(ct, b, st, sep = "_")] <- cost_val
          }
        }
      }
    }
    
    # Totals across stages for each bucket and cost type
    for (ct in cost_types) {
      for (b in buckets) {
        st_cols <- paste(ct, b, stages, sep = "_")
        total_col <- paste(ct, b, "total", sep = "_")
        cost_df[[total_col]] <- rowSums(cost_df[st_cols])
      }
    }
    
    # s_dx = s_e_s
    # i_dx = s_e_i + s_ne_i
    for (ct in cost_types) {
      for (st in stages) {
        cost_df[[paste(ct, "s_dx", st, sep = "_")]] <- cost_df[[paste(ct, "s_e_s", st, sep = "_")]]
        cost_df[[paste(ct, "i_dx", st, sep = "_")]] <- cost_df[[paste(ct, "s_e_i", st, sep = "_")]] + 
          cost_df[[paste(ct, "s_ne_i", st, sep = "_")]]
      }
      # Totals across stages for s_dx and i_dx
      cost_df[[paste(ct, "s_dx_total", sep = "_")]] <- rowSums(cost_df[paste(ct, "s_dx", stages, sep = "_")])
      cost_df[[paste(ct, "i_dx_total", sep = "_")]] <- rowSums(cost_df[paste(ct, "i_dx", stages, sep = "_")])
    }
    
    # Total cost across cost types for s_dx and i_dx
    cost_df$total_imaging_per_cycle <- cost_df$imaging_s_dx_total + cost_df$imaging_i_dx_total
    cost_df$total_diagnostic_per_cycle <- cost_df$diagnostic_s_dx_total + cost_df$diagnostic_i_dx_total
    cost_df$total_invasive_per_cycle <- cost_df$invasive_s_dx_total + cost_df$invasive_i_dx_total
    
    cost_df$total_s_dx_per_cycle <- cost_df$imaging_s_dx_total + cost_df$diagnostic_s_dx_total + cost_df$invasive_s_dx_total
    cost_df$total_i_dx_per_cycle <- cost_df$imaging_i_dx_total + cost_df$diagnostic_i_dx_total + cost_df$invasive_i_dx_total
    cost_df$total_cost_per_cycle <- cost_df$total_s_dx_per_cycle + cost_df$total_i_dx_per_cycle
    
    cost_df
  }
  
  # Compute costs per cancer
  for (cancer_type in cancers) {
    if (!"Cancer_cases_screening" %in% names(final_result_list[[arm]][[cancer_type]])) next
    df <- final_result_list[[arm]][[cancer_type]]$Cancer_cases_screening
    cost_df <- compute_cancer_costs(df, cancer_type)
    final_result_list[[arm]][[cancer_type]]$Cancer_costs <- cost_df
  }
  
  # Compute overall summary
  # Sum across all cancers
  cost_dfs <- lapply(cancers, function(ct) final_result_list[[arm]][[ct]]$Cancer_costs)
  cost_dfs <- cost_dfs[!sapply(cost_dfs, is.null)]
  
  if (length(cost_dfs) > 0) {
    overall_cols <- colnames(cost_dfs[[1]])
    overall_cost_summary <- data.frame(Cycle = cost_dfs[[1]]$Cycle)
    for (colnm in overall_cols[-1]) {
      vals <- sapply(cost_dfs, function(x) if (colnm %in% colnames(x)) x[[colnm]] else 0)
      overall_cost_summary[[colnm]] <- rowSums(vals)
    }
    final_result_list[["Overall"]][["Cancer_costs_summary"]] <- overall_cost_summary
  } else {
    final_result_list[["Overall"]][["Cancer_costs_summary"]] <- NULL
  }
  
  return(final_result_list)
}

mced_assign_cost <- function(final_result_list, 
                             screening_data,
                             screening_cancers,
                             arm = "MCED",
                             discount_factor = d_c) {
  
  # Helper: discount factor by cycle
  discount_fn <- function(cycle) {
    (1 / ((1 + discount_factor)^(cycle - 1)))
  }
  
  # Identify cancers processed
  cancers <- names(final_result_list[[arm]])
  
  stages <- c("I", "II", "III", "IV")
  buckets <- c("s_e_m", "s_e_s", "s_e_i", "s_ne_i")
  cost_types <- c("imaging", "diagnostic", "invasive")
  
  # Create a zero-filled cost data frame template
  create_cost_df <- function(n) {
    df <- data.frame(Cycle = 1:n)
    for (ct in cost_types) {
      for (b in buckets) {
        for (st in stages) {
          df[[paste(ct, b, st, sep = "_")]] <- 0
        }
      }
    }
    df
  }
  
  # Assign costs given the category and whether it's a screening cancer
  get_bucket_costs <- function(is_screening_cancer, is_prostate, Screening_cost, FP_cost, incidental_cost, invasive_cost) {
    # s_e_m
    if (is_screening_cancer) {
      if (is_prostate) {
        s_e_m_costs <- list(imaging = 0, diagnostic = FP_cost, invasive = invasive_cost)
      } else {
      s_e_m_costs <- list(imaging = Screening_cost, diagnostic = FP_cost, invasive = invasive_cost) 
      }
    } else {
      s_e_m_costs <- list(imaging = incidental_cost, diagnostic = 0, invasive = invasive_cost)
    }
    
    # s_e_s
    s_e_s_costs <- list(imaging = Screening_cost, diagnostic = FP_cost, invasive = invasive_cost)
    
    # s_e_i
    s_e_i_costs <- list(imaging = incidental_cost, diagnostic = FP_cost, invasive = invasive_cost)
    
    # s_ne_i
    if (is_screening_cancer) {
      s_ne_i_costs <- list(imaging = incidental_cost, diagnostic = FP_cost, invasive = invasive_cost)
    } else {
      s_ne_i_costs <- list(imaging = incidental_cost, diagnostic = 0, invasive = invasive_cost)
    }
    
    list(s_e_m = s_e_m_costs,
         s_e_s = s_e_s_costs,
         s_e_i = s_e_i_costs,
         s_ne_i = s_ne_i_costs)
  }
  
  compute_cancer_costs <- function(df, cancer_type) {
    # df is Cancer_cases_screening for that cancer
    
    # Lookup costs in screening_data
    cost_info <- screening_data[screening_data$Cancer == cancer_type, ]
    if (nrow(cost_info) == 0) {
      Screening_cost <- 0
      FP_cost <- 0
      incidental_cost <- 0
      invasive_cost <- 0
    } else {
      Screening_cost <- cost_info$Screening_cost
      FP_cost <- cost_info$FP_cost
      incidental_cost <- cost_info$incidental_cost
      invasive_cost <- cost_info$invasive_cost
    }
    
    is_screening_cancer <- cancer_type %in% screening_cancers
    is_prostate <- cancer_type == "prostate"
    
    bucket_costs <- get_bucket_costs(is_screening_cancer, is_prostate, Screening_cost, FP_cost, incidental_cost, invasive_cost)
    
    n <- nrow(df)
    cost_df <- create_cost_df(n)
    
    # For each cycle and stage, apply costs
    for (i in 1:n) {
      d_factor <- discount_fn(i)
      for (b in buckets) {
        for (st in stages) {
          count_col <- paste(b, st, sep = "_")
          count_val <- df[i, count_col]
          if (count_val > 0) {
            # Apply costs
            for (ct in cost_types) {
              cost_val <- bucket_costs[[b]][[ct]] * count_val * d_factor
              cost_df[i, paste(ct, b, st, sep = "_")] <- cost_val
            }
          }
        }
      }
    }
    
    # Totals across stages for each bucket and cost type
    for (ct in cost_types) {
      for (b in buckets) {
        st_cols <- paste(ct, b, stages, sep = "_")
        cost_df[[paste(ct, b, "total", sep = "_")]] <- rowSums(cost_df[st_cols])
      }
    }
    
    # m_dx = s_e_m
    # s_dx = s_e_s
    # i_dx = s_e_i + s_ne_i
    # Add these by stage
    for (ct in cost_types) {
      for (st in stages) {
        cost_df[[paste(ct, "m_dx", st, sep = "_")]] <- cost_df[[paste(ct, "s_e_m", st, sep = "_")]]
        cost_df[[paste(ct, "s_dx", st, sep = "_")]] <- cost_df[[paste(ct, "s_e_s", st, sep = "_")]]
        cost_df[[paste(ct, "i_dx", st, sep = "_")]] <- cost_df[[paste(ct, "s_e_i", st, sep = "_")]] +
          cost_df[[paste(ct, "s_ne_i", st, sep = "_")]]
      }
      
      # Totals across stages for m_dx, s_dx, i_dx
      cost_df[[paste(ct, "m_dx_total", sep = "_")]] <- rowSums(cost_df[paste(ct, "m_dx", stages, sep = "_")])
      cost_df[[paste(ct, "s_dx_total", sep = "_")]] <- rowSums(cost_df[paste(ct, "s_dx", stages, sep = "_")])
      cost_df[[paste(ct, "i_dx_total", sep = "_")]] <- rowSums(cost_df[paste(ct, "i_dx", stages, sep = "_")])
    }
    
    # Total cost across categories (m_dx, s_dx, i_dx) for each cost type
    # imaging: imaging_m_dx_total + imaging_s_dx_total + imaging_i_dx_total
    cost_df$total_imaging_per_cycle <- cost_df$imaging_m_dx_total + cost_df$imaging_s_dx_total + cost_df$imaging_i_dx_total
    cost_df$total_diagnostic_per_cycle <- cost_df$diagnostic_m_dx_total + cost_df$diagnostic_s_dx_total + cost_df$diagnostic_i_dx_total
    cost_df$total_invasive_per_cycle <- cost_df$invasive_m_dx_total + cost_df$invasive_s_dx_total + cost_df$invasive_i_dx_total
    
    # Total cost across cost types for all categories
    cost_df$total_m_dx_per_cycle <- cost_df$imaging_m_dx_total + cost_df$diagnostic_m_dx_total + cost_df$invasive_m_dx_total
    cost_df$total_s_dx_per_cycle <- cost_df$imaging_s_dx_total + cost_df$diagnostic_s_dx_total + cost_df$invasive_s_dx_total
    cost_df$total_i_dx_per_cycle <- cost_df$imaging_i_dx_total + cost_df$diagnostic_i_dx_total + cost_df$invasive_i_dx_total
    
    # Overall total per cycle (sum of all cost types and categories)
    cost_df$total_cost_per_cycle <- cost_df$total_m_dx_per_cycle + cost_df$total_s_dx_per_cycle + cost_df$total_i_dx_per_cycle
    
    cost_df
  }
  
  # Compute costs per cancer
  cost_dfs <- list()
  for (cancer_type in cancers) {
    if (!"MCED_Cancer_cases" %in% names(final_result_list[[arm]][[cancer_type]])) next
    df <- final_result_list[[arm]][[cancer_type]]$MCED_Cancer_cases
    cost_df <- compute_cancer_costs(df, cancer_type)
    final_result_list[[arm]][[cancer_type]]$MCED_Cancer_costs <- cost_df
    cost_dfs[[cancer_type]] <- cost_df
  }
  
  # Compute overall summary
  if (length(cost_dfs) > 0) {
    # All have the same columns, sum across cancers
    overall_cols <- colnames(cost_dfs[[1]])
    overall_cost_summary <- data.frame(Cycle = cost_dfs[[1]]$Cycle)
    
    for (colnm in overall_cols[-1]) {
      vals <- sapply(cost_dfs, function(x) if (colnm %in% colnames(x)) x[[colnm]] else 0)
      overall_cost_summary[[colnm]] <- rowSums(vals)
    }
    final_result_list[["Overall"]][["MCED_Cancer_costs_summary"]] <- overall_cost_summary
  } else {
    final_result_list[["Overall"]][["MCED_Cancer_costs_summary"]] <- NULL
  }
  
  return(final_result_list)
}

