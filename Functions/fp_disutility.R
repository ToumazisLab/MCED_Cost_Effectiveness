apply_fp_disutility <- function(final_result_list, 
                                expanded_utility,
                                fp_count = 1, 
                                fp_disutility_val,
                                tn_disutility_val) {
  expanded_utility_ordered <- expanded_utility[order(expanded_utility$Age), ]
  
  fp_disutility <- fp_count * fp_disutility_val
  tn_disutility <- fp_count * tn_disutility_val
  
  overall_summary <- final_result_list[["Overall"]][["FP_TN_Screening_Summary"]]
  num_cycles <- nrow(overall_summary)
  
  new_cols <- c("f_utility", "m_utility", "total_utility",
                "total_fp_disutility", "total_invasive_disutility", "total_tn_disutility",
                "total_disutility", "adjusted_utility")
  
  for (col in new_cols) {
    overall_summary[[col]] <- 0
  }
  
  # Population-level utility for each cycle in Overall
  f_discounted_vec <- expanded_utility_ordered$DiscountedFemale[1:num_cycles]
  m_discounted_vec <- expanded_utility_ordered$DiscountedMale[1:num_cycles]
  
  for (i in seq_len(num_cycles)) {
    # Female & male population in this cycle
    f_pop <- overall_summary$Female_pop[i]
    m_pop <- overall_summary$Male_pop[i]
    
    # Utility = pop * discounted_utility
    f_util <- f_pop * f_discounted_vec[i]
    m_util <- m_pop * m_discounted_vec[i]
    total_util <- f_util + m_util
    
    overall_summary$f_utility[i]    <- f_util
    overall_summary$m_utility[i]    <- m_util
    overall_summary$total_utility[i]<- total_util
  }
  
  #    Then compute disutilities within each cancer’s data frame.
  #    Meanwhile, accumulate the female & male disutilities into Overall.
  
  for (cancer_name in names(final_result_list$usual)) {
    cancer_df <- final_result_list$usual[[cancer_name]]$FP_TN_Screening_Summary
    num_cycles_cancer <- nrow(cancer_df)
    
    # If mismatch in cycles, handle or stop:
    if (num_cycles_cancer != num_cycles) {
      stop("Cancer '", cancer_name, "' has a different number of cycles than Overall.")
    }
    
    for (col in new_cols) {
      cancer_df[[col]] <- 0
    }
    
    f_discounted_vec_cancer <- expanded_utility_ordered$DiscountedFemale[1:num_cycles_cancer]
    m_discounted_vec_cancer <- expanded_utility_ordered$DiscountedMale[1:num_cycles_cancer]
    
    # Compute utility & disutility per cycle
    for (i in seq_len(num_cycles_cancer)) {
      # (a) Utility
      f_pop <- cancer_df$Female_pop[i]
      m_pop <- cancer_df$Male_pop[i]
      
      f_util <- f_pop * f_discounted_vec_cancer[i]
      m_util <- m_pop * m_discounted_vec_cancer[i]
      total_util <- f_util + m_util
      
      cancer_df$f_utility[i]    <- f_util
      cancer_df$m_utility[i]    <- m_util
      cancer_df$total_utility[i]<- total_util
      
      # (b) Disutility:
      #     For false positives + invasive, use fp_disutility
      #     For true negatives, use tn_disutility
      f_fp_dis <- cancer_df$FP_f[i] * fp_disutility
      m_fp_dis <- cancer_df$FP_m[i] * fp_disutility
      
      f_invasive_dis <- cancer_df$Invasive_f[i] * fp_disutility
      m_invasive_dis <- cancer_df$Invasive_m[i] * fp_disutility

      f_tn_dis <- cancer_df$TN_f[i] * tn_disutility
      m_tn_dis <- cancer_df$TN_m[i] * tn_disutility
      
      
      total_fp_dis <- f_fp_dis + m_fp_dis
      total_invasive_dis <- f_invasive_dis + m_invasive_dis
      total_tn_dis <- f_tn_dis + m_tn_dis
      
      
      # (c) Adjusted utility for female & male
      cancer_df$total_fp_disutility[i] <- total_fp_dis
      cancer_df$total_invasive_disutility[i] <- total_invasive_dis
      cancer_df$total_tn_disutility[i] <- total_tn_dis
      
      
      # (d) Overall totals for each row
      cancer_df$total_disutility[i]    <- total_fp_dis + total_invasive_dis + total_tn_dis
      cancer_df$adjusted_utility[i]    <- total_util - cancer_df$total_disutility[i]
      
      # (e) Accumulate into Overall’s total_disutility_f and total_disutility_m
      overall_summary$total_fp_disutility[i] <- overall_summary$total_fp_disutility[i] + total_fp_dis
      overall_summary$total_invasive_disutility[i] <- overall_summary$total_invasive_disutility[i] + total_invasive_dis
      overall_summary$total_tn_disutility[i] <- overall_summary$total_tn_disutility[i] + total_tn_dis
    }
    
    final_result_list$usual[[cancer_name]]$FP_TN_Screening_Summary <- cancer_df
  }
  
  # Compute final adjusted utilities in Overall:
  for (i in seq_len(num_cycles)) {
    overall_summary$total_disutility[i] <- overall_summary$total_fp_disutility[i] +
      overall_summary$total_invasive_disutility[i] + overall_summary$total_tn_disutility[i]

    overall_summary$adjusted_utility[i] <- overall_summary$total_utility[i] -
      overall_summary$total_disutility[i]
  }
  
  final_result_list[["Overall"]][["FP_TN_Screening_Summary"]] <- overall_summary
  
  # 9. Return
  return(final_result_list)
}
