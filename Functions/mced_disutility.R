calculate_mced_utility <- function(final_result_list, expanded_utility, fp_disutility,
                                   female_cancers, male_cancers, female_prop) {
  
  # A helper function to ensure a column exists, else set to 0
  ensure_col <- function(df, colname) {
    if (!colname %in% colnames(df)) df[[colname]] <- 0
    df
  }
  
  # Process each cancer in MCED
  mced_cancers <- names(final_result_list[["MCED"]])
  
  for (cancer in mced_cancers) {
    mced_tn_fp_df <- final_result_list[["MCED"]][[cancer]][["mced_TN_FP"]]
    if (is.null(mced_tn_fp_df)) {
      stop("there is an issue with mced_tn_fp_df")
    }
    # Ensure needed columns exist, else 0
    mced_tn_fp_df <- ensure_col(mced_tn_fp_df, "H")
    mced_tn_fp_df <- ensure_col(mced_tn_fp_df, "FP_fu_compliant")
    mced_tn_fp_df <- ensure_col(mced_tn_fp_df, "FP_fu_compliant_f")
    mced_tn_fp_df <- ensure_col(mced_tn_fp_df, "FP_fu_compliant_m")
    mced_tn_fp_df <- ensure_col(mced_tn_fp_df, "m_invasive")
    mced_tn_fp_df <- ensure_col(mced_tn_fp_df, "u_FP")
    mced_tn_fp_df <- ensure_col(mced_tn_fp_df, "nc_FP")
    mced_tn_fp_df <- ensure_col(mced_tn_fp_df, "u_invasive")
    mced_tn_fp_df <- ensure_col(mced_tn_fp_df, "nc_Invasive")
    
    n_cycles <- nrow(mced_tn_fp_df)
    mced_tn_fp_df$Cycle <- seq_len(n_cycles)
    mced_tn_fp_df$Age <- mced_tn_fp_df$Cycle + 39
    
    # Merge with expanded_utility by Age
    merged_df <- merge(mced_tn_fp_df, expanded_utility, by = "Age", all.x = TRUE)
    # Handle any missing discount values
    merged_df$DiscountedFemale[is.na(merged_df$DiscountedFemale)] <- 0
    merged_df$DiscountedMale[is.na(merged_df$DiscountedMale)] <- 0
    merged_df$Discounted_avg[is.na(merged_df$Discounted_avg)] <- 0

    # Grab columns for calculations
    H <- merged_df$H
    FP_fu <- merged_df$FP_fu_compliant
    FP_fu_f <- merged_df$FP_fu_compliant_f
    FP_fu_m <- merged_df$FP_fu_compliant_m
    m_invasive <- merged_df$m_invasive
    Total_u_FP <- merged_df$u_FP
    Total_nc_FP <- merged_df$nc_FP
    Total_u_Invasive <- merged_df$u_invasive
    Total_nc_Invasive <- merged_df$nc_Invasive
    Total_usual_FP <- Total_u_FP + Total_nc_FP
    
    # Determine which utility columns to use based on cancer type
    if (cancer %in% female_cancers) {
      Female_population <- H * female_prop
      Total_population_utility <- Female_population * merged_df$DiscountedFemale
      screening_disutility <- FP_fu_f * fp_disutility
      diagnosis_disutility <- Total_usual_FP * fp_disutility
      m_invasive_disutility <- m_invasive * fp_disutility
      invasive_disutility <- Total_u_Invasive * fp_disutility
    } else if (cancer %in% male_cancers) {
      Male_population <- H * (1 - female_prop)
      Total_population_utility <- Male_population * merged_df$DiscountedMale
      screening_disutility <- FP_fu_m * fp_disutility
      diagnosis_disutility <- Total_usual_FP * fp_disutility
      m_invasive_disutility <- m_invasive * fp_disutility
      invasive_disutility <- Total_u_Invasive * fp_disutility
    } else {
      # Non-gender-specific cancer
      Total_population_utility <- H * merged_df$Discounted_avg
      screening_disutility <- FP_fu * fp_disutility
      diagnosis_disutility <- Total_usual_FP * fp_disutility
      m_invasive_disutility <- m_invasive * fp_disutility
      invasive_disutility <- Total_u_Invasive * fp_disutility
    }
    
    Total_disutility <- screening_disutility + diagnosis_disutility + invasive_disutility + m_invasive_disutility
    Adjusted_population_utility <- Total_population_utility - Total_disutility

    # Store results
    merged_df$Total_population_utility <- Total_population_utility
    merged_df$screening_disutility <- screening_disutility
    merged_df$diagnosis_disutility <- diagnosis_disutility
    merged_df$m_invasive_disutility <- m_invasive_disutility
    merged_df$invasive_disutility <- invasive_disutility
    merged_df$Total_disutility <- Total_disutility
    merged_df$Adjusted_population_utility <- Adjusted_population_utility
    # Drop columns not needed
    drop_cols <- c("DiscountedFemale", "DiscountedMale", "Discounted_avg",
                   "Female", "Male", "utility_ave", "Time")
    merged_df <- merged_df[ , !(colnames(merged_df) %in% drop_cols)]
    
    final_result_list[["MCED"]][[cancer]][["mced_TN_FP"]] <- merged_df
  }
  
  # Overall calculation:
  overall_df <- final_result_list[["Overall"]][["mced_TN_FP"]]
  if (is.null(overall_df)) {
    stop("You do not have mced_TN_FP for calculating mced_disutility.")
  }
  
  # Ensure needed columns exist
  overall_df <- ensure_col(overall_df, "H")
  overall_df <- ensure_col(overall_df, "FP_fu_compliant")
  overall_df <- ensure_col(overall_df, "m_invasive")
  overall_df <- ensure_col(overall_df, "u_FP")
  overall_df <- ensure_col(overall_df, "nc_FP")
  overall_df <- ensure_col(overall_df, "u_invasive")
  overall_df <- ensure_col(overall_df, "nc_Invasive")
  
  # Create Cycle if not present
  if (!"Cycle" %in% colnames(overall_df)) {
    overall_df$Cycle <- seq_len(nrow(overall_df))
  }
  overall_df$Age <- overall_df$Cycle + 39
  
  # Merge with expanded_utility by Age
  merged_overall <- merge(overall_df, expanded_utility, by = "Age", all.x = TRUE)
  merged_overall$DiscountedFemale[is.na(merged_overall$DiscountedFemale)] <- 0
  merged_overall$DiscountedMale[is.na(merged_overall$DiscountedMale)] <- 0
  merged_overall$Discounted_avg[is.na(merged_overall$Discounted_avg)] <- 0

  # Population utility uses Discounted_avg overall
  H_overall <- merged_overall$H
  Total_population_utility_overall <- H_overall * merged_overall$Discounted_avg
  
  # For disutilities, separate them the same way as above:
  # 1) screening disutility = MCED FP follow-ups
  Total_FP_fu_overall <- merged_overall$FP_fu_compliant
  Total_m_invasive_overall <- merged_overall$m_invasive
  
  # 2) diagnosis disutility = usual follow-ups (FP + non-compliance, etc.)
  Total_u_FP_overall <- merged_overall$u_FP
  Total_nc_FP_overall <- merged_overall$nc_FP
  Total_usual_FP_overall <- Total_u_FP_overall + Total_nc_FP_overall
  
  # 3) invasive disutility = sum of invasive procedures
  Total_u_Invasive_overall <- merged_overall$u_invasive
  Total_nc_Invasive_overall <- merged_overall$nc_Invasive
  Total_Invasive_overall <- Total_u_Invasive_overall + Total_nc_Invasive_overall
  
  # Calculate them individually
  screening_disutility_overall <- Total_FP_fu_overall * fp_disutility
  diagnosis_disutility_overall <- Total_usual_FP_overall * fp_disutility
  m_invasive_disutility_overall <- Total_m_invasive_overall * fp_disutility
  invasive_disutility_overall  <- Total_Invasive_overall * fp_disutility
  
  # Sum to get total
  Total_disutility_overall <- screening_disutility_overall +
    diagnosis_disutility_overall + m_invasive_disutility_overall +
    invasive_disutility_overall
  
  # Adjusted population utility
  Adjusted_population_utility_overall <- Total_population_utility_overall - Total_disutility_overall
  
  # Store in merged_overall
  merged_overall$Total_population_utility <- Total_population_utility_overall
  merged_overall$screening_disutility <- screening_disutility_overall
  merged_overall$diagnosis_disutility <- diagnosis_disutility_overall
  merged_overall$invasive_disutility <- invasive_disutility_overall
  merged_overall$m_invasive_disutility <- m_invasive_disutility_overall
  merged_overall$Total_disutility <- Total_disutility_overall
  merged_overall$Adjusted_population_utility <- Adjusted_population_utility_overall
  
  # Remove unneeded columns
  drop_cols <- c("DiscountedFemale", "DiscountedMale", "Discounted_avg",
                 "Female", "Male", "utility_ave", "Time")
  merged_overall <- merged_overall[ , !(colnames(merged_overall) %in% drop_cols)]
  
  final_result_list[["Overall"]][["mced_TN_FP"]] <- merged_overall
  
  return(final_result_list)
}
