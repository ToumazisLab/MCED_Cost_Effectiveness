calculate_overall_eligibility <- function(cancer_list, 
                                          cycle, 
                                          screening_data, 
                                          H, 
                                          female_prop) {
  age <- cycle + 39
  screening_cancers <- c("lung and bronchus", "breast", "cervix", "prostate", "colon and rectum")
  
  female_pop <- H * female_prop
  male_pop   <- H * (1 - female_prop)
  
  # Get only included screening cancers
  screening_info <- screening_data[
    screening_data$Cancer %in% cancer_list &
      screening_data$Cancer %in% screening_cancers, 
  ]
  
  # Keep only cancers in-range for this age
  eligible_cancers <- screening_info[
    screening_info$Start_age <= age & screening_info$end_age >= age, 
  ]
  
  if (nrow(eligible_cancers) == 0) {
    return(0)
  }
  
  gender_flags <- eligible_cancers$gender  # 0 = both, 1 = female, 2 = male
  
  if (any(gender_flags == 0)) {
    return(H)  # Any non-gendered cancer → whole population is eligible
  }
  
  if (any(gender_flags == 1) && any(gender_flags == 2)) {
    return(H)  # Both male-only and female-only cancers → whole population eligible
  }
  
  if (all(gender_flags == 1)) {
    return(female_pop)
  }
  
  if (all(gender_flags == 2)) {
    return(male_pop)
  }
  
  return(0)
}




usual_calculate_tn_fp <- function(final_result_list, 
                                  cancer_types, 
                                  incident_matrix, 
                                  screening_data, 
                                  female_cancers, 
                                  male_cancers, 
                                  female_prop, 
                                  arm = "usual", 
                                  tn_skip_years = 0) {
  # Number of cycles
  num_cycles <- nrow(incident_matrix)
  tn_skip <- tn_skip_years + 1
  
  # Prepare an Overall results data frame, by cycle.
  # We'll fill it by summing across all cancers later.
  overall_cols <- c(
    "Cycle", "Age", "H", 
    "Male_pop", "Female_pop", 
    # We'll store placeholders for columns we want to sum across
    "Eligible", "Non_Eligible",
    "Compliant", "Non_Compliant",
    "TN_f", "TN_m", "FP_f", "FP_m",
    "Invasive_f", "Invasive_m",
    "TN", "FP", 
    "Invasive"
  )
  
  overall_df <- data.frame(
    matrix(0, nrow = num_cycles, ncol = length(overall_cols)),
    stringsAsFactors = FALSE
  )
  names(overall_df) <- overall_cols
  
  # Fill the cycle-based columns in overall_df
  overall_df$Cycle <- 1:num_cycles
  overall_df$Age   <- overall_df$Cycle + 39
  overall_df$H     <- incident_matrix[, "H"]
  # Precompute total female vs. male population each cycle
  overall_df$Female_pop <- overall_df$H * female_prop
  overall_df$Male_pop   <- overall_df$H * (1 - female_prop)
  
  # For each cancer, initialize FP_TN_Screening_Summary with the new columns
  for (cancer in cancer_types) {
    final_result_list[[arm]][[cancer]]$FP_TN_Screening_Summary <- data.frame(
      Cycle           = 1:num_cycles,
      Age             = 1:num_cycles + 39,
      H               = incident_matrix[, "H"],
      Male_pop        = incident_matrix[, "H"] * (1 - female_prop),
      Female_pop      = incident_matrix[, "H"] * female_prop,
      Gender_eligible = numeric(num_cycles),  # depends on female-only, male-only, or both
      eligibility_rate= numeric(num_cycles),
      Eligible        = numeric(num_cycles),
      Non_Eligible    = numeric(num_cycles),
      compliance_rate = numeric(num_cycles),
      Compliant       = numeric(num_cycles),
      Non_Compliant   = numeric(num_cycles),
      TN_f            = numeric(num_cycles),
      TN_m            = numeric(num_cycles),
      FP_f            = numeric(num_cycles),
      FP_m            = numeric(num_cycles),
      Invasive_f      = numeric(num_cycles),
      Invasive_m      = numeric(num_cycles),
      TN              = numeric(num_cycles),
      FP              = numeric(num_cycles),
      invasive_rate   = numeric(num_cycles),
      Invasive        = numeric(num_cycles),
      stringsAsFactors = FALSE
    )
  }
  
  # Iterate through cycles to compute for each cancer
  for (cycle in 1:num_cycles) {
    # Overall population and male/female subpop
    H <- incident_matrix[cycle, "H"]
    female_pop <- H * female_prop
    male_pop   <- H * (1 - female_prop)
    age        <- cycle + 39
    
    # Temporary accumulators for the overall sum across all cancers this cycle
    # We won't accumulate columns we do not sum.
    cycle_sums <- setNames(as.list(rep(0, length(overall_cols))), overall_cols)
    
    cycle_sums$Cycle <- cycle
    cycle_sums$Age   <- age
    cycle_sums$H     <- H
    cycle_sums$Female_pop <- female_pop
    cycle_sums$Male_pop   <- male_pop
    
    for (cancer in cancer_types) {
      screening_info <- screening_data[screening_data$Cancer == cancer, ]
      
      # If the age is not in the screening range, skip
      if (age < screening_info$Start_age || age > screening_info$end_age) {
        # Just ensure the row is populated in the cancer’s data frame (mostly zeros)
        # We already have the row in place, but we can confirm or overwrite with zeros if needed.
        # For clarity, we do so here:
        final_result_list[[arm]][[cancer]]$FP_TN_Screening_Summary[cycle, 
                                                                   c("Gender_eligible","eligibility_rate","Eligible","Non_Eligible",
                                                                     "compliance_rate","Compliant","Non_Compliant",
                                                                     "TN_f","TN_m","FP_f","FP_m","Invasive_f","Invasive_m",
                                                                     "TN","FP","invasive_rate","Invasive")] <- 0
        next
      }
      
      # We are within the age range for screening:
      gender_flag <- screening_info$gender   # 0=both, 1=female, 2=male
      eligibility_rate <- screening_info$Eligibility
      compliance_rate  <- screening_info$Compliance_frequency
      specificity      <- screening_info$Specificity
      invasive_rate    <- screening_info$Invasive_procedure
      
      # -- Determine how many are "gender eligible" --
      #    That is the portion of the population relevant to this cancer
      if (gender_flag == 1) {
        # female-only
        gender_eligible_pop <- female_pop
      } else if (gender_flag == 2) {
        # male-only
        gender_eligible_pop <- male_pop
      } else {
        # 0 => non-gendered => entire population
        gender_eligible_pop <- H
      }
      
      # -- Eligible = fraction (eligibility_rate) of gender-eligible pop
      #    We'll split out how many are actually female vs. male
      #    for non-gendered cancers. For gendered, one side is 0.
      
      # female_eligible
      if (gender_flag == 1) {
        female_eligible <- female_pop * eligibility_rate
      } else if (gender_flag == 0) {
        female_eligible <- female_pop * eligibility_rate
      } else {
        female_eligible <- 0
      }
      
      # male_eligible
      if (gender_flag == 2) {
        male_eligible <- male_pop * eligibility_rate
      } else if (gender_flag == 0) {
        male_eligible <- male_pop * eligibility_rate
      } else {
        male_eligible <- 0
      }
      
      total_eligible   <- female_eligible + male_eligible
      non_eligible     <- H - total_eligible
      
      # -- Compliant = fraction (compliance_rate) of the eligible
      female_compliant <- female_eligible * compliance_rate
      male_compliant   <- male_eligible   * compliance_rate
      total_compliant  <- female_compliant + male_compliant
      non_compliant    <- total_eligible - total_compliant
      
      # -- TN and FP by gender
      female_TN <- female_compliant * specificity
      male_TN   <- male_compliant   * specificity
      female_FP <- female_compliant * (1 - specificity)
      male_FP   <- male_compliant   * (1 - specificity)
      
      total_TN <- female_TN + male_TN
      total_FP <- female_FP + male_FP
      
      # -- Invasive
      female_inv <- female_FP * invasive_rate
      male_inv   <- male_FP   * invasive_rate
      total_inv  <- female_inv + male_inv
      
      # Populate the row in the cancer’s FP_TN_Screening_Summary
      final_result_list[[arm]][[cancer]]$FP_TN_Screening_Summary[cycle, ] <- list(
        Cycle           = cycle,
        Age             = age,
        H               = H,
        Male_pop        = male_pop,
        Female_pop      = female_pop,
        Gender_eligible = gender_eligible_pop,
        eligibility_rate= eligibility_rate,
        Eligible        = total_eligible,
        Non_Eligible    = non_eligible,
        compliance_rate = compliance_rate,
        Compliant       = total_compliant,
        Non_Compliant   = non_compliant,
        TN_f            = female_TN,
        TN_m            = male_TN,
        FP_f            = female_FP,
        FP_m            = male_FP,
        Invasive_f      = female_inv,
        Invasive_m      = male_inv,
        TN              = total_TN,
        FP              = total_FP,
        invasive_rate   = invasive_rate,
        Invasive        = total_inv
      )
      
      # Accumulate into cycle_sums (for Overall)
      # We do NOT sum "Cycle", "Age", "H", "Male_pop", "Female_pop" repeatedly.

      cycle_sums$Compliant     <- cycle_sums$Compliant     + total_compliant
      cycle_sums$Non_Compliant <- cycle_sums$Non_Compliant + non_compliant
      cycle_sums$TN_f          <- cycle_sums$TN_f          + female_TN
      cycle_sums$TN_m          <- cycle_sums$TN_m          + male_TN
      cycle_sums$FP_f          <- cycle_sums$FP_f          + female_FP
      cycle_sums$FP_m          <- cycle_sums$FP_m          + male_FP
      cycle_sums$Invasive_f    <- cycle_sums$Invasive_f    + female_inv
      cycle_sums$Invasive_m    <- cycle_sums$Invasive_m    + male_inv
      cycle_sums$TN            <- cycle_sums$TN            + total_TN
      cycle_sums$FP            <- cycle_sums$FP            + total_FP
      cycle_sums$Invasive      <- cycle_sums$Invasive      + total_inv
    }  # end loop over cancers
    
    # Place these sums into overall_df
    # (We leave alone the columns we don’t sum repeatedly: Cycle, Age, H, Male_pop, Female_pop)
    overall_df[cycle, "Eligible"] <- calculate_overall_eligibility(
      cancer_list      = cancer_types,
      cycle            = cycle,
      screening_data   = screening_data,
      H                = H,
      female_prop      = female_prop
    )
    overall_df[cycle, "Non_Eligible"] <- H - overall_df$Eligible[cycle]
    overall_df[cycle, "Compliant"]     <- cycle_sums$Compliant
    overall_df[cycle, "Non_Compliant"] <- cycle_sums$Non_Compliant
    overall_df[cycle, "TN_f"]          <- cycle_sums$TN_f
    overall_df[cycle, "TN_m"]          <- cycle_sums$TN_m
    overall_df[cycle, "FP_f"]          <- cycle_sums$FP_f
    overall_df[cycle, "FP_m"]          <- cycle_sums$FP_m
    overall_df[cycle, "Invasive_f"]    <- cycle_sums$Invasive_f
    overall_df[cycle, "Invasive_m"]    <- cycle_sums$Invasive_m
    overall_df[cycle, "TN"]            <- cycle_sums$TN
    overall_df[cycle, "FP"]            <- cycle_sums$FP
    overall_df[cycle, "Invasive"]      <- cycle_sums$Invasive
  }  # end loop over cycles
  

  
  
  
  
  # Store the Overall summary
  final_result_list$Overall$FP_TN_Screening_Summary <- overall_df
  
  # Return updated final_result_list
  return(final_result_list)
}

