# Helper function to classify cases
classify_cases <- function(total_cases, cancer_type, cycle, screening_data, screening_cancers, case_type = "NA",
                           reduced_comply = NULL) {
  # Extract screening info
  info <- get_screening_info(cancer_type, screening_data, case_type, reduced_comply)
  start_age <- info$start_age
  end_age <- info$end_age
  eligibility <- info$eligibility
  compliance <- info$compliance
  is_screening_cancer <- cancer_type %in% screening_cancers
  
  
  # Age
  age <- cycle + 39
  
  # Initialize results
  result <- list(
    s_e_s = 0,   # Eligible & Compliant
    s_e_i = 0,   # Eligible & Non-compliant
    s_ne_i = 0   # Not Eligible / outside age range / no screening
  )
  
  if (total_cases == 0) {
    return(result)
  }
  
  if (!is_screening_cancer) {
    # Not a screening cancer
    result$s_ne_i <- total_cases
    return(result)
  }
  
  # Check age eligibility
  within_screening_age <- !is.na(start_age) && !is.na(end_age) && start_age != 0 && end_age != 0 && age >= start_age && age <= end_age
  
  if (!within_screening_age || is.na(eligibility) || is.na(compliance) || eligibility == 0) {
    # Outside screening age range or no eligibility
    result$s_ne_i <- total_cases
    return(result)
  }
  
  # Within screening age and eligible fraction
  total_e <- total_cases * eligibility
  total_ne <- total_cases - total_e   # not eligible fraction
  
  # Among eligible, compliance split
  total_e_c <- total_e * compliance       # Eligible & Compliant
  total_e_nc <- total_e - total_e_c       # Eligible & Non-compliant
  
  # Assign results
  result$s_e_s <- total_e_c
  result$s_e_i <- total_e_nc
  result$s_ne_i <- total_ne
  
  return(result)
}

create_cancer_df <- function(rows, arm = arm) {
  data.frame(
    Cycle = seq_len(rows),
    I_calculated = 0, II_calculated = 0, III_calculated = 0, IV_calculated = 0,
    s_e_m_I = 0, s_e_m_II = 0, s_e_m_III = 0, s_e_m_IV = 0,
    s_e_s_I = 0, s_e_s_II = 0, s_e_s_III = 0, s_e_s_IV = 0,
    s_e_i_I = 0, s_e_i_II = 0, s_e_i_III = 0, s_e_i_IV = 0,
    s_ne_i_I = 0, s_ne_i_II = 0, s_ne_i_III = 0, s_ne_i_IV = 0,
    
    # Totals
    Total_cases = 0,
    s_e_m_total = 0,
    s_e_s_total = 0,
    s_e_i_total = 0,
    s_ne_i_total = 0,
    
    # Derived
    m_dx_I = 0, m_dx_II = 0, m_dx_III = 0, m_dx_IV = 0, m_dx_total = 0,
    s_dx_I = 0, s_dx_II = 0, s_dx_III = 0, s_dx_IV = 0, s_dx_total = 0,
    i_dx_I = 0, i_dx_II = 0, i_dx_III = 0, i_dx_IV = 0, i_dx_total = 0
  )
}

get_screening_info <- function(cancer_type, screening_data, case_type, reduced_comply) {
  # Filter the screening data for the specified cancer type
  screening_info <- screening_data[screening_data$Cancer == cancer_type, ]
  
  # If no relevant data is found, return default values
  if (nrow(screening_info) == 0) {
    return(list(start_age = NA, end_age = NA, eligibility = FALSE, compliance = 0))
  }
  
  # Extract necessary details
  start_age <- screening_info$Start_age
  end_age <- screening_info$end_age
  eligibility <- screening_info$Eligibility
  compliance <- screening_info$Compliance_frequency
  
  # Apply the compliance rules
  if (case_type == "FN" && reduced_comply == 0) {
      compliance <- 0
    } else if (case_type == "FN" && reduced_comply != 0) {
      compliance <- compliance * (1/reduced_comply)
    }
  
  # Return the screening information as a list
  return(list(
    start_age = start_age,
    end_age = end_age,
    eligibility = eligibility,
    compliance = compliance
  ))
}
