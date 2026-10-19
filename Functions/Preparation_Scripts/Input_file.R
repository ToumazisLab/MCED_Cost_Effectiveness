update_utility <- function(d_e = 0.03, 
                           over_all_data_folder,
                           female_prop = def_param_list$female_prop) {
  age_utility <- read_csv(file.path(over_all_data_folder, "Age_adjusted_utility.csv"), show_col_types = FALSE)
  age_utility$utility_avg <- with(age_utility, female_prop * Female + (1 - female_prop) * Male)
  
  # Define age ranges
  age_ranges <- list(
    '40-49' = 40:49,
    '50-59' = 50:59,
    '60-69' = 60:69,
    '70-79' = 70:79,
    '80+'   = 80:100)
  
  # Expand the dataset by individual ages
  expanded_utility <- do.call(rbind, lapply(seq_len(nrow(age_utility)), function(i) {
    age_range <- age_ranges[[age_utility$Age[i]]]
    data.frame(
      Age       = age_range,
      Female    = rep(age_utility$Female[i], length(age_range)),
      Male      = rep(age_utility$Male[i], length(age_range)),
      utility_ave = rep(age_utility$utility_avg[i], length(age_range))
    )
  }))
  
  min_age <- min(expanded_utility$Age)
  expanded_utility$Time <- expanded_utility$Age - min_age
  
  # Apply the discount factor
  expanded_utility$DiscountedFemale <- expanded_utility$Female / ((1 + d_e) ^ expanded_utility$Time)
  expanded_utility$DiscountedMale   <- expanded_utility$Male   / ((1 + d_e) ^ expanded_utility$Time)
  expanded_utility$Discounted_avg   <- expanded_utility$utility_ave / ((1 + d_e) ^ expanded_utility$Time)
  
  
  return(expanded_utility)
}

processCancerData <- function(cancer_mort_seer, all_cause, cancer_list) {
  # Step 1: Filter rows where the 'Cancer' column is in the provided cancer_list
  filtered_cancer <- subset(cancer_mort_seer, Cancer %in% cancer_list)
  
  # Step 2: Aggregate the cancer mortality probabilities by Age
  total_cancer_mortality_by_age <- aggregate(probability ~ Age, 
                                             data = filtered_cancer, 
                                             sum)
  colnames(total_cancer_mortality_by_age)[2] <- "Total_cancer_probability"
  
  # Step 3: Merge the aggregated cancer mortality data with the all-cause mortality data
  all_cancer_mort <- merge(all_cause, total_cancer_mortality_by_age, 
                           by = "Age", all.x = TRUE)
  
  # Handle ages with no cancer data by replacing NA values with 0
  all_cancer_mort$Total_cancer_probability[is.na(all_cancer_mort$Total_cancer_probability)] <- 0
  
  # Calculate the non-cancer all-cause mortality probability
  all_cancer_mort$Non_cancer_all_cause_probability <- 
    all_cancer_mort$All_cause_mortality - all_cancer_mort$Total_cancer_probability
  
  # Return the resulting data frame
  return(all_cancer_mort)
}


load_data <- function(input_param_list,
                      def_param_list,
                      mda_data_folder = "Data/default9",
                      over_all_data_folder = "Data/overall",
                      raw_data_folder = "Data/raw",
                      cancer_list = cancer_types, 
                      cancer_types_MDA = cancer_types_MDA, 
                      cancer_types_overall = cancer_types,
                      run_type = "base") {
  
  # Always load raw mortality because other variables depend on it
  cancer_mort     <- read_csv(file.path(over_all_data_folder, "cancer_mort.csv"), show_col_types = FALSE)
  all_cause_raw   <- read_csv(file.path(over_all_data_folder, "all_cause_mortality_40.csv"), show_col_types = FALSE)
  all_cause       <- processCancerData(cancer_mort, all_cause_raw, cancer_list)
  
  # Conditionally load other elements
  if (run_type %in% c("base", "shiny")) {
    seer_draw               <- read_csv(file.path(over_all_data_folder, "seer_draw.csv"), show_col_types = FALSE)
    cost_per_year           <- read_csv(file.path(over_all_data_folder, "cost_per_year_2.csv"), show_col_types = FALSE)
    dwell_time              <- read_csv(file.path(over_all_data_folder, "dwell_time.csv"), show_col_types = FALSE)
    utility_long            <- read_csv(file.path(over_all_data_folder, "utility_long.csv"), show_col_types = FALSE)
    incidence_final         <- read_csv(file.path(over_all_data_folder, "incidence_final.csv"), show_col_types = FALSE)
    screening_fp_diag_costs <- read_csv(file.path(over_all_data_folder, "screening_fp_diag_costs_new.csv"), show_col_types = FALSE)
    survival_expanded       <- read_csv(file.path(over_all_data_folder, "survival_expanded_smooth.csv"), show_col_types = FALSE)
    transposed_survival     <- read_csv(file.path(over_all_data_folder, "transposed_survival_smoothed_4.14.25.csv"), show_col_types = FALSE)
    dwell_model_all_cancers <- read_csv(file.path(over_all_data_folder, "dwell_model_all_cancers.csv"), show_col_types = FALSE)
  }
  
  if (run_type == "shiny") {
    expanded_utility <- update_utility(0.03, over_all_data_folder, 
                                       female_prop = def_param_list$female_prop)
  } else {
    expanded_utility <- update_utility(input_param_list$d_e, over_all_data_folder, 
                                       female_prop = def_param_list$female_prop)
  }

  
  data_list <- list(
    seer_draw               = if (exists("seer_draw")) seer_draw else NULL,
    all_cancer_mort         = all_cause,
    cost_per_year           = if (exists("cost_per_year")) cost_per_year else NULL,
    dwell_time              = if (exists("dwell_time")) dwell_time else NULL,
    utility_long            = if (exists("utility_long")) utility_long else NULL,
    incidence_final         = if (exists("incidence_final")) incidence_final else NULL,
    screening_fp_diag_costs = if (exists("screening_fp_diag_costs")) screening_fp_diag_costs else NULL,
    survival_expanded       = if (exists("survival_expanded")) survival_expanded else NULL,
    transposed_survival     = if (exists("transposed_survival")) transposed_survival else NULL,
    expanded_utility        = expanded_utility,
    dwell_model_all_cancers = if (exists("dwell_model_all_cancers")) dwell_model_all_cancers else NULL
  )
  
  return(data_list)
}


