calculate_union_probability <- function(compliance_probs) {
  # If there are no eligible cancers, return 0
  if (length(compliance_probs) == 0) {
    return(0)
  }
  
  n <- length(compliance_probs)
  total_prob <- 0
  
  # Loop over all subset sizes (1 to n) to apply gender mutual exclusivity
  for (k in 1:n) {
    combinations <- combn(compliance_probs, k, simplify = FALSE)
    
    sum_combinations <- sum(sapply(combinations, prod))
    
    total_prob <- total_prob + (-1)^(k + 1) * sum_combinations
  }
  
  return(total_prob)
}

calculate_screened_population <- function(final_result_list, incident_matrix, cancer_types, female_prop, screening_data, female_cancers, male_cancers) {
  cycle_results <- list(
    cycle = numeric(),
    eligible_cancers = character(),
    total_population = numeric(),
    total_elig = numeric(),
    screened_population = numeric(),
    no_screened_combined = numeric(),
    prop_elig = numeric(), 
    prop_screened = numeric(),
    total_female = numeric(),
    female_elig_1 = numeric(),
    screened_females_1 = numeric(),
    screened_females_all = numeric(),
    no_screened_females = numeric(),
    tn_females_1 = numeric(),
    tn_females_all = numeric(),
    tn_females_none = numeric(),
    fp_females_1 = numeric(),
    fp_females_all = numeric(),
    fp_females_none = numeric(),
    invasive_females_1 = numeric(),
    invasive_females_all = numeric(),
    invasive_females_none = numeric(),
    total_male = numeric(),
    male_elig_1 = numeric(),
    screened_males_1 = numeric(),
    screened_males_all = numeric(),
    no_screened_males = numeric(),
    tn_males_1 = numeric(),
    tn_males_all = numeric(),
    fp_males_1 = numeric(),
    fp_males_all = numeric(),
    fp_males_none = numeric(),
    invasive_males_1 = numeric(),
    invasive_males_all = numeric(),
    invasive_males_none = numeric()
  )
  
  # Loop over each cycle
  for (cycle in 1:nrow(incident_matrix)) {
    
    # total population from the incident_matrix
    total_population <- incident_matrix[cycle, "H"]
    
    # Split the population by gender
    female_population <- total_population * female_prop
    male_population <- total_population * (1 - female_prop)
    
    age <- cycle + 39
    
    # Filter eligible cancers for this cycle based on age
    eligible_cancers <- screening_data %>%
      filter(Cancer %in% cancer_types) %>% 
      filter(Start_age <= age & end_age >= age)
    
    # Calculate eligibility and compliance probabilities for each cancer group
    female_eligible <- eligible_cancers %>% filter(Cancer %in% female_cancers)
    male_eligible <- eligible_cancers %>% filter(Cancer %in% male_cancers)
    non_gender_eligible <- eligible_cancers %>% filter(!(Cancer %in% c(female_cancers, male_cancers)))
    
    # Calculate eligibility union probabilities for females and males
    female_elig_probs <- c(female_eligible$Eligibility, non_gender_eligible$Eligibility)
    male_elig_probs <- c(male_eligible$Eligibility, non_gender_eligible$Eligibility)
    
    # Probability that females/males are eligible for at least one cancer
    p_female_elig_1 <- calculate_union_probability(female_elig_probs)
    p_male_elig_1 <- calculate_union_probability(male_elig_probs)
    
    # Calculate counts for individuals eligible for at least one cancer
    female_elig_1 <- female_population * p_female_elig_1
    male_elig_1 <- male_population * p_male_elig_1
    
    # Calculate with compliant
    female_probs_compliant <- c(female_eligible$Eligibility * female_eligible$Compliance_frequency,
                                non_gender_eligible$Eligibility * non_gender_eligible$Compliance_frequency)
    p_fc_all <- if (prod(female_probs_compliant) != 1) prod(female_probs_compliant) else 0
    p_fc_union <- calculate_union_probability(female_probs_compliant)
    p_fc_none <- 1 - p_fc_union
    screened_females_1 <- female_population * p_fc_union
    screened_females_all <- female_population * p_fc_all
    no_screened_females <- female_population * p_fc_none
    
    male_probs_compliant <- c(male_eligible$Eligibility * male_eligible$Compliance_frequency,
                              non_gender_eligible$Eligibility * non_gender_eligible$Compliance_frequency)
    p_mc_all <- if (prod(male_probs_compliant) != 1) prod(male_probs_compliant) else 0
    p_mc_union <- calculate_union_probability(male_probs_compliant)
    p_mc_none <- 1 - p_mc_union
    screened_males_1 <- male_population * p_mc_union
    screened_males_all <- male_population * p_mc_all
    no_screened_males <- male_population * p_mc_none
    
    # Total screened population for at least one cancer occurrence
    screened_population <- screened_females_1 + screened_males_1
    total_elig <- female_elig_1 + male_elig_1
    prop_elig <- total_elig / total_population
    no_screened_combined <- no_screened_females + no_screened_males
    prop_screened <- screened_population / total_population
    
    # True Negatives with Specificity
    female_probs_tn <- c(female_eligible$Eligibility * female_eligible$Compliance_frequency * female_eligible$Specificity,
                         non_gender_eligible$Eligibility * non_gender_eligible$Compliance_frequency * non_gender_eligible$Specificity)
    p_fc_tn_all <- if (prod(female_probs_tn) != 1) prod(female_probs_tn) else 0
    p_fc_tn_union <- calculate_union_probability(female_probs_tn)
    tn_females_1 <- female_population * p_fc_tn_union
    tn_females_all <- female_population * p_fc_tn_all
    tn_females_none <- female_population * (1 - p_fc_tn_union)
    
    male_probs_tn <- c(male_eligible$Eligibility * male_eligible$Compliance_frequency * male_eligible$Specificity,
                       non_gender_eligible$Eligibility * non_gender_eligible$Compliance_frequency * non_gender_eligible$Specificity)
    p_mc_tn_all <- if (prod(male_probs_tn) != 1) prod(male_probs_tn) else 0
    p_mc_tn_union <- calculate_union_probability(male_probs_tn)
    tn_males_1 <- male_population * p_mc_tn_union
    tn_males_all <- male_population * p_mc_tn_all
    tn_males_none <- male_population * (1 - p_mc_tn_union)
    
    # False Positives with (1 - Specificity)
    female_probs_fp <- c(female_eligible$Eligibility * female_eligible$Compliance_frequency * (1 - female_eligible$Specificity),
                         non_gender_eligible$Eligibility * non_gender_eligible$Compliance_frequency * (1 - non_gender_eligible$Specificity))
    p_fc_fp_all <- if (prod(female_probs_fp) != 1) prod(female_probs_fp) else 0
    p_fc_fp_union <- calculate_union_probability(female_probs_fp)
    fp_females_1 <- female_population * p_fc_fp_union
    fp_females_all <- female_population * p_fc_fp_all
    fp_females_none <- female_population * (1 - p_fc_fp_union)
    
    male_probs_fp <- c(male_eligible$Eligibility * male_eligible$Compliance_frequency * (1 - male_eligible$Specificity),
                       non_gender_eligible$Eligibility * non_gender_eligible$Compliance_frequency * (1 - non_gender_eligible$Specificity))
    p_mc_fp_all <- if (prod(male_probs_fp) != 1) prod(male_probs_fp) else 0
    p_mc_fp_union <- calculate_union_probability(male_probs_fp)
    fp_males_1 <- male_population * p_mc_fp_union
    fp_males_all <- male_population * p_mc_fp_all
    fp_males_none <- male_population * (1 - p_mc_fp_union)
    
    # Invasive Procedures for False Positives
    female_probs_in <- c(female_eligible$Eligibility * female_eligible$Compliance_frequency * (1 - female_eligible$Specificity) * female_eligible$Invasive_procedure,
                         non_gender_eligible$Eligibility * non_gender_eligible$Compliance_frequency * (1 - non_gender_eligible$Specificity) * non_gender_eligible$Invasive_procedure)
    p_fc_fp_in_all <- if (prod(female_probs_in) != 1) prod(female_probs_in) else 0
    p_fc_fp_in_union <- calculate_union_probability(female_probs_in)
    in_females_1 <- female_population * p_fc_fp_in_union
    in_females_all <- female_population * p_fc_fp_in_all
    in_females_none <- female_population * (1 - p_fc_fp_in_union)
    
    male_probs_in <- c(male_eligible$Eligibility * male_eligible$Compliance_frequency * (1 - male_eligible$Specificity) * male_eligible$Invasive_procedure,
                       non_gender_eligible$Eligibility * non_gender_eligible$Compliance_frequency * (1 - non_gender_eligible$Specificity) * non_gender_eligible$Invasive_procedure)
    p_mc_fp_in_all <- if (prod(male_probs_in) != 1) prod(male_probs_in) else 0
    p_mc_fp_in_union <- calculate_union_probability(male_probs_in)
    in_males_1 <- male_population * p_mc_fp_in_union
    in_males_all <- male_population * p_mc_fp_in_all
    in_males_none <- male_population * (1 - p_mc_fp_in_union)
    
    # Append results to the list for this cycle
    cycle_results$cycle <- c(cycle_results$cycle, cycle)
    cycle_results$eligible_cancers <- c(cycle_results$eligible_cancers, paste(eligible_cancers$Cancer, collapse = ", "))
    cycle_results$screened_population <- c(cycle_results$screened_population, screened_population)
    cycle_results$no_screened_combined <- c(cycle_results$no_screened_combined, no_screened_combined)
    cycle_results$prop_screened <- c(cycle_results$prop_screened, prop_screened)
    cycle_results$total_population <- c(cycle_results$total_population, total_population)
    cycle_results$total_female <- c(cycle_results$total_female, female_population)
    cycle_results$total_male <- c(cycle_results$total_male, male_population)
    cycle_results$screened_females_1 <- c(cycle_results$screened_females_1, screened_females_1)
    cycle_results$screened_females_all <- c(cycle_results$screened_females_all, screened_females_all)
    cycle_results$no_screened_females <- c(cycle_results$no_screened_females, no_screened_females)
    cycle_results$screened_males_1 <- c(cycle_results$screened_males_1, screened_males_1)
    cycle_results$screened_males_all <- c(cycle_results$screened_males_all, screened_males_all)
    cycle_results$no_screened_males <- c(cycle_results$no_screened_males, no_screened_males)
    cycle_results$female_elig_1 <- c(cycle_results$female_elig_1, female_elig_1)
    cycle_results$male_elig_1 <- c(cycle_results$male_elig_1, male_elig_1)
    cycle_results$total_elig <- c(cycle_results$total_elig, total_elig)
    cycle_results$prop_elig <- c(cycle_results$prop_elig, prop_elig)
    cycle_results$tn_females_1 <- c(cycle_results$tn_females_1, tn_females_1)
    cycle_results$tn_females_all <- c(cycle_results$tn_females_all, tn_females_all)
    cycle_results$tn_females_none <- c(cycle_results$tn_females_none, tn_females_none)
    cycle_results$fp_females_1 <- c(cycle_results$fp_females_1, fp_females_1)
    cycle_results$fp_females_all <- c(cycle_results$fp_females_all, fp_females_all)
    cycle_results$fp_females_none <- c(cycle_results$fp_females_none, fp_females_none)
    cycle_results$invasive_females_1 <- c(cycle_results$invasive_females_1, in_females_1)
    cycle_results$invasive_females_all <- c(cycle_results$invasive_females_all, in_females_all)
    cycle_results$invasive_females_none <- c(cycle_results$invasive_females_none, in_females_none)
    cycle_results$tn_males_1 <- c(cycle_results$tn_males_1, tn_males_1)
    cycle_results$tn_males_all <- c(cycle_results$tn_males_all, tn_males_all)
    cycle_results$tn_males_none <- c(cycle_results$tn_males_none, tn_males_none)
    cycle_results$fp_males_1 <- c(cycle_results$fp_males_1, fp_males_1)
    cycle_results$fp_males_all <- c(cycle_results$fp_males_all, fp_males_all)
    cycle_results$fp_males_none <- c(cycle_results$fp_males_none, fp_males_none)
    cycle_results$invasive_males_1 <- c(cycle_results$invasive_males_1, in_males_1)
    cycle_results$invasive_males_all <- c(cycle_results$invasive_males_all, in_males_all)
    cycle_results$invasive_males_none <- c(cycle_results$invasive_males_none, in_males_none)
  }
  
  final_result_list$Overall$Usual_population_Screening_counts <- as.data.frame(cycle_results)
  return(final_result_list)
}