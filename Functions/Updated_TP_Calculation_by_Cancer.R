calculate_transition_probabilities <- function(cancer_type, age, incidence_final, all_cause, female_prop, male_prop, female_cancers) {
  # Filter data for the specific cancer and age
  incidence_data <- subset(incidence_final, Cancer == cancer_type & Age == age)
  all_cause_data <- subset(all_cause, Age == age)
  
  # Initialize the transition matrix
  tp_matrix <- matrix(0, nrow = 6, ncol = 6)
  colnames(tp_matrix) <- rownames(tp_matrix) <- c("H", "I", "II", "III", "IV", "D")
  
  # # Adjust incidence rates based on gender-specific cancers
  # if (cancer_type %in% female_cancers) {
  #   incidence_data$Rate <- incidence_data$rate * female_prop
  # } else if (cancer_type == "prostate") {
  #   incidence_data$Rate <- incidence_data$rate * male_prop
  # }
  
  # Assign H_D first (all-cause death probability)
  if (nrow(all_cause_data) > 0) {
    tp_matrix["H", "D"] <- all_cause_data$Non_cancer_all_cause_probability
  } else {
    tp_matrix["H", "D"] <- 0
  }
  
  # Calculate remaining probability after H_D
  remaining_prob <- 1 - tp_matrix["H", "D"]
  
  # Assign H_I to H_IV directly from the remaining probability
  stages <- c("I", "II", "III", "IV")
  for (stage in stages) {
    stage_data <- incidence_data[incidence_data$Stage == stage, ]
    if (nrow(stage_data) > 0) {
      tp_matrix["H", stage] <- remaining_prob * stage_data$Rate
    }
  }
  
  # Assign H_H as the final remaining probability
  tp_matrix["H", "H"] <- 1 - sum(tp_matrix["H", ])
  
  # Set absorbing states
  for (stage in c("I", "II", "III", "IV", "D")) {
    tp_matrix[stage, stage] <- 1
  }
  
  return(tp_matrix)
}

# Generate the transition probability matrices
generate_tp_matrices <- function(cancers, incidence_final, all_cause, female_prop, male_prop, female_cancers) {
  age_groups <- 40:100
  tp_matrices <- list()
  
  for (cancer in cancers) {
    for (age in age_groups) {
      tp_matrix <- calculate_transition_probabilities(cancer, age, incidence_final, all_cause, female_prop, male_prop, female_cancers)
      matrix_name <- paste(age, cancer, sep = ".")
      tp_matrices[[matrix_name]] <- tp_matrix
    }
  }
  
  return(tp_matrices)
}
