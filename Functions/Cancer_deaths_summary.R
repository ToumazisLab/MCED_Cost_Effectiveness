summarize_cancer_deaths_old <- function(mced_incident_matrix, survival_expanded) {
  # Extract stage and cancer from the matrix column names
  stage_cancer_split <- strsplit(colnames(mced_incident_matrix), "_")
  stages <- sapply(stage_cancer_split, function(x) x[1])
  cancers <- sapply(stage_cancer_split, function(x) paste(x[-1], collapse = "_"))
  
  # Initialize results list
  result_list <- list()
  
  # Loop over all matrix columns (Stage_Cancer combinations)
  for (col_idx in seq_along(colnames(mced_incident_matrix))) {
    stage <- stages[col_idx]
    cancer <- cancers[col_idx]
    
    # Get case counts per cycle
    cycles <- seq_len(nrow(mced_incident_matrix))
    ages <- cycles + 39
    cases <- mced_incident_matrix[, col_idx]
    
    # Filter survival data for this cancer and stage
    surv_data <- survival_expanded %>%
      filter(Cancer == cancer, Stage == stage)
    
    # Sum cancer deaths across all ages
    total_cancer_deaths <- 0
    
    for (i in seq_along(ages)) {
      age <- ages[i]
      n_cases <- cases[i]
      
      # Check if Cancer_death == 1 for this cancer, stage, and age
      match_row <- surv_data %>%
        filter(Age == age, Cancer_deaths == 1)
      
      if (nrow(match_row) > 0) {
        total_cancer_deaths <- total_cancer_deaths + n_cases
      }
    }
    
    # Store result
    result_list[[length(result_list) + 1]] <- data.frame(
      Cancer = cancer,
      Stage = stage,
      Cancer_dead_cases = total_cancer_deaths
    )
  }
  
  # Combine all results into one dataframe
  result_df <- bind_rows(result_list)
  
  return(result_df)
}

library(data.table)

summarize_cancer_deaths <- function(mced_incident_matrix, survival_expanded) {
  # Step 1: Parse Stage and Cancer from column names
  split_names <- tstrsplit(colnames(mced_incident_matrix), "_", fixed = TRUE)
  stages <- split_names[[1]]
  cancers <- do.call(paste, c(split_names[-1], sep = "_"))
  
  # Step 2: Melt the matrix into long format
  dt_cases <- as.data.table(mced_incident_matrix)
  dt_cases[, Age := .I + 39]  # add age column
  dt_cases_long <- melt(dt_cases, id.vars = "Age", variable.name = "Stage_Cancer", value.name = "Cases")
  
  # Add stage and cancer columns using pre-parsed names
  dt_cases_long[, Stage := stages[match(Stage_Cancer, colnames(mced_incident_matrix))]]
  dt_cases_long[, Cancer := cancers[match(Stage_Cancer, colnames(mced_incident_matrix))]]
  
  # Step 3: Filter survival_expanded for Cancer_deaths == 1
  dt_surv <- as.data.table(survival_expanded)
  dt_surv_deaths <- dt_surv[Cancer_deaths == 1, .(Cancer, Stage, Age)]
  
  # Step 4: Join cases with survival deaths
  setkey(dt_cases_long, Cancer, Stage, Age)
  setkey(dt_surv_deaths, Cancer, Stage, Age)
  dt_joined <- dt_cases_long[dt_surv_deaths, nomatch = 0]
  
  # Step 5: Aggregate total cancer deaths
  result_df <- dt_joined[, .(Cancer_dead_cases = sum(Cases)), by = .(Cancer, Stage)]
  return(as.data.frame(result_df))
}
