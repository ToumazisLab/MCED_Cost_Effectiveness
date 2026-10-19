decision_tree <- function(incident_matrix_trimmed, transposed_survival,
                               cost_per_year_2, discount_rate, cancer_list) {

  # --- Setup -----------------------------------------------------------------
  incident_df <- as.data.frame(incident_matrix_trimmed)
  cancer_pattern <- paste(cancer_list, collapse = "|")
  columns_to_include <- grepl(cancer_pattern, colnames(transposed_survival))
  survival_df <- as.data.frame(transposed_survival[, columns_to_include])

  common_columns <- intersect(colnames(incident_df), colnames(survival_df))
  incident_df <- incident_df[, common_columns, drop = FALSE]
  survival_df <- survival_df[, common_columns, drop = FALSE]
  survival_df[] <- lapply(survival_df, as.numeric)

  # Convert to matrices for faster indexing
  incident_mat <- as.matrix(incident_df)
  survival_mat <- as.matrix(survival_df)
  
  stage_names <- colnames(survival_mat)
  n_cycles <- nrow(incident_mat)
  n_years <- 5
  n_cols <- ncol(survival_mat)

  # --- Preprocess: Cost lookup list and discount matrix ---------------------
  cost_list <- cost_per_year_2 %>%
    mutate(Year = as.integer(gsub("Year ", "", Year))) %>%
    group_by(Cancer_Stage) %>%
    summarize(cost_vec = list(setNames(Cost, Year)), .groups = "drop") %>%
    deframe()

  discount_matrix <- outer(
    0:(n_cycles - 1), 0:(n_years - 1),
    function(cycle, year) 1 / ((1 + discount_rate) ^ (cycle + year))
  )

  # --- Initialize output matrix ---------------------------------------------
  total_costs <- matrix(0, nrow = n_cycles, ncol = n_cols)

  # --- Main computation loop ------------------------------------------------
  for (j in seq_len(n_cols)) {
    stage_type <- stage_names[j]
    cost_vec <- cost_list[[stage_type]]

    for (i in seq_len(n_cycles)) {
      surv_years <- survival_mat[i, j]

      if (!is.na(surv_years)) {
        full_years <- min(floor(surv_years), 5)
        partial_year <- surv_years - full_years

        # Get full year costs (only if year <= 5 and exists)
        full_costs <- unlist(cost_vec[as.character(seq_len(full_years))])
        
        if (length(full_costs) > 0) {
          discounted <- full_costs * discount_matrix[i, seq_along(full_costs)]
        } else {
          discounted <- 0
        }

        # Partial year (if survival not an integer and <5 years)
        if (partial_year > 0 && full_years < n_years) {
          part_cost <- cost_vec[[as.character(full_years + 1)]]
          if (!is.null(part_cost)) {
            disc_part <- (i - 1) + full_years
            discounted_partial <- part_cost * partial_year / ((1 + discount_rate) ^ disc_part)
          } else {
            discounted_partial <- 0
          }
        } else {
          discounted_partial <- 0
        }

        total_costs[i, j] <- sum(discounted) + discounted_partial
      }
    }
  }

  colnames(total_costs) <- stage_names

  # --- Combine with incidence ------------------------------------------------
  survival_calcs <- incident_mat * survival_mat
  total_costs_incidence <- total_costs * incident_mat

  total_costs_incidence_df <- as.data.frame(total_costs_incidence)
  survival_df <- as.data.frame(survival_calcs)
  
  colnames(total_costs_incidence_df) <- colnames(total_costs)
  colnames(survival_df) <- colnames(total_costs)
  
  return(list(
    total_costs_incidence = total_costs_incidence_df,
    survival_df           = survival_df
  ))
}


update_post_diagnosis_survival_deaths_utility <- function(final_result_list, 
                                                          u_post_diagnosis, 
                                                          u_survival_years, 
                                                          u_deaths, 
                                                          u_utility_cancers_wide, 
                                                          arm = "usual", group = "TP") {
  # Define stages and determine group suffix based on arm type
  stages <- c("I", "II", "III", "IV")
  group_suffix <- if (arm == "usual") "" else paste0("_", group)
  
  # Helper: Retrieve cancer screening dataframe for a specific cancer type
  get_cancer_screening <- function(final_result_list, arm, cancer_type) {
    if (arm == "MCED") {
      final_result_list[[arm]][[cancer_type]]$MCED_Cancer_cases
    } else if (arm == "usual") {
      final_result_list[[arm]][[cancer_type]]$Cancer_cases_screening
    } else {
      stop("Invalid arm specified. Expected 'MCED' or 'usual'.")
    }
  }
  
  # Helper: Update cancer screening dataframe in final_result_list for a given cancer type
  update_cancer_screening <- function(final_result_list, arm, cancer_type, screening_df) {
    if (arm == "MCED") {
      final_result_list[[arm]][[cancer_type]]$MCED_Cancer_cases <- screening_df
    } else if (arm == "usual") {
      final_result_list[[arm]][[cancer_type]]$Cancer_cases_screening <- screening_df
    } else {
      stop("Invalid arm specified. Expected 'MCED' or 'usual'.")
    }
    final_result_list
  }
  
  # ----- PART 1: Update cancer-specific screening data ----- #
  for (col_name in colnames(u_post_diagnosis)) {
    # Extract stage and cancer type from column name
    stage <- sub("_.*", "", col_name)
    cancer_type <- sub("^[^_]+_", "", col_name)
    
    # Get the corresponding cancer screening dataframe
    cancer_screening <- get_cancer_screening(final_result_list, arm, cancer_type)
    
    # Initialize required columns for every stage if they are missing
    for (s in stages) {
      required_cols <- paste0(s, c("_post_cost", "_survival_years", "_deaths", "_post_utility"), group_suffix)
      for (col in required_cols) {
        if (!col %in% colnames(cancer_screening)) {
          cancer_screening[[col]] <- 0
        }
      }
    }
    
    # Populate the screening dataframe with values from the corresponding matrices for each cycle
    for (cycle in seq_len(nrow(u_post_diagnosis))) {
      cancer_screening[[paste0(stage, "_post_cost", group_suffix)]][cycle] <-
        u_post_diagnosis[cycle, col_name]
      cancer_screening[[paste0(stage, "_survival_years", group_suffix)]][cycle] <-
        u_survival_years[cycle, col_name]
      cancer_screening[[paste0(stage, "_deaths", group_suffix)]][cycle] <-
        u_deaths[cycle, col_name]
      cancer_screening[[paste0(stage, "_post_utility", group_suffix)]][cycle] <-
        u_utility_cancers_wide[cycle, col_name]
    }
    
    # Save the updated screening dataframe back into final_result_list
    final_result_list <- update_cancer_screening(final_result_list, arm, cancer_type, cancer_screening)
  }
  
  # ----- PART 2: Compute Totals for Each Cancer Type ----- #
  for (cancer_type in names(final_result_list[[arm]])) {
    cancer_screening <- get_cancer_screening(final_result_list, arm, cancer_type)
    totals <- c("post_cost", "survival_years", "deaths", "post_utility")
    
    for (total in totals) {
      total_col <- paste0("Total_", total, group_suffix)
      stage_cols <- paste0(stages, "_", total, group_suffix)
      # Create total column if it does not exist
      if (!total_col %in% colnames(cancer_screening)) {
        cancer_screening[[total_col]] <- rowSums(cancer_screening[, stage_cols, drop = FALSE], na.rm = TRUE)
      }
    }
    
    final_result_list <- update_cancer_screening(final_result_list, arm, cancer_type, cancer_screening)
  }
  
  # ----- PART 3: Update Overall Summary ----- #
  overall_summary_name <- if (arm == "usual") "Cancer_cases_screening_summary" else "MCED_Cancer_costs_summary"
  overall_summary <- final_result_list[["Overall"]][[overall_summary_name]]
  
  # Determine all required columns from the cancer-specific screening data
  cancer_specific_cols <- unique(unlist(lapply(names(final_result_list[[arm]]), function(cancer_type) {
    screening <- get_cancer_screening(final_result_list, arm, cancer_type)
    colnames(screening)
  })))
  
  # Initialize missing columns in the overall summary
  for (col in cancer_specific_cols) {
    if (!col %in% colnames(overall_summary)) {
      overall_summary[[col]] <- 0
    }
  }
  
  # Aggregate each cancer-specific column into the overall summary by summing over cancers
  for (col in cancer_specific_cols) {
    overall_summary[[col]] <- rowSums(sapply(names(final_result_list[[arm]]), function(cancer_type) {
      screening <- get_cancer_screening(final_result_list, arm, cancer_type)
      if (col %in% colnames(screening)) screening[[col]] else 0
    }), na.rm = TRUE)
  }
  
  # Add totals over stages to the overall summary
  for (total in c("post_cost", "survival_years", "deaths", "post_utility")) {
    total_col <- paste0("Total_", total, group_suffix)
    stage_cols <- paste0(stages, "_", total, group_suffix)
    if (!total_col %in% colnames(overall_summary)) {
      overall_summary[[total_col]] <- rowSums(overall_summary[, stage_cols, drop = FALSE], na.rm = TRUE)
    }
  }
  
  # Update the overall summary in final_result_list
  final_result_list[["Overall"]][[overall_summary_name]] <- overall_summary
  
  return(final_result_list)
}





