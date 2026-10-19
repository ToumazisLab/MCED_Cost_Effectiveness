process_sensitivity <- function(cancer_list, sensitivity_long, seer_draw, xStage, active_slip,
                                dwell_model_df) {
  dwell_standard_model <- dwell_model_df[dwell_model_df$Cancer %in% cancer_list, ]
    
  # Create synthetic `c` and `n` columns
  sensitivity_long$n <- 1000
  sensitivity_long <- sensitivity_long %>%
    mutate(original_sens = Sensitivity,
           c = Sensitivity * n)
  
  # Fix sensitivity with isotonic regression within each cancer type
  sensitivity_long <- sensitivity_long %>%
    group_by(Cancer) %>%
    mutate(sens = isotone_fix(original_sens, Stage, n)) %>%
    ungroup() %>%
    select(Cancer, Stage, c, n, original_sens, sens)
  
  # Join sensitivity_long with SEER data
  iso_sens_joined <- sensitivity_long
  seer_draw <- seer_draw[seer_draw$SEER_Draw %in% unique(sensitivity_long$Cancer), ]
  
  iso_join_seer <- seer_draw %>%
    left_join(sensitivity_long %>%
                mutate(SEER_Draw = Cancer)) %>%
    mutate(Cancer = replace_na(Cancer, "NotFound"),
           c = replace_na(c, 0),
           n = replace_na(n, 0),
           sens = replace_na(original_sens, 0.0),
           original_sens = replace_na(original_sens, 0.0)) %>%
    select(SEER_Draw, Stage, IR, Survival, c, n, original_sens, sens) %>%
    rename(Cancer = SEER_Draw)
  
  # Prepare data for intercept modeling
  incidence_sens_source <- iso_join_seer %>%
    filter(Stage != "NotStaged") %>%
    mutate(iso_sens = sens)
  
  incidence_excluded_source <- iso_join_seer %>%
    filter(Stage == "NotStaged")
  
  # Generate slip rates
  dwell_slip_rate <- exact_slip_rate_from_dwell(dwell_standard_model, screen_interval = 1, weibull_shape = 1)
  long_interval <- 100
  dwell_prevalent_rate <- exact_slip_rate_from_dwell(dwell_standard_model, screen_interval = long_interval, weibull_shape = 1)
  dwell_no_rate <- dwell_slip_rate %>%
    mutate(slip = 1.0, slip_clinical = 1.0)
  
  # Generate performance scenarios
  my_list <- vector("list", 4 * 2 + 2)
  k <- 1
  for (dw_scenario in 1:4) {
    local_performance <- run_intercept_model(
      incidence_sens_source,
      dwell_slip_rate %>% filter(scenario == dw_scenario),
      active_slip_clinical = active_slip
    )
    local_excluded <- run_excluded_model(incidence_excluded_source)
    local_performance <- bind_rows(local_performance, local_excluded) %>%
      mutate(screen_interval = 1, dw_scenario = dw_scenario, scan = "incident")
    
    my_list[[k]] <- local_performance
    k <- k + 1
    
    # Prevalent performance
    prevalent_performance <- run_intercept_model(
      incidence_sens_source,
      dwell_prevalent_rate %>% filter(scenario == dw_scenario),
      active_slip_clinical = active_slip
    )
    prevalent_performance <- prevalent_performance %>%
      filter(found_clinical == 1) %>%
      mutate(caught = caught * long_interval) %>%
      bind_rows(local_performance %>% filter(found_clinical == 2)) %>%
      mutate(screen_interval = 1, dw_scenario = dw_scenario, scan = "prevalent")
    my_list[[k]] <- prevalent_performance
    k <- k + 1
  }
  
  optimal_performance <- run_intercept_model(
    incidence_sens_source,
    dwell_slip_rate %>% filter(scenario == 0),
    active_slip_clinical = active_slip, xStage = xStage
  )
  optimal_excluded <- run_excluded_model(incidence_excluded_source)
  optimal_performance <- bind_rows(optimal_performance, optimal_excluded) %>%
    mutate(opt = "0", dw_scenario = 0, scan = "incident")
  
  no_screening_performance <- run_intercept_model(
    incidence_sens_source,
    dwell_no_rate %>% filter(scenario == 0),
    active_slip_clinical = active_slip, xStage = xStage
  )
  no_screening_performance <- bind_rows(no_screening_performance, optimal_excluded) %>%
    mutate(opt = "NO", dw_scenario = 0, scan = "no")
  
  all_options_df <- bind_rows(my_list, .id = "opt") %>%
    bind_rows(optimal_performance) %>%
    bind_rows(no_screening_performance)
  
  # Generate probability matrices
  Values <- list()
  for (opt_val in unique(all_options_df$opt)) {
    AllData_opt <- subset(all_options_df, opt == opt_val)
    Values[[opt_val]] <- list()
    for (cancer_name in unique(AllData_opt$Cancer)) {
      AllData_cancertype <- subset(AllData_opt, Cancer == cancer_name & found_clinical == 1)
      M <- matrix(0, nrow = 4, ncol = 4)
      if (AllData_cancertype$caught[3] == 0) M[2, 2] <- 1
      if (AllData_cancertype$caught[6] == 0) M[3, 3] <- 1
      if (AllData_cancertype$caught[10] == 0) M[4, 4] <- 1
      M[1, 1] <- 1
      if (AllData_cancertype$caught[3] != 0) {
        M[2, 2] <- AllData_cancertype$caught[3] / sum(AllData_cancertype$caught[2] + AllData_cancertype$caught[3])
      }
      if (AllData_cancertype$caught[6] != 0) {
        M[3, 3] <- AllData_cancertype$caught[6] / sum(AllData_cancertype$caught[4:6])
      }
      if (AllData_cancertype$caught[10] != 0) {
        M[4, 4] <- AllData_cancertype$caught[10] / sum(AllData_cancertype$caught[7:10])
      }
      # Other calculations for matrix entries
      M[1, 2] <- 1 - M[2, 2]
      M[1, 3] <- M[1, 2]
      M[1, 4] <- M[1, 2]
      
      M[2, 3] <- 1 - M[3, 3] - M[1, 3]
      if (M[2, 3] < 0 & !is.na(M[2, 3])) {
        M[3, 3] <- M[3, 3] + M[2, 3]
      }
      M[2, 3] <- 1 - M[3, 3] - M[1, 3]
      M[2, 4] <- M[2, 3]
      
      M[3, 4] <- 1 - M[1, 4] - M[2, 4] - M[4, 4]
      if (M[3, 4] < 0 & !is.na(M[3, 4])) {
        M[4, 4] <- M[4, 4] + M[3, 4]
      }
      M[3, 4] <- 1 - M[1, 4] - M[2, 4] - M[4, 4]
      
      # Store the matrix for this specific Cancer under this specific Opt
      Values[[opt_val]][[cancer_name]] <- M
    }
  }
  
  return(list(sensitivity_long = sensitivity_long, Values = Values, all_options_df = all_options_df))
}

# sensitivity_original <- readr::read_csv("Data/overall/sensitivity_long_overall.csv")
# seer_draw_original <- readr::read_csv("Data/overall/seer_draw.csv")
# dwell_time_original <- readr::read_csv("Data/raw/dwell_model_all_cancers.csv")
# 
# sensitivity_overall_average <- process_sensitivity(cancer_types,
#                                                    sensitivity_original,
#                                                    seer_draw = seer_draw_original,
#                                                    xStage = c("I", "II", "III", "IV"),
#                                                    active_slip = FALSE,
#                                                    dwell_model_df = dwell_time_original)
