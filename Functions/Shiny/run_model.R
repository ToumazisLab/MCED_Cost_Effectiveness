run_usual_model <- function(cancer_list, data_list, def_param_list, input_params, sensitivity_long, 
                            specificity_data, Values, cancer_types_overall) {
  
  
  #### Set up ####
  screening_cancers <- intersect(cancer_list, def_param_list$screening_cancers)
  
  # sensitivity_results <- process_sensitivity(
  #   cancer_list      = sens_cancers,
  #   sensitivity_long = sensitivity_long,
  #   seer_draw        = data_list$seer_draw,
  #   xStage           = def_param_list$xStage,
  #   active_slip      = FALSE,
  #   dwell_model_df   = data_list$dwell_model_all_cancers
  # )
  # 
  # sensitivity_processed <- sensitivity_results$sensitivity_long %>% 
  #   rename(Sensitivity = sens)
  # 
  # Values <- sensitivity_results$Values
  
  sensitivity_processed <- sensitivity_long
  
  tp_matrices <- generate_tp_matrices(cancer_list, 
                                      data_list$incidence_final, 
                                      data_list$all_cancer_mort, 
                                      def_param_list$female_prop, 
                                      def_param_list$male_prop, 
                                      def_param_list$female_cancers)
  
  final_cumulative_matrix <- run_cancer_model_pipeline(cancer_list, 
                                                       data_list$incidence_final, 
                                                       data_list$all_cancer_mort, 
                                                       def_param_list$female_prop, 
                                                       def_param_list$male_prop, 
                                                       def_param_list$female_cancers, 
                                                       def_param_list$v_names_states, 
                                                       def_param_list$num_people)
  
  result <- calculate_incidents(final_cumulative_matrix)
  incident_matrix <- result$incident_matrix
  incident_matrix_trimmed <- result$incident_matrix_trimmed
  
  final_result_list <- initialize_cancer_metrics(cancer_list)
  
  #### USUAL ARM #### 
  final_result_list <- usual_calculate_tn_fp(final_result_list, cancer_list, incident_matrix, 
                                             data_list$screening_fp_diag_costs, def_param_list$female_cancers, def_param_list$male_cancers, 
                                             def_param_list$female_prop)
  
  final_result_list <- apply_usual_costs(final_result_list, 
                                         data_list$screening_fp_diag_costs, 
                                         discount_factor = input_params$d_c)
  
  final_result_list <- apply_fp_disutility(final_result_list, 
                                           data_list$expanded_utility, 
                                           fp_count = 1, 
                                           fp_disutility_val = input_params$fp_disutility_val, 
                                           tn_disutility_val = input_params$tn_disutility_val)
  
  tp_incidentals <- count_tp_incidentals(incident_matrix_trimmed, 
                                         data_list$screening_fp_diag_costs, 
                                         screening_cancers, 
                                         final_result_list, 
                                         arm = "usual")
  
  final_result_list <- tp_incidentals$final_result_list
  
  final_result_list <- assign_costs(final_result_list, 
                                    data_list$screening_fp_diag_costs, 
                                    screening_cancers, 
                                    discount_factor = input_params$d_c)
  
  u_post_diagnosis_list <- decision_tree(incident_matrix_trimmed, 
                                         data_list$transposed_survival, 
                                         data_list$cost_per_year, 
                                         discount_rate = input_params$d_c, 
                                         cancer_list)
  
  u_post_diagnosis <- as.data.frame(u_post_diagnosis_list$total_costs_incidence)
  u_survival_years <- as.data.frame(u_post_diagnosis_list$survival_df)
  
  u_deaths <- as.matrix(calculate_deaths_by_cancer_stage(incident_matrix_trimmed, 
                                                         data_list$transposed_survival))
  
  u_utility_cancers <- calculate_total_qalys(incident_matrix_trimmed, 
                                             data_list$transposed_survival, 
                                             data_list$utility_long, 
                                             data_list$expanded_utility)
  
  u_utility_cancers_wide <- pivot_wider(u_utility_cancers, names_from = Cancer_Stage, values_from = QALYs, values_fill = list(QALYs = 0))
  u_utility_cancers_wide <- as.data.frame(u_utility_cancers_wide[, -(c(1))])
  

  final_result_list <- update_post_diagnosis_survival_deaths_utility(final_result_list,
                                                                     u_post_diagnosis,
                                                                     u_survival_years,
                                                                     u_deaths,
                                                                     u_utility_cancers_wide,
                                                                     arm = "usual",
                                                                     group = "TP")
  
  
  uc_cancer_deaths <- summarize_cancer_deaths(incident_matrix_trimmed, data_list$survival_expanded)
  #### MCED Arm ####
  results <- calculate_mced_tp_fn_compliance(incident_matrix_trimmed, 
                                  sensitivity_processed, 
                                  final_result_list, 
                                  screening_data = data_list$screening_fp_diag_costs,
                                  cycles = input_params$mced_start_cycle:input_params$mced_end_cycle, 
                                  mced_compliance_rate = input_params$mced_compliance_rate,
                                  two_step = input_params$two_step,
                                  overall_sens = input_params$overall_sens)
  
  final_result_list <- results$final_result_list
  m_count_FN <- as.matrix(results$m_count_FN)
  m_count_TP <- as.matrix(results$m_count_TP)
  m_count_mced_NA <- as.matrix(results$m_count_mced_NA)
  
  overdiagnosis_results <- calculate_overdiagnosis(m_count_TP, 
                                                   incident_matrix, 
                                                   sensitivity_processed, 
                                                   overdiagnosis_factor = input_params$overdiagnosis_factor,
                                                   mced_compliance_rate = input_params$mced_compliance_rate,
                                                   cancer_list,
                                                   all_cancer_list = cancer_types_overall,
                                                   mced_start_cycle = input_params$mced_start_cycle, 
                                                   mced_end_cycle = input_params$mced_end_cycle, 
                                                   final_result_list,
                                                   df = data_list$incidence_final,
                                                   two_step = input_params$two_step,
                                                   overall_sens = input_params$overall_sens)
  
  final_result_list <- overdiagnosis_results$final_result_list
  m_count_TP_overdiagnosis <- overdiagnosis_results$m_count_TP
  m_count_overdiagnosis <- overdiagnosis_results$overdiagnosed_matrix
  
  ts_matrix <- generate_ts_matrix(cancer_list, dwell_data = data_list$dwell_time)
  ## DEPENDENCY: input_params: ss_scenario, ts_matrix
  # III. Shift forward and backward
  results_tp <- shift_backward(cancer_list, 
                               m_count_TP_overdiagnosis, 
                               ss_matrix_list = Values,
                               opt_val = input_params$ss_scenario, 
                               ts_matrix_list = ts_matrix, 
                               input_params$ts_matrix, 
                               input_params$mced_start_cycle, 
                               input_params$mced_end_cycle, 
                               final_result_list)
  
  final_result_list <- results_tp$final_result_list
  m_count_shifted_TP <- results_tp$updated_matrix
  m_count_TP_shifted_only <- results_tp$shifted_cases_matrix
  sankey_df <- results_tp$sankey_df
  tp_transitions      <- results_tp$tp_transitions 
  
  sankey_df_all  <- expand_sankey_with_FN_NA(
    sankey_df, m_count_FN, m_count_mced_NA, cancer_list)
  
  tp_mced <- mced_cancer_categories(m_count_shifted_TP, 
                                    m_count_FN,
                                    m_count_mced_NA, 
                                    data_list$screening_fp_diag_costs,
                                    screening_cancers, 
                                    final_result_list, 
                                    arm = "MCED",
                                    reduced_comply = input_params$reduced_comply)
  
  final_result_list <- tp_mced$final_result_list
  
  final_screen_df <- tp_mced$final_screen_df
  final_screen_df_stage <- tp_mced$final_screen_df_stage
  
  final_result_list <- mced_assign_cost(final_result_list, 
                                        data_list$screening_fp_diag_costs, 
                                        screening_cancers,
                                        arm = "MCED", 
                                        input_params$d_c)
  
  #### X. Decision Tree ####
  m_post_diagnosis_list_tp <- decision_tree(m_count_shifted_TP, 
                                            data_list$transposed_survival, 
                                            data_list$cost_per_year, 
                                            discount_rate = input_params$d_c, 
                                            cancer_list = cancer_list)
  
  m_post_diagnosis_tp <- as.data.frame(m_post_diagnosis_list_tp[["total_costs_incidence"]])
  m_survival_years_tp <- as.data.frame(m_post_diagnosis_list_tp[["survival_df"]])
  
  m_post_diagnosis_list_fn <- decision_tree(m_count_FN,
                                            data_list$transposed_survival, 
                                            data_list$cost_per_year, 
                                            discount_rate = input_params$d_c, 
                                            cancer_list = cancer_list)
  
  m_post_diagnosis_fn <- as.data.frame(m_post_diagnosis_list_fn[["total_costs_incidence"]])
  m_survival_years_fn <- as.data.frame(m_post_diagnosis_list_fn[["survival_df"]])
  
  m_post_diagnosis_list_na <- decision_tree(m_count_mced_NA, 
                                            data_list$transposed_survival, 
                                            data_list$cost_per_year, 
                                            discount_rate = input_params$d_c, 
                                            cancer_list = cancer_list)
  
  m_post_diagnosis_na <- as.data.frame(m_post_diagnosis_list_na[["total_costs_incidence"]])
  m_survival_years_na <- as.data.frame(m_post_diagnosis_list_na[["survival_df"]])
  
  # XI. Deaths
  # m_deaths_tp <- as.matrix(calculate_deaths_by_cancer_stage(m_count_shifted_TP, 
  #                                                           data_list$transposed_survival))
  
  tp_pd <- calculate_tp_post_diagnosis_floored(
    tp_transitions      = results_tp$tp_transitions,
    transposed_survival = data_list$transposed_survival,
    cost_per_year       = data_list$cost_per_year,
    utility_long        = data_list$utility_long,
    expanded_utility    = data_list$expanded_utility,
    discount_rate       = input_params$d_c,
    n_cycles            = nrow(m_count_shifted_TP),
    col_names           = colnames(m_count_shifted_TP),
    survival_method     = input_params$survival_method    # Approach 1; "add_shift" = Approach 2
  )
  m_post_diagnosis_tp       <- tp_pd$total_costs_incidence
  m_survival_years_tp       <- tp_pd$survival_df
  m_deaths_tp               <- tp_pd$deaths
  m_utility_cancers_wide_tp <- tp_pd$qaly_wide
  transition_table          <- tp_pd$transition_table
  
  # m_utility_cancers_tp <- calculate_total_qalys(m_count_shifted_TP, 
  #                                               data_list$transposed_survival,
  #                                               data_list$utility_long, 
  #                                               data_list$expanded_utility)
  # 
  # m_utility_cancers_wide_tp <- pivot_wider(m_utility_cancers_tp, names_from = Cancer_Stage, values_from = QALYs, values_fill = list(QALYs = 0))
  # m_utility_cancers_wide_tp <- as.data.frame(m_utility_cancers_wide_tp[, -(c(1))])
  
  m_deaths_fn <- as.matrix(calculate_deaths_by_cancer_stage(m_count_FN, 
                                                            data_list$transposed_survival))
  
  m_utility_cancers_fn <- calculate_total_qalys(m_count_FN, 
                                                data_list$transposed_survival, 
                                                data_list$utility_long, 
                                                data_list$expanded_utility)
  
  m_utility_cancers_wide_fn <- pivot_wider(m_utility_cancers_fn, names_from = Cancer_Stage, values_from = QALYs, values_fill = list(QALYs = 0))
  m_utility_cancers_wide_fn <- as.data.frame(m_utility_cancers_wide_fn[, -(c(1))])
  
  m_deaths_na <- as.matrix(calculate_deaths_by_cancer_stage(m_count_mced_NA, 
                                                            data_list$transposed_survival))
  
  m_utility_cancers_na <- calculate_total_qalys(m_count_mced_NA, 
                                                data_list$transposed_survival, 
                                                data_list$utility_long, 
                                                data_list$expanded_utility)
  
  m_utility_cancers_wide_na <- pivot_wider(m_utility_cancers_na, names_from = Cancer_Stage, values_from = QALYs, values_fill = list(QALYs = 0))
  m_utility_cancers_wide_na <- as.data.frame(m_utility_cancers_wide_na[, -(c(1))])
  
  
  final_result_list <- update_post_diagnosis_survival_deaths_utility(final_result_list, 
                                                                     m_post_diagnosis_tp, 
                                                                     m_survival_years_tp, 
                                                                     m_deaths_tp,
                                                                     m_utility_cancers_wide_tp, 
                                                                     arm = "MCED", 
                                                                     group = "TP")
  
  final_result_list <- update_post_diagnosis_survival_deaths_utility(final_result_list, 
                                                                     m_post_diagnosis_fn, 
                                                                     m_survival_years_fn, 
                                                                     m_deaths_fn,
                                                                     m_utility_cancers_wide_fn, 
                                                                     arm = "MCED", 
                                                                     group = "FN")
  
  final_result_list <- update_post_diagnosis_survival_deaths_utility(final_result_list, 
                                                                     m_post_diagnosis_na, 
                                                                     m_survival_years_na, 
                                                                     m_deaths_na,
                                                                     m_utility_cancers_wide_na, 
                                                                     arm = "MCED", 
                                                                     group = "NA")
  
  mced_incident_matrix <- m_count_shifted_TP + m_count_FN + m_count_mced_NA
  
  m_cancer_deaths <- summarize_cancer_deaths(mced_incident_matrix, data_list$survival_expanded)
  m_cancer_deaths_TP <- summarize_cancer_deaths(m_count_shifted_TP, data_list$survival_expanded)
  m_cancer_deaths_FN <- summarize_cancer_deaths(m_count_FN, data_list$survival_expanded)
  m_cancer_deaths_NC <- summarize_cancer_deaths(m_count_mced_NA, data_list$survival_expanded)
  #### MCED: TN/FP ####
  ## DEPENDENCY: input_params: mced_compliance_rate
  # MCED TNs and FPs using probabilities
  H_compliant <- generate_compliant_vector(incident_matrix, 
                                           input_params$mced_compliance_rate, 
                                           input_params$mced_start_cycle, 
                                           input_params$mced_end_cycle)
  H_vector <- H_compliant$H_vector
  compliant_vector <- H_compliant$compliant_vector
  
  ## DEPENDENCY: if specificity != data_list$specificity
  ## input_params: fp_track_num
  updated_results <- run_all_cycles(
    compliant_vector = compliant_vector,
    H_vector = H_vector,
    female_prop = def_param_list$female_prop,
    fp_track_num = input_params$fp_track_num,
    specificity_df = specificity_data,
    cancer_list = cancer_list,
    female_cancers = def_param_list$female_cancers,
    male_cancers = def_param_list$male_cancers,
    final_result_list = final_result_list,
    two_step = input_params$two_step,
    overall_specificity = input_params$overall_specificity,
    incident_matrix_trimmed = mced_incident_matrix
  )
  
  final_result_list <- updated_results$final_result_list
  
  ## DEPENDENCY: input_params: mced_fp_compliance
  final_result_list <- mced_fp_usual_compliance(final_result_list, 
                                                cancer_list, 
                                                mced_fp_compliance = input_params$mced_fp_compliance)
  
  ## DEPENDENCY: input_params: tn_skip_years
  # IV. Calculate MCED Usual Care
  final_result_list <- apply_mced_usual_care_and_costs(final_result_list, 
                                                       data_list$screening_fp_diag_costs,
                                                       discount_factor = input_params$d_c, 
                                                       screening_cancers,
                                                       tn_skip_factor = input_params$reduced_comply,
                                                       invasive_multiplier = input_params$invasive_multi)
  
  final_result_list <- summarize_overall_results(final_result_list)
  
  ## MCED Non-compliant (1-mced_compliance_rate * usual care values for FP, TN, Invasive, FP_imaging, FP_diagnostic, TN_imaging, Invasive costs)
  final_result_list <- apply_mced_nc_calculations(final_result_list, 
                                                  data_list$screening_fp_diag_costs, 
                                                  input_params$d_c, 
                                                  input_params$mced_compliance_rate)
  
  ## MCED Test Cost
  final_result_list <- apply_mced_test_cost(final_result_list,
                                            input_params$mced_test_cost,
                                            input_params$mced_secondary_test_cost,
                                            input_params$d_c)
  
  ### MCED TN/FP Utility
  final_result_list <- calculate_mced_utility(final_result_list, 
                                              data_list$expanded_utility,
                                              fp_disutility = input_params$fp_disutility_val,
                                              def_param_list$female_cancers, 
                                              def_param_list$male_cancers, 
                                              def_param_list$female_prop)
  
  
  
  final_result_list <- sum_mced_usual_costs(final_result_list)
  final_result_list <- calculate_icer(final_result_list)
  
  overall_summary <- generate_full_long_df(final_result_list, cancer_list,
                                           m_cancer_deaths, uc_cancer_deaths, final_screen_df)
  
  return(list(final_result_list = final_result_list, mced_incident_matrix = mced_incident_matrix, 
              incident_matrix_trimmed = incident_matrix_trimmed, sankey_df_all = sankey_df_all,
              overall_summary = overall_summary, final_screen_df_stage = final_screen_df_stage,
              m_cancer_deaths_TP = m_cancer_deaths_TP, m_cancer_deaths_FN = m_cancer_deaths_FN, 
              m_cancer_deaths_NC = m_cancer_deaths_NC, m_cancer_deaths = m_cancer_deaths, 
              uc_cancer_deaths = uc_cancer_deaths, final_screen_df = final_screen_df, transition_table = transition_table))
  
}



