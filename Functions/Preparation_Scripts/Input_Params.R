# Change using names in cancer_types (loaded list)
cancer_list <- c("colon and rectum", "breast: hormone receptor-negative", 
                 "lung and bronchus", "pancreas", "esophagus", "liver and intrahepatic bile duct", 
                 "ovarian")


# Only Change this List when rerunning the model
input_params <- list(
  mced_start_cycle = 6,
  mced_end_cycle = 41,
  mced_compliance_rate = 0.75,
  overdiagnosis_factor = 0.05,
  mced_test_cost = 500,
  dwell_scenario = "Scenario_2",
  ss_scenario = 3,
  ts_matrix = "time_shift_matrix_2",
  reduced_comply = 1,
  fp_track_num = 3,
  invasive_multi = 0.1, # Must be between 0-1
  mced_fp_compliance = 1, 
  fp_disutility_val = 0.0125,
  tn_disutility_val = 0,
  d_c = 0.03,
  d_e = 0.03,
  active_slip = FALSE,
  two_step = FALSE,
  overall_sens = 0,
  overall_specificity = 0,
  mced_secondary_test_cost = 0,
  survival_method = "floor_max"
)
# 
# 
# ### ORIGINAL PARAMETERS
# original_params <- list(
#   len_cancer_list = length(cancer_types_MDA),
#   mced_start_cycle = 1,
#   mced_end_cycle = 30,
#   mced_compliance_rate = 0.9,
#   overdiagnosis_factor = 0.05,
#   mced_test_cost = 150,
#   dwell_scenario = "Scenario_1",
#   ss_scenario = 1,
#   ts_matrix = "time_shift_matrix_1",
#   reduced_comply = 2,
#   fp_track_num = 3,
#   invasive_multi = 0.1,
#   mced_fp_compliance = 1, 
#   fp_disutility_val = 0.0125,
#   tn_disutility_val = 0,
#   sens = "normal",
#   spec = "normal",
#   d_c = 0.03,
#   d_e = 0.03,
#   active_slip = FALSE
# )
