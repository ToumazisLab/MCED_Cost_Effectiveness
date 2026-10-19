# setwd("MCED_Cost_Effectiveness")
set.seed(42)
n_psa <- 10000
spec_val <- 0.99 # Change this to adjust MCED specificity

base_dir <- "Sensitivity_analysis"
date_str <- format(Sys.Date(), "%Y-%m-%d")
prefix   <- file.path(base_dir, paste0("Multivar_", date_str))

# Function to find the next available directory
get_next_dir <- function(prefix) {
  i <- 1
  while (dir.exists(paste0(prefix, "_", i))) {
    i <- i + 1
  }
  paste0(prefix, "_", i)
}

# Final path and creation
output_dir <- get_next_dir(prefix)
dir.create(output_dir, recursive = TRUE)

input_params <- list(
  mced_start_cycle = 5,
  mced_end_cycle = 40,
  mced_compliance_rate = 0.75,
  overdiagnosis_factor = 0.05,
  mced_test_cost = 150,
  dwell_scenario = "Scenario_2",
  ss_scenario = 3,
  ts_matrix = "time_shift_matrix_2",
  reduced_comply = 1,
  fp_track_num = 3,
  invasive_multi = 0.1,
  mced_fp_compliance = 1, 
  fp_disutility_val = 0.0125,
  tn_disutility_val = 0,
  d_c = 0.03,
  d_e = 0.03,
  active_slip = FALSE,
  two_step = FALSE,
  overall_sens = 0,
  overall_specificity = 0,
  mced_secondary_test_cost = 0
)

source("Functions/Shiny/Probabilistic.R", echo=TRUE)
