source("Global.R")
# source("Functions/Preparation_Scripts/Source_file.R", echo=TRUE)
source("Functions/Preparation_Scripts/Def_params.R", echo=TRUE)
source("Functions/Preparation_Scripts/Source_file.R", echo=TRUE)


# Only Change this List when rerunning the model
input_params <- list(
  mced_start_cycle = 6,
  mced_end_cycle = 41,
  mced_compliance_rate = 0.75,
  overdiagnosis_factor = 0.05,
  mced_test_cost = 500,
  dwell_scenario = "Scenario_2", # change to 1 and 3
  ss_scenario = 3, # change to 1 and 5
  ts_matrix = "time_shift_matrix_2", # change to 1 and 3
  reduced_comply = 1,
  fp_track_num = 3,
  invasive_multi = 0.1, # Must be between 0-1
  mced_fp_compliance = 1, 
  fp_disutility_val = 0.0125,
  tn_disutility_val = 0,
  d_c = 0.03,
  d_e = 0.03,
  active_slip = TRUE,
  two_step = FALSE,
  overall_sens = 0,
  overall_specificity = 0,
  mced_secondary_test_cost = 0
)

data_loaded <- load_data(input_params, def_param_list)
# For the sensitivity analysis we will work with all 16 cancers, and vary input parameters one at a time (+/- 25% of their base-value 
# 
# Fixed screening period to 45-80
# Compliance rate to mced (input_params$mced_compliance_rate - bounded)
# Overdiagnosis rate (input_params$overdiagnosis_factor bounded)
# Mced test cost (input_params$mced_test_cost - manual: 100, 1000)
# Dwell scenarios 1 and 3 Keep this the same
# Cancer and stage specific sensitivities (all cancers and stages +/- 25%; not all combinations) (increase/decrease all sensitivities ensuring not to go over 1)
# Reduce compliance to UC after an MCED result (input_params$mced_fp_compliance)
# Change number of FP to act to 1-5 (input_params$fp_track_num)
# Invasive_multi (input_params$invasive_multi)
# All the disutilities/utilities (input_params$fp_disutility_val (not to go below zero) and the utility data)
# Discount factors (input_params$d_e, input_params$d_c)
# Cost screening (data_loaded$screening_fp_diag_costs$Screening_cost FOR SCREENING CANCERS ONLY)
# Cost of treatment (all combined - by cancer, stage and year) (cost_per_year data where columns are Cancer_Stage, Year, Cost - long format)
# Cost of invasive procedures (data_loaded$screening_fp_diag_costs$Invasive_procedure to not go over 1)
# Change in survival for all cancers (transposed_survival, columns are named with convention Stage_Cancer - wide format)
# Change in UC compliance (data_loaded$screening_fp_diag_costs$Compliance_frequency - ONLY FOR SCREENING CANCERS)
# Cost of FP and incidental (data_loaded$screening_fp_diag_costs$FP_cost & data_loaded$screening_fp_diag_costs$incidental_cost to be adjusted up/down TOGETHER)


sens_cancers <- cancer_types_overall <- cancer_types
# data_loaded <- load_data(input_params, def_param_list)
df_sensitivity <- read_csv("Data/overall/sensitivity_long_overall.csv")
original_specificity <- data.frame(
  Cancer = cancer_types,
  # Default starting value for specificity
  Specificity = rep(0.98, 16)
)

sensitivity_results <- process_sensitivity(
  cancer_list      = sens_cancers,
  sensitivity_long = df_sensitivity,
  seer_draw        = data_loaded$seer_draw,
  xStage           = def_param_list$xStage,
  active_slip      = TRUE,
  dwell_model_df   = data_loaded$dwell_model_all_cancers
)

sensitivity_processed <- sensitivity_results$sensitivity_long
sensitivity_matrix <- sensitivity_results$Values



library(dplyr)
library(tidyr)
library(purrr)
library(glue)

# --------------------------------------------------
# helper: two‑row ±25 % design for scalars
# --------------------------------------------------
two_rows_bounded <- function(name, base, pct = 0.25) {
  tibble(
    param_name  = name,
    change_type = c("minus","plus"),
    # e.g. base * (1 - 0.25) and base * (1 + 0.25),
    # but clamp to [0,1]
    new_value   = c(base*(1-pct), base*(1+pct)) |>
      pmin(1) |> pmax(0)
  )
}

two_rows_unbounded <- function(name, base, pct = 0.25) {
  tibble(
    param_name  = name,
    change_type = c("minus","plus"),
    # no upper clamp; possibly clamp at 0 if negative doesn't make sense
    new_value   = c(base*(1-pct), base*(1+pct))
  )
}

sens_design <- tibble(
  param_name  = c("all_sens", "all_sens"),
  change_type = c("minus",            "plus"),
  new_value   = c(0.75,               1.25)
)

spec_design <- tibble(
  param_name  = c("all_spec", "all_spec"),
  change_type = c("minus", "plus"),
  new_value   = c(0.75, 0.999)
)

# 1A. scalar params you want to vary
mced_fp_design <- tibble(
  param_name  = c("mced_fp_compliance", "mced_fp_compliance"),
  change_type = c("minus", "plus"),
  new_value   = c((input_params$mced_fp_compliance * 0.75), 1)
)

scalar_design <- bind_rows(
  two_rows_bounded("mced_compliance_rate", input_params$mced_compliance_rate),
  # two_rows_bounded("overdiagnosis_factor", input_params$overdiagnosis_factor),
  
  tibble(
    param_name = "overdiagnosis_factor", 
    change_type = c("minus", "plus"),
    new_value = c(0.01, 0.1)
  ),
  
  # Explicit values for mced_test_cost instead of +/- 25%
  tibble(
    param_name  = "mced_test_cost",
    change_type = c("minus", "plus"),
    new_value   = c(50, 1000)
  ),
  
  tibble(
    param_name  = "fp_track_num",
    change_type = c("minus", "plus"),
    new_value   = c(1, 5)
  ),
  
  tibble(
    param_name  = "invasive_multi_mult",
    change_type = c("minus", "plus"),
    new_value   = c(0.075, 0.0125)
  ),
  
  two_rows_unbounded("fp_disutility_val",   input_params$fp_disutility_val),
  two_rows_bounded("d_c",                 input_params$d_c),
  two_rows_bounded("d_e",                 input_params$d_e)
)


# 1B. dwell‑scenario values you want to test
dwell_design <- tibble(
  param_name  = "dwell_scenario",
  change_type = c("plus", "minus"),
  new_value   = c(1, 2)  # Use scalar identifiers
)
# 1C. screening‑cost ±25% (applied as a % change)
screen_design <- tribble(
  ~param_name,         ~change_type, ~new_value,
  "screen_cost_pct",   "minus",      -0.25,
  "screen_cost_pct",   "plus",        0.25
)

# USPSTF‑screened cancers
screen_cancers <- c("lung and bronchus","colon and rectum",
                    "breast","prostate","cervix")

# ---- tweak utility (Years 1‑4 only) ----------
util_mult <- function(df, mult) {
  df |> mutate(
    utility = ifelse(
      Year == "Year 5+",
      utility,
      pmin(utility * mult, 1)
    )
  )
}

cpy_extreme_by_cancer <- function(df, type = c("plus", "minus")) {
  type <- match.arg(type)
  
  df %>%
    separate(Cancer_Stage, into = c("Stage", "Cancer"), sep = "_", remove = FALSE) %>%
    group_by(Cancer, Year) %>%
    mutate(Cost = case_when(
      type == "plus" ~ max(Cost, na.rm = TRUE),
      type == "minus" ~ min(Cost, na.rm = TRUE)
    )) %>%
    ungroup() %>%
    select(-Stage, -Cancer)
}

# ---- tweak survival --------------------------
surv_mult <- function(df, mult) df * mult      # matrix or wide df is fine

# ---- tweak cost‑per‑year (Years 1‑4 only) ----
cpy_mult <- function(df, mult) {
  df |> mutate(
    Cost = ifelse(Year == "Year 5+", Cost, Cost * mult)
  )
}

cpy_mult_one_stage <- function(df, target_stage, mult = 1.25) {
  df %>%
    tidyr::separate(Cancer_Stage, into = c("Stage", "Cancer"), sep = "_", remove = FALSE) %>%
    dplyr::mutate(
      Cost = ifelse(
        Stage == target_stage & Year != "Year 5+",
        Cost * mult,
        Cost
      )
    ) %>%
    dplyr::select(-Stage, -Cancer)
}

cpy_extreme_design <- tibble(
  param_name = "cpy_extreme",
  change_type = c("plus", "minus"),
  new_value = c(1, 2)  
)

cpy_by_stage_design <- expand.grid(
  stage = c("I", "II", "III", "IV"),
  change_type = c("minus", "plus"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    param_name  = paste0("stage_cost_", stage),
    new_value   = ifelse(change_type == "minus", 0.75, 1.25)
  ) %>%
  select(param_name, change_type, new_value)

inv_cost_by_cancer_design <- data_loaded$screening_fp_diag_costs %>%
  select(Cancer, Invasive_procedure) %>%
  distinct() %>%
  filter(!is.na(Invasive_procedure)) %>%
  mutate(
    param_name_minus = paste0("invasive_proc_cost_", Cancer, "_minus"),
    param_name_plus  = paste0("invasive_proc_cost_", Cancer, "_plus"),
    new_value_minus  = Invasive_procedure * 0.75,
    new_value_plus   = Invasive_procedure * 1.25
  ) %>%
  select(Cancer, starts_with("param_name"), starts_with("new_value")) %>%
  pivot_longer(cols = starts_with("param_name"), names_to = "change_type", values_to = "param_name") %>%
  mutate(
    change_type = ifelse(str_detect(change_type, "minus"), "minus", "plus"),
    new_value = ifelse(change_type == "minus", new_value_minus, new_value_plus)
  ) %>%
  filter(
    new_value != 0
  ) %>% 
  select(param_name, change_type, new_value)

compliance_freq_mult <- function(df, pct, screen_cancers) {
  m <- 1 + pct
  df |> mutate(
    Compliance_frequency = if_else(Cancer %in% screen_cancers,
                                   Compliance_frequency * m,
                                   Compliance_frequency)
  )
}


adjust_all_cause_mort <- function(all_cause_raw, cancer_mort, cancer_list,
                                  mult_all_cause = 1, mult_total_cancer = 1) {
  all_cause_adj <- all_cause_raw
  all_cause_adj$All_cause_mortality <- all_cause_adj$All_cause_mortality * mult_all_cause
  
  filtered_cancer <- subset(cancer_mort, Cancer %in% cancer_list)
  total_cancer_mort <- aggregate(probability ~ Age, data = filtered_cancer, sum)
  total_cancer_mort$probability <- total_cancer_mort$probability * mult_total_cancer
  colnames(total_cancer_mort)[2] <- "Total_cancer_probability"
  
  merged <- merge(all_cause_adj, total_cancer_mort, by = "Age", all.x = TRUE)
  merged$Total_cancer_probability[is.na(merged$Total_cancer_probability)] <- 0
  
  merged$Non_cancer_all_cause_probability <-
    merged$All_cause_mortality - merged$Total_cancer_probability
  
  return(merged)
}

util_design <- two_rows_unbounded("util_mult", 1, pct = 0.25)

cost_design <- tribble(
  ~param_name,        ~change_type, ~new_value,
  "fp_incidental_cost_pct", "minus", -0.25,
  "fp_incidental_cost_pct", "plus",   0.25,
  "invasive_cost_pct","minus",       -0.25,
  "invasive_cost_pct","plus",         0.25
)

# 1D. survival & cost‑per‑year multipliers (±25 %)
surv_cpy_design <- bind_rows(
  two_rows_unbounded("survival_mult",   1, pct = 0.25),
  two_rows_unbounded("cost_per_year_mult", 1, pct = 0.25)
)

# ---- tweak individual screening‑cost cells ---
cost_component_pct <- function(df, pct, col) {
  m <- 1 + pct
  df |> mutate("{col}" := .data[[col]] * m)
}
# ---- tweak screening cost table ----------------
screen_cost_pct <- function(df, pct) {
  m <- 1 + pct
  df |> mutate(
    Screening_cost = if_else(Cancer %in% screen_cancers,
                             Screening_cost*m, Screening_cost)
  )
}

mortality_design <- bind_rows(
  two_rows_unbounded("all_cause_mortality_mult", 1, pct = 0.25),
  two_rows_unbounded("total_cancer_mortality_mult", 1, pct = 0.25)
)

# bind everything together
dsa_design <- bind_rows(
  scalar_design, inv_cost_by_cancer_design, dwell_design, screen_design, mced_fp_design, cpy_by_stage_design,
  cost_design, surv_cpy_design, mortality_design, sens_design, spec_design, util_design) |>
  mutate(
    run_id = row_number(),
    label = paste0(param_name, "_", change_type)
  )

# --------------------------------------------------
tmp_ip <- input_params
tmp_row <- dsa_design %>% slice(1)

run_one <- function(design_row,
                    base_params,
                    base_data,
                    df_sensitivity,
                    original_specificity,
                    sensitivity_matrix,
                    def_param_list,
                    cancer_types_overall,
                    result_dir) {
  
  cat(glue::glue("\nNow running scenario: {design_row$run_id} ",
                 "({design_row$param_name} – {design_row$change_type})\n"))
  
  ip <- base_params
  dl <- base_data
  get_scalar <- function(x) as.numeric(x[[1]])
  
  multiplier <- as.numeric(design_row$new_value)
  
  # By default, pass the original data frames
  local_sens <- df_sensitivity
  local_spec <- original_specificity
  
  # If the row is about sensitivity:
  if (design_row$label == "all_sens_plus") {
    
    # Here 'multiplier' is 0.75 or 1.25
    local_sens <- local_sens %>%
      mutate(
        Sensitivity = Sensitivity * multiplier,
        Sensitivity = pmin(1, pmax(0, Sensitivity))
      )
  }
  
  if (design_row$label == "all_sens_minus") {
    
    # Here 'multiplier' is 0.75 or 1.25
    local_sens <- local_sens %>%
      mutate(
        Sensitivity = Sensitivity * multiplier,
        Sensitivity = pmin(1, pmax(0, Sensitivity))
      )
  }
  
  if (grepl("^stage_cost_", design_row$param_name)) {
    # Extract stage from param_name (e.g., "stage_cost_II" → "II")
    target_stage <- sub("stage_cost_", "", design_row$param_name)
    
    dl$cost_per_year <- cpy_mult_one_stage(
      df = dl$cost_per_year,
      target_stage = target_stage,
      mult = as.numeric(design_row$new_value)
    )
  }
  
  if (design_row$label == "all_spec_minus") {
    # Suppose 'multiplier' = 0.75
    local_spec <- local_spec %>%
      mutate(
        Specificity = Specificity * multiplier,
        Specificity = pmin(1, pmax(0, Specificity))
      )
  }
  
  if (design_row$label == "all_spec_plus") {
    # Suppose 'multiplier' = 0.75
    local_spec <- local_spec %>%
      mutate(
        Specificity = Specificity * multiplier,
        Specificity = pmin(1, pmax(0, Specificity))
      )
  }
  
  # ---- 3A. scalar replacements (unchanged) ----
  if (design_row$param_name %in%
      c("mced_compliance_rate","overdiagnosis_factor","mced_test_cost",
        "fp_track_num","mced_fp_compliance","fp_disutility_val",
        "d_c","d_e")) {
    ip[[design_row$param_name]] <- get_scalar(design_row$new_value)
  }
  
  if (design_row$param_name == "Compliance_frequency") {
    dl$screening_fp_diag_costs <- compliance_freq_mult(
      dl$screening_fp_diag_costs,
      get_scalar(design_row$new_value),
      screen_cancers
    )
  }
  
  
  if (design_row$param_name == "cpy_by_stage") {
    dl$cost_per_year <- cpy_mult_one_stage(
      df = dl
    )
  } 
  
  if (design_row$param_name == "cpy_extreme") {
    type <- if (design_row$new_value == 1) "plus" else "minus"
    dl$cost_per_year <- cpy_extreme_by_cancer(
      df = dl$cost_per_year,
      type = type
    )
  }
  
  # invasive_multi uses a multiplier
  if (design_row$param_name == "invasive_multi_mult") {
    ip$invasive_multi <- get_scalar(design_row$new_value)
  }
  
  # ---- 3B. linked triplet (unchanged) ---------
  if (design_row$param_name == "dwell_scenario") {
    if (design_row$new_value == 1) {
      ip$dwell_scenario <- "Scenario_1"
      ip$ss_scenario    <- 1
    } else if (design_row$new_value == 2) {
      ip$dwell_scenario <- "Scenario_3"
      ip$ss_scenario    <- 5
    } else {
      stop("Unknown dwell_scenario value")
    }
  }
  
  # ---- 3C. cost tables ------------------------
  if (design_row$param_name == "screen_cost_pct") {
    pct <- get_scalar(design_row$new_value)
    dl$screening_fp_diag_costs <-
      screen_cost_pct(dl$screening_fp_diag_costs, pct)
  }
  
  if (design_row$param_name == "fp_incidental_cost_pct") {
    pct <- get_scalar(design_row$new_value)
    dl$screening_fp_diag_costs <- dl$screening_fp_diag_costs |>
      cost_component_pct(pct, "FP_cost") |>
      cost_component_pct(pct, "incidental_cost")
  }
  if (design_row$param_name == "invasive_cost_pct") {
    dl$screening_fp_diag_costs <-
      cost_component_pct(dl$screening_fp_diag_costs,
                         get_scalar(design_row$new_value), "invasive_cost")
  }
  
  # ---- 3D. utilities, survival, yearly costs --
  if (design_row$param_name == "util_mult") {
    dl$utility_long <- util_mult(dl$utility_long, get_scalar(design_row$new_value))
  }
  if (design_row$param_name == "survival_mult") {
    dl$transposed_survival <- surv_mult(dl$transposed_survival,
                                        get_scalar(design_row$new_value))
  }
  if (design_row$param_name == "cost_per_year_mult") {
    dl$cost_per_year <- cpy_mult(dl$cost_per_year,
                                 get_scalar(design_row$new_value))
  }
  
  # If d_e is modified, recalculate expanded_utility
  if (design_row$param_name == "d_e") {
    dl$expanded_utility <- update_utility(
      d_e = get_scalar(design_row$new_value),
      raw_data_folder = "Data/raw",
      female_prop = def_param_list$female_prop
    )
  }
  
  # If mortality is modified, adjust all_cause_mort
  if (design_row$param_name == "all_cause_mortality_mult") {
    dl$all_cancer_mort <- adjust_all_cause_mort(
      all_cause_raw = read_csv("Data/raw/all_cause_mortality_40.csv"),
      cancer_mort   = read_csv("Data/raw/cancer_mort.csv"),
      cancer_list   = cancer_types,
      mult_all_cause = get_scalar(design_row$new_value),
      mult_total_cancer = 1  # leave cancer unchanged
    )
  }
  
  if (design_row$param_name == "total_cancer_mortality_mult") {
    dl$all_cancer_mort <- adjust_all_cause_mort(
      all_cause_raw = read_csv("Data/raw/all_cause_mortality_40.csv"),
      cancer_mort   = read_csv("Data/raw/cancer_mort.csv"),
      cancer_list   = cancer_types,
      mult_all_cause = 1,
      mult_total_cancer = get_scalar(design_row$new_value)
    )
  }
  
  # ---- 3E. run the model ----------------------
  res <- run_usual_model(
    cancer_list          = cancer_types_overall,
    data_list            = dl,
    def_param_list       = def_param_list,
    input_params         = ip,
    original_sensitivity = local_sens,         # <--- use local_sens
    specificity_data     = local_spec,         # <--- use local_spec
    Values               = sensitivity_matrix,
    cancer_types_overall = cancer_types_overall
  )
  
  # ---- 3F. pull ICERs -------------------------
  icer_df <- res$final_result_list$Overall$final_icer_result
  out_tibble <- tibble(
    run_id              = design_row$run_id,
    param_name          = design_row$param_name,
    change_type         = design_row$change_type,
    new_value           = design_row$new_value,
    Incremental_Cost    = icer_df$Value[icer_df$Variable=="Incremental_Cost"],
    Incremental_Utility = icer_df$Value[icer_df$Variable=="Incremental_Utility"],
    ICER                = icer_df$Value[icer_df$Variable=="ICER"]
  )
  
  message(
    "\n---> Results for run ",
    design_row$run_id,
    " (", design_row$param_name, " – ", design_row$change_type, "):\n"
  )
  print(out_tibble)
  
  # ---- 3G. save final_result_list to disk ----
  # e.g. "final_results_1.rds", "final_results_2.rds", etc.
  out_rds_file <- file.path(
    result_dir,
    paste0("final_results_", design_row$run_id, ".rds")
  )
  saveRDS(res$final_result_list, out_rds_file)
  
  # Return your summary tibble
  return(out_tibble)
}

run_sensitivity_analysis <- function(dsa_design, 
                                     input_params, 
                                     data_loaded,
                                     df_sensitivity,
                                     original_specificity,
                                     sensitivity_matrix,
                                     def_param_list,
                                     cancer_types_overall,
                                     folder_suffix = "_1") {
  
  # 1. Create top-level folder if it doesn't exist
  dir.create("Sensitivity_analysis", showWarnings = FALSE)
  
  # 2. Create subfolder named with today's date (plus optional suffix)
  #    e.g. "2023-05-12_1"
  date_str   <- as.character(Sys.Date())  # e.g. "2023-05-12"
  result_dir <- file.path("Sensitivity_analysis", paste0(date_str, folder_suffix))
  dir.create(result_dir, showWarnings = FALSE)
  
  # 3. Write the design CSV to that subfolder
  readr::write_csv(dsa_design, file.path(result_dir, "dsa_design.csv"))
  
  # 4. Run each row of the design, returning a combined data frame
  dsa_results <- purrr::map_dfr(
    split(dsa_design, dsa_design$run_id),
    run_one,
    base_params         = input_params,
    base_data           = data_loaded,
    df_sensitivity      = df_sensitivity,
    original_specificity = original_specificity,
    sensitivity_matrix  = sensitivity_matrix,
    def_param_list      = def_param_list,
    cancer_types_overall = cancer_types_overall,
    result_dir          = result_dir  # pass folder path
  )
  
  # 5. Optionally, also save the final aggregated results
  readr::write_csv(dsa_results, file.path(result_dir, "dsa_results.csv"))
  
  return(invisible(dsa_results))
}

# ------------------------------------------------
# Load static data once
# ------------------------------------------------
data_loaded <- load_data(input_params, def_param_list)

# ------------------------------------------------
# Run every row in the design matrix
# ------------------------------------------------
results_07_02_2026 <- run_sensitivity_analysis(
  dsa_design          = dsa_design,
  input_params        = input_params,
  data_loaded         = data_loaded,
  df_sensitivity      = df_sensitivity,
  original_specificity = original_specificity,
  sensitivity_matrix  = sensitivity_matrix,
  def_param_list      = def_param_list,
  cancer_types_overall = cancer_types_overall,
  folder_suffix       = "_1"  # or "_2", etc.
)

# ------------------------------------------------
# Save & (optionally) plot
# ------------------------------------------------

baseline <- run_usual_model(cancer_types_overall,
                            data_loaded,
                            def_param_list,
                            input_params = input_params,
                            df_sensitivity,
                            original_specificity,
                            sensitivity_matrix,
                            cancer_types_overall)

# quick tornado
library(ggplot2)
tornado <- results %>%
  select(param_name, change_type, ICER) %>%
  pivot_wider(
    names_from = change_type, 
    values_from = ICER,
    names_prefix = "ICER_",
    values_fn = mean  # Change this to sum, median, etc. if needed.
  ) %>%
  mutate(delta = ICER_plus - ICER_minus) %>%
  arrange(desc(abs(delta)))

ggplot(tornado,
       aes(x = reorder(param_name, abs(delta)),
           ymin = pmin(ICER_minus, ICER_plus),
           ymax = pmax(ICER_minus, ICER_plus))) +
  geom_linerange(size = 2) +
  coord_flip() +
  labs(x = "Parameter", y = "ICER",
       title = "Univariate Sensitivity Analysis")

tornado_cost <- results %>%
  select(param_name, change_type, Incremental_Cost) %>%
  pivot_wider(
    names_from = change_type, 
    values_from = Incremental_Cost,
    names_prefix = "ICER_",
    values_fn = mean  # Change this to sum, median, etc. if needed.
  ) %>%
  mutate(delta = ICER_plus - ICER_minus) %>%
  arrange(desc(abs(delta)))

ggplot(tornado_cost,
       aes(x = reorder(param_name, abs(delta)),
           ymin = pmin(ICER_minus, ICER_plus),
           ymax = pmax(ICER_minus, ICER_plus))) +
  geom_linerange(size = 2) +
  coord_flip() +
  labs(x = "Parameter", y = "Incremental Cost",
       title = "Univariate Sensitivity Analysis - Cost")

tornado_util <- results %>%
  select(param_name, change_type, Incremental_Utility) %>%
  pivot_wider(
    names_from = change_type, 
    values_from = Incremental_Utility,
    names_prefix = "ICER_",
    values_fn = mean  # Change this to sum, median, etc. if needed.
  ) %>%
  mutate(delta = ICER_plus - ICER_minus) %>%
  arrange(desc(abs(delta)))

ggplot(tornado_util,
       aes(x = reorder(param_name, abs(delta)),
           ymin = pmin(ICER_minus, ICER_plus),
           ymax = pmax(ICER_minus, ICER_plus))) +
  geom_linerange(size = 2) +
  coord_flip() +
  labs(x = "Parameter", y = "Incremental Utility",
       title = "Univariate Sensitivity Analysis - Utility")
