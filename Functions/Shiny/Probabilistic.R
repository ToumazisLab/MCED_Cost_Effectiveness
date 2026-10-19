# setwd("MCED_Cost_Effectiveness")
source("Global.R")
# source("Functions/Preparation_Scripts/Source_file.R", echo=TRUE)
source("Functions/Preparation_Scripts/Def_params.R", echo=TRUE)
source("Functions/Preparation_Scripts/Source_file.R", echo=TRUE)


# Change this List when rerunning the model


data_loaded <- load_data(input_params, def_param_list)
sens_cancers <- cancer_types_overall <- cancer_types
# data_loaded <- load_data(input_params, def_param_list)
df_sensitivity <- read_csv("Data/overall/sensitivity_long_overall.csv")
original_specificity <- data.frame(
  Cancer = cancer_types,
  # Default starting value for specificity
  Specificity = rep(spec_val, 16)
)

sensitivity_results <- process_sensitivity(
  cancer_list      = sens_cancers,
  sensitivity_long = df_sensitivity,
  seer_draw        = data_loaded$seer_draw,
  xStage           = def_param_list$xStage,
  active_slip      = FALSE,
  dwell_model_df   = data_loaded$dwell_model_all_cancers
)

sensitivity_processed <- sensitivity_results$sensitivity_long
sensitivity_matrix <- sensitivity_results$Values

# 2A. Utility multipliers (Years 1‑4 only)
util_mult <- function(df, mult) {
  df |> mutate(
    utility = ifelse(
      Year == "Year 5+",
      utility,
      pmin(utility * mult, 1)
    )
  )
}

# 2B. Cost-per-year extremes by cancer (replaces cost by min or max)
cpy_extreme_by_cancer <- function(df, type = c("plus", "minus")) {
  type <- match.arg(type)
  
  df %>%
    separate(Cancer_Stage, into = c("Stage", "Cancer"), sep = "_", remove = FALSE) %>%
    group_by(Cancer, Year) %>%
    mutate(Cost = case_when(
      type == "plus"  ~ max(Cost, na.rm = TRUE),
      type == "minus" ~ min(Cost, na.rm = TRUE)
    )) %>%
    ungroup() %>%
    select(-Stage, -Cancer)
}

cost_extreme_max <- cpy_extreme_by_cancer(data_loaded$cost_per_year, type = "plus")
cost_extreme_min <- cpy_extreme_by_cancer(data_loaded$cost_per_year, type = "minus")

# 2C. Survival multiplier
surv_mult <- function(df, mult) {
  df * mult  # works for matrix or wide df
}

# 2D. Cost-per-year multiplier (Years 1‑4 only)
cpy_mult <- function(df, mult) {
  df |> mutate(
    Cost = ifelse(Year == "Year 5+", Cost, Cost * mult)
  )
}

# 2E. Cost-per-year multiplier for a single stage
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

# 2F. Adjust compliance frequencies (for screening)
compliance_freq_mult <- function(df, pct, screen_cancers) {
  m <- 1 + pct
  df |> mutate(
    Compliance_frequency = if_else(
      Cancer %in% screen_cancers,
      Compliance_frequency * m,
      Compliance_frequency
    )
  )
}

# 2G. Adjust all-cause mortality to reflect changes in total cancer mortality
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

# 2H. Cost multipliers for screening-cost tables
cost_component_pct <- function(df, pct, col) {
  m <- 1 + pct
  df %>% 
    mutate(
      "{col}" := case_when(
        .data[[col]] == 0 ~ .data[[col]],
        TRUE              ~ .data[[col]] * m
      )
    )
}


screen_cost_pct <- function(df, pct) {
  m <- 1 + pct
  df |> mutate(
    Screening_cost = if_else(Cancer %in% screen_cancers,
                             Screening_cost*m, 
                             Screening_cost)
  )
}

##############################################################################
# 3. DESIGN TIBBLES FOR SENSITIVITY ANALYSIS
##############################################################################

# 3A. Pre-defined sensitivity & specificity changes
sens_design <- tibble(
  param_name  = c("all_sens", "all_sens"),
  change_type = c("minus",    "plus"),
  new_value   = c(0.75,       1.25)
)

spec_design <- tibble(
  param_name  = c("all_spec", "all_spec"),
  change_type = c("minus",    "plus"),
  new_value   = c(0.75,       0.999)
)

# 3B. Example: scalar params
mced_fp_design <- tibble(
  param_name  = c("mced_fp_compliance", "mced_fp_compliance"),
  change_type = c("minus", "plus"),
  new_value   = c((input_params$mced_fp_compliance * 0.75), 1)
)

scalar_design <- bind_rows(
  two_rows_bounded("mced_compliance_rate", input_params$mced_compliance_rate),
  two_rows_bounded("overdiagnosis_factor", input_params$overdiagnosis_factor),
  
  # Explicit values for mced_test_cost instead of +/- 25%
  tibble(
    param_name  = "mced_test_cost",
    change_type = c("minus", "plus"),
    new_value   = c(100, 1000)
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
  
  two_rows_unbounded("fp_disutility_val", input_params$fp_disutility_val),
  two_rows_bounded("d_c",                 input_params$d_c),
  two_rows_bounded("d_e",                 input_params$d_e)
)

# 3C. Dwell-scenario variation
dwell_design <- tibble(
  param_name  = "dwell_scenario",
  change_type = c("plus", "minus", "base"),
  new_value   = c(1, 2, 3)  # Use scalar identifiers
)

# 3D. Screening-cost ±25 % (applied as a % change)
screen_design <- tribble(
  ~param_name,       ~change_type, ~new_value,
  "screen_cost_pct", "minus",      -0.25,
  "screen_cost_pct", "plus",        0.25
)

# 3E. Cost per year (extreme) and by stage
cpy_extreme_design <- tibble(
  param_name  = "cpy_extreme",
  change_type = c("plus", "minus"),
  new_value   = c(1, 2)  
)

cpy_by_stage_design <- expand.grid(
  stage       = c("I", "II", "III", "IV"),
  change_type = c("minus", "plus"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    param_name = paste0("stage_cost_", stage),
    new_value  = ifelse(change_type == "minus", 0.75, 1.25)
  ) %>%
  select(param_name, change_type, new_value)

# 3F. Additional cost and survival designs
util_design <- two_rows_unbounded("util_mult", 1, pct = 0.25)

cost_design <- tribble(
  ~param_name,               ~change_type, ~new_value,
  "fp_incidental_cost_pct",  "minus",      -0.25,
  "fp_incidental_cost_pct",  "plus",        0.25,
  "invasive_cost_pct",       "minus",      -0.25,
  "invasive_cost_pct",       "plus",        0.25
)

surv_cpy_design <- bind_rows(
  two_rows_unbounded("survival_mult",         1, pct = 0.25),
  two_rows_unbounded("cost_per_year_mult",    1, pct = 0.25)
)

surv_design <- bind_rows(
  two_rows_unbounded("survival_mult",         1, pct = 0.25)
)

mortality_design <- bind_rows(
  two_rows_unbounded("all_cause_mortality_mult",    1, pct = 0.25),
  two_rows_unbounded("total_cancer_mortality_mult", 1, pct = 0.25)
)

# 3G. Cancer types that get USPSTF screening
screen_cancers <- c("lung and bronchus", "colon and rectum",
                    "breast", "prostate", "cervix")


dsa_design <- bind_rows(
  scalar_design, 
  dwell_design, 
  screen_design, 
  mced_fp_design, 
  cpy_by_stage_design,
  cost_design, 
  surv_design, 
  mortality_design, 
  sens_design, 
  spec_design, 
  util_design
) |>
  mutate(
    run_id = row_number(),
    label  = paste0(param_name, "_", change_type)
  )

param_bounds <- dsa_design |>
  filter(!is.na(new_value)) |>
  group_by(param_name) |>
  summarise(
    min    = min(new_value),   # “minus” value
    max    = max(new_value),   # “plus”  value
    .groups = "drop"
  )

apply_change <- function(ip, dl, sens, spec,
                         param_name, value,
                         screen_cancers, def_param_list) {
  
  #------ 1. Sensitivity / specificity multipliers ---------------------------
  if (param_name == "all_sens") {
    sens <- sens |>
      mutate(Sensitivity = pmin(0.99, pmax(0.0000000001, Sensitivity * value)))
  }
  
  if (param_name == "all_spec") {
    spec <- spec |>
      mutate(Specificity = pmin(0.9999999999, pmax(0.0000000001, Specificity * value)))
  }
  
  #------ 2. Cost‑per‑year by stage -----------------------------------------
  if (grepl("^stage_cost_", param_name)) {
    target_stage <- sub("^stage_cost_", "", param_name)
    dl$cost_per_year <- cpy_mult_one_stage(
      df           = dl$cost_per_year,
      target_stage = target_stage,
      mult         = value
    )
  }
  
  #------ 3. Scalar parameters ----------------------------------------------
  if (param_name %in% c("mced_compliance_rate", "overdiagnosis_factor",
                        "mced_test_cost", "fp_track_num",
                        "mced_fp_compliance", "fp_disutility_val",
                        "d_c", "d_e")) {
    ip[[param_name]] <- value
  }
  
  #------ 4. Compliance frequency multiplier --------------------------------
  if (param_name == "Compliance_frequency") {
    dl$screening_fp_diag_costs <- compliance_freq_mult(
      dl$screening_fp_diag_costs,
      pct            = value,          # already a ±% value (e.g. –0.25, +0.25)
      screen_cancers = screen_cancers
    )
  }
  
  #------ 5. Invasive‑test multiplier ---------------------------------------
  if (param_name == "invasive_multi_mult") {
    ip$invasive_multi <- value
  }
  
  #------ 6. Dwell‑time scenario (discrete 1,2,3) ---------------------------
  if (param_name == "dwell_scenario") {
    ip$dwell_scenario <- paste0("Scenario_", value)
    ip$ss_scenario    <- c(`1` = 1, `2` = 3, `3` = 5)[as.character(value)]
  }
  
  #------ 7. Screening‑cost table tweaks ------------------------------------
  if (param_name == "screen_cost_pct") {
    dl$screening_fp_diag_costs <- screen_cost_pct(
      dl$screening_fp_diag_costs, pct = value)
  }
  
  if (param_name == "fp_incidental_cost_pct") {
    dl$screening_fp_diag_costs <- dl$screening_fp_diag_costs |>
      cost_component_pct(value, "FP_cost") |>
      cost_component_pct(value, "incidental_cost")
  }
  
  if (param_name == "invasive_cost_pct" && !dplyr::near(value, 0)) {
    dl$screening_fp_diag_costs <-
      cost_component_pct(dl$screening_fp_diag_costs, value, "invasive_cost")
  }
  
  #------ 8. Utility / survival / cpy multipliers ---------------------------
  if (param_name == "util_mult")        dl$utility_long        <- util_mult(dl$utility_long,  value)
  if (param_name == "survival_mult")    dl$transposed_survival <- surv_mult(dl$transposed_survival, value)
  # if (param_name == "cost_per_year_mult") dl$cost_per_year     <- cpy_mult(dl$cost_per_year,  value)
  
  #------ 9. Re‑derive expanded utilities if d_e changed --------------------
  if (param_name == "d_e") {
    dl$expanded_utility <- update_utility(
      d_e             = value,
      raw_data_folder = "Data/raw",
      female_prop     = def_param_list$female_prop
    )
  }
  
  #------10. Mortality multipliers ------------------------------------------
  if (param_name == "all_cause_mortality_mult") {
    dl$all_cancer_mort <- adjust_all_cause_mort(
      all_cause_raw      = read_csv("Data/raw/all_cause_mortality_40.csv"),
      cancer_mort        = read_csv("Data/raw/cancer_mort.csv"),
      cancer_list        = cancer_types,
      mult_all_cause     = value,
      mult_total_cancer  = 1)
  }
  
  if (param_name == "total_cancer_mortality_mult") {
    dl$all_cancer_mort <- adjust_all_cause_mort(
      all_cause_raw      = read_csv("Data/raw/all_cause_mortality_40.csv"),
      cancer_mort        = read_csv("Data/raw/cancer_mort.csv"),
      cancer_list        = cancer_types,
      mult_all_cause     = 1,
      mult_total_cancer  = value)
  }
  
  #---------------------------------------------------------------------------
  list(ip = ip, dl = dl, sens = sens, spec = spec)
}



psa_design <- tidyr::crossing(run_id = seq_len(n_psa), param_bounds) %>%
  group_by(run_id, param_name) %>%
  
  # for each run_id & param, draw exactly one sample
  summarise(
    sample_value = switch(
      param_name,
      "dwell_scenario" = sample(1:3, 1),
      "fp_track_num"  = round(runif(1, min, max)),
      runif(1, min, max)
    ),
    .groups = "drop"
  )


run_psa_one <- function(sample_df,
                        base_params, base_data,
                        df_sensitivity, original_specificity,
                        sensitivity_matrix, def_param_list,
                        cancer_types_overall) {
  
  ip  <- base_params
  dl   <- base_data
  sens  <- df_sensitivity
  spec  <- original_specificity
  
  # iteratively apply each parameter draw
  for (i in seq_len(nrow(sample_df))) {
    res <- apply_change(ip, dl, sens, spec,
                        param_name = sample_df$param_name[i],
                        value = sample_df$sample_value[i],
                        screen_cancers = screen_cancers,
                        def_param_list = def_param_list)
    ip <- res$ip
    dl <- res$dl
    sens <- res$sens
    spec <- res$spec
  }
  
  # ----- run the model exactly once with the full draw -----
  out <- run_usual_model(
    cancer_list = cancer_types_overall,
    data_list   = dl,
    def_param_list = def_param_list,
    input_params    = ip,
    original_sensitivity = sens,
    specificity_data     = spec,
    Values        = sensitivity_matrix,
    cancer_types_overall = cancer_types_overall)
  
  keep <- c("Incremental_Cost",  "Incremental_Utility", "ICER",
            "mced_cost", "mced_utility",
            "usual_cost", "usual_utility")
  
  icer_row <- out$final_result_list$Overall$final_icer_result |>
    dplyr::filter(Variable %in% keep) |>
    dplyr::select(Variable, Value) |> 
    tidyr::pivot_wider(names_from = Variable,
                       values_from = Value) |>
    dplyr::mutate(run_id = unique(sample_df$run_id)) |>
    dplyr::select(run_id, dplyr::all_of(keep))  
}

plan(multisession, workers = parallel::detectCores() - 1)

psa_results <- split(psa_design, psa_design$run_id)

psa_out <- psa_results |>                     # → list of 10 000 small dfs
  furrr::future_map_dfr(
    run_psa_one,
    base_params         = input_params,
    base_data           = data_loaded,
    df_sensitivity      = df_sensitivity,
    original_specificity= original_specificity,
    sensitivity_matrix  = sensitivity_matrix,
    def_param_list      = def_param_list,
    cancer_types_overall= cancer_types_overall,
    .progress = TRUE
  ) |>
  dplyr::bind_rows()

write_csv(psa_out, file.path(output_dir, "psa_results.csv"))


willingness_to_pay <- seq(0, 200000, 1000)
ceac <- purrr::map_dfr(willingness_to_pay, ~
                         collapsed |>
                         summarise(p_accept = mean(ICER < .x)) |>
                         mutate(wtp = .x)
)
write_csv(ceac, file.path(output_dir, "ceac_results.csv"))
