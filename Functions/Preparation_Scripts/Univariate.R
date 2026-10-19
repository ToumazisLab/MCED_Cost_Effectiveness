##############################################################################
# 1. HELPER FUNCTIONS FOR BUILDING TIBBLES (BOUNDED/UNBOUNDED ROWS)
##############################################################################
source("Global.R")
source("Functions/Preparation_Scripts/Def_params.R", echo=FALSE)
source("Functions/Preparation_Scripts/Source_file.R", echo=FALSE)
source("Functions/Run_Univariate_Model.R", echo=TRUE)

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
  active_slip = FALSE,
  two_step = FALSE,
  overall_sens = 0,
  overall_specificity = 0,
  mced_secondary_test_cost = 0
)

data_loaded <- load_data(input_params, def_param_list)

two_rows_bounded <- function(name, base, pct = 0.25) {
  tibble(
    param_name  = name,
    change_type = c("minus","plus"),
    new_value   = c(base*(1-pct), base*(1+pct)) |>
      pmin(1) |> pmax(0)  # clamp to [0,1]
  )
}

two_rows_unbounded <- function(name, base, pct = 0.25) {
  tibble(
    param_name  = name,
    change_type = c("minus","plus"),
    new_value   = c(base*(1-pct), base*(1+pct))
  )
}

##############################################################################
# 2. HELPER FUNCTIONS FOR TRANSFORMING/UPDATING DATA
##############################################################################

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
  new_value   = c(0.95,       0.9999)
)

# 3B. Example: scalar params
mced_fp_design <- tibble(
  param_name  = c("mced_fp_compliance", "mced_fp_compliance"),
  change_type = c("minus", "plus"),
  new_value   = c((input_params$mced_fp_compliance * 0.75), 1)
)

scalar_design <- bind_rows(
  two_rows_bounded("mced_compliance_rate", input_params$mced_compliance_rate),
  
  # tibble(
  #   param_name = "mced_compliance_rate",
  #   change_type = c("minus", "plus"),
  #   new_value = c(0.5, 1)
  # ), 
  
  two_rows_bounded("overdiagnosis_factor", input_params$overdiagnosis_factor),
  
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
    new_value   = c(0.075, 0.125)
  ),
  
  two_rows_unbounded("fp_disutility_val", input_params$fp_disutility_val),
  two_rows_bounded("d_c",                 input_params$d_c),
  two_rows_bounded("d_e",                 input_params$d_e)
)

invasive_mutli_design <- tibble(
  param_name  = "invasive_multi_mult",
  change_type = c("minus", "plus"),
  new_value   = c(0.075, 0.125)
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

eligible_cancers <- data_loaded$screening_fp_diag_costs %>%
  group_by(Cancer) %>%
  summarise(max_inv = max(Invasive_procedure, na.rm = TRUE)) %>%
  filter(max_inv > 0) %>%
  pull(Cancer)

# Build DSA rows for each eligible cancer
inv_cancer_design <- expand.grid(
  Cancer      = eligible_cancers,
  change_type = c("minus", "plus"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    param_name = paste0("invasive_rate_", Cancer),
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

##############################################################################
# 4. COMBINE ALL DESIGNS INTO A SINGLE DSA DESIGN TABLE
##############################################################################

dsa_design <- bind_rows(
  sens_design,
  spec_design,
  scalar_design,
  dwell_design,
  screen_design,
  mced_fp_design,
  cpy_by_stage_design,
  cost_design,
  surv_cpy_design,
  mortality_design,
  util_design,
  inv_cancer_design
) |>
  mutate(
    run_id = row_number(),
    label  = paste0(param_name, "_", change_type)
  )


##############################################################################
# 5. MAIN RUN FUNCTION
##############################################################################



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
  
  # 5A. Copy base params/data locally
  ip <- base_params
  dl <- base_data
  get_scalar <- function(x) as.numeric(x[[1]])
  
  # 5B. Prepare local versions of sens/spec
  multiplier <- as.numeric(design_row$new_value)
  local_sens <- df_sensitivity
  local_spec <- original_specificity
  
  # 5C. Adjust sensitivity if needed
  if (design_row$label == "all_sens_plus" ||
      design_row$label == "all_sens_minus") {
    
    local_sens <- local_sens %>%
      mutate(
        Sensitivity = Sensitivity * multiplier,
        Sensitivity = pmin(0.9999, pmax(0, Sensitivity))
      )
  }
  
  # 5D. Adjust specificity if needed
  if (design_row$label == "all_spec_minus" ||
      design_row$label == "all_spec_plus") {
    
    local_spec <- local_spec %>%
      mutate(
        Specificity = multiplier,
        Specificity = pmin(0.9999, pmax(0, Specificity))
      )
  }
  
  # 5E. Adjust cost per year by stage
  if (grepl("^stage_cost_", design_row$param_name)) {
    target_stage <- sub("stage_cost_", "", design_row$param_name)
    dl$cost_per_year <- cpy_mult_one_stage(
      df = dl$cost_per_year,
      target_stage = target_stage,
      mult = multiplier
    )
  }
  
  if (startsWith(design_row$param_name, "invasive_rate_")) {
    target_cancer <- sub("invasive_rate_", "", design_row$param_name)
    
    # Only modify if this cancer is in the data
    dl$screening_fp_diag_costs <- dl$screening_fp_diag_costs %>%
      mutate(
        Invasive_procedure = if_else(
          Cancer == target_cancer & Invasive_procedure > 0,
          Invasive_procedure * get_scalar(design_row$new_value),
          Invasive_procedure
        )
      )
  }
  
  # 5F. Modify scalar params directly
  if (design_row$param_name %in%
      c("mced_compliance_rate", "overdiagnosis_factor", "mced_test_cost",
        "fp_track_num", "mced_fp_compliance", "fp_disutility_val",
        "d_c", "d_e")) {
    ip[[design_row$param_name]] <- get_scalar(design_row$new_value)
  }
  
  # 5G. Compliance frequency multiplier (for screening)
  if (design_row$param_name == "Compliance_frequency") {
    dl$screening_fp_diag_costs <- compliance_freq_mult(
      dl$screening_fp_diag_costs,
      get_scalar(design_row$new_value),
      screen_cancers
    )
  }
  
  # 5H. Extreme or by-stage cost per year
  if (design_row$param_name == "cpy_extreme") {
    type <- if (design_row$new_value == 1) "plus" else "minus"
    dl$cost_per_year <- cpy_extreme_by_cancer(
      df   = dl$cost_per_year,
      type = type
    )
  }
  
  # 5I. Invasive multiplier
  if (design_row$param_name == "invasive_multi_mult") {
    ip$invasive_multi <- get_scalar(design_row$new_value)
  }
  
  # 5J. Dwell scenario
  if (design_row$param_name == "dwell_scenario") {
    if (design_row$new_value == 1) {
      ip$dwell_scenario <- "Scenario_1"
      ip$ss_scenario    <- 1
    } else if (design_row$new_value == 2) {
      ip$dwell_scenario <- "Scenario_3"
      ip$ss_scenario    <- 5
    } else {
      ip$dwell_scenario <- "Scenario_2"
      ip$ss_scenario    <- 3
    }
  }
  
  # 5K. Adjust screening-cost table
  if (design_row$param_name == "screen_cost_pct") {
    dl$screening_fp_diag_costs <- screen_cost_pct(dl$screening_fp_diag_costs,
                                                  get_scalar(design_row$new_value))
  }
  
  if (design_row$param_name == "fp_incidental_cost_pct") {
    dl$screening_fp_diag_costs <- dl$screening_fp_diag_costs |>
      cost_component_pct(get_scalar(design_row$new_value), "FP_cost") |>
      cost_component_pct(get_scalar(design_row$new_value), "incidental_cost")
  }
  
  if (design_row$param_name == "invasive_cost_pct") {
    dl$screening_fp_diag_costs <-
      cost_component_pct(dl$screening_fp_diag_costs,
                         get_scalar(design_row$new_value),
                         "invasive_cost")
  }
  
  # 5L. Adjust utilities/survival/cost-per-year
  if (design_row$param_name == "util_mult") {
    dl$utility_long <- util_mult(dl$utility_long, get_scalar(design_row$new_value))
  }
  
  if (design_row$param_name == "survival_mult") {
    dl$transposed_survival <- surv_mult(dl$transposed_survival, get_scalar(design_row$new_value))
  }
  
  if (design_row$param_name == "cost_per_year_mult") {
    dl$cost_per_year <- cpy_mult(dl$cost_per_year, get_scalar(design_row$new_value))
  }
  
  # 5M. If d_e changed, recalc expanded utility
  if (design_row$param_name == "d_e") {
    dl$expanded_utility <- update_utility(
      d_e            = get_scalar(design_row$new_value),
      over_all_data_folder = "Data/overall",
      female_prop    = def_param_list$female_prop
    )
  }
  
  # 5N. Adjust mortality
  if (design_row$param_name == "all_cause_mortality_mult") {
    dl$all_cancer_mort <- adjust_all_cause_mort(
      all_cause_raw = read_csv("Data/overall/all_cause_mortality_40.csv", show_col_types = FALSE),
      cancer_mort   = read_csv("Data/overall/cancer_mort.csv", show_col_types = FALSE),
      cancer_list   = cancer_types,
      mult_all_cause = get_scalar(design_row$new_value),
      mult_total_cancer = 1  # leave cancer as is
    )
  }
  
  if (design_row$param_name == "total_cancer_mortality_mult") {
    dl$all_cancer_mort <- adjust_all_cause_mort(
      all_cause_raw    = read_csv("Data/overall/all_cause_mortality_40.csv", show_col_types = FALSE),
      cancer_mort      = read_csv("Data/overall/cancer_mort.csv", show_col_types = FALSE),
      cancer_list      = cancer_types,
      mult_all_cause   = 1,
      mult_total_cancer= get_scalar(design_row$new_value)
    )
  }
  
  # 5O. RUN THE MODEL
  res <- run_usual_model(
    cancer_list          = cancer_types_overall,
    data_list            = dl,
    def_param_list       = def_param_list,
    input_params         = ip,
    original_sensitivity = local_sens,
    specificity_data     = local_spec,
    Values               = sensitivity_matrix,
    cancer_types_overall = cancer_types_overall
  )
  
  # 5P. EXTRACT & RETURN RESULTS
  icer_df <- res$final_result_list$Overall$final_icer_result
  
  out_tibble <- tibble(
    run_id              = design_row$run_id,
    param_name          = design_row$param_name,
    change_type         = design_row$change_type,
    new_value           = design_row$new_value,
    Incremental_Cost    = icer_df$Value[icer_df$Variable=="Incremental_Cost"],
    Incremental_Utility = icer_df$Value[icer_df$Variable=="Incremental_Utility"],
    ICER                = icer_df$Value[icer_df$Variable=="ICER"],
    Mced_cost = icer_df$Value[icer_df$Variable == "mced_cost"],
    Mced_utility = icer_df$Value[icer_df$Variable == "mced_utility"],
    Usual_cost = icer_df$Value[icer_df$Variable == "usual_cost"],
    Usual_utility = icer_df$Value[icer_df$Variable == "usual_utility"]
  )
  
  message(
    "\n---> Results for run ",
    design_row$run_id,
    " (", design_row$param_name, " – ", design_row$change_type, "):\n"
  )
  print(out_tibble)
  
  # 5Q. SAVE RESULTS
  out_rds_file <- file.path(
    result_dir,
    paste0("final_results_", design_row$run_id, ".rds")
  )
  # saveRDS(res$final_result_list, out_rds_file)
  
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

data_loaded <- load_data(input_params, def_param_list)

sens_cancers <- cancer_types_overall <- cancer_types
# data_loaded <- load_data(input_params, def_param_list)
df_sensitivity <- read_csv("Data/overall/sensitivity_long_overall.csv")
sensitivity_results <- process_sensitivity(
  cancer_list      = sens_cancers,
  sensitivity_long = df_sensitivity,
  seer_draw        = data_loaded$seer_draw,
  xStage           = def_param_list$xStage,
  active_slip      = TRUE,
  dwell_model_df   = data_loaded$dwell_model_all_cancers
)

original_specificity <- data.frame(
  Cancer = cancer_types,
  # Default starting value for specificity
  Specificity = rep(0.98, 16)
)
sensitivity_processed <- sensitivity_results$sensitivity_long
sensitivity_matrix <- sensitivity_results$Values

results <- run_sensitivity_analysis(
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


baseline <- run_usual_model(
  cancer_list          = cancer_types_overall,
  data_list            = data_loaded,
  def_param_list       = def_param_list,
  input_params         = input_params,
  original_sensitivity = df_sensitivity,
  specificity_data     = original_specificity,
  Values               = sensitivity_matrix,
  cancer_types_overall = cancer_types_overall
)


make_tornado <- function(df, value_col, baseline, xlab) {
  
  # ---------------- data prep (unchanged) ---------------------------
  df2 <- df %>%
    filter(param_name != "baseline") %>%
    group_by(plot_label) %>%
    filter(max(.data[[value_col]]) != min(.data[[value_col]])) %>%
    ungroup()
  
  reorder_order <- df2 %>%
    group_by(plot_label) %>%
    summarise(span = max(.data[[value_col]]) - min(.data[[value_col]]),
              .groups = "drop") %>%
    arrange(desc(abs(span))) %>%
    pull(plot_label)
  
  df2 <- df2 %>%
    mutate(plot_label = factor(plot_label, levels = rev(reorder_order)))
  
  ranges <- df2 %>%
    group_by(plot_label) %>%
    summarise(
      xmin = min(.data[[value_col]]),
      xmax = max(.data[[value_col]]),
      .groups = "drop"
    )
  
  endpoints <- df2 %>%
    select(plot_label, plot_new_value, value = !!sym(value_col)) %>%
    left_join(ranges, by = "plot_label") %>%
    mutate(
      hjust      = ifelse(value == xmin, 1.02, -0.02),  # keep labels inside panel
      label_text = plot_new_value
    )
  
  # ---------------- dynamic x-limits --------------------------------
  min_x <- min(ranges$xmin)
  max_x <- max(ranges$xmax)
  pad   <- 0.2 * (max_x - min_x)              # a bit more padding
  xlims <- c(min_x - pad, max_x + pad)
  
  # ---------------- build plot --------------------------------------
  p <- ggplot(ranges, aes(y = plot_label)) +
    geom_segment(aes(x = xmin, xend = xmax, yend = plot_label),
                 linewidth = 8, colour = "#69b3a2") +
    geom_point(data = endpoints, aes(x = value),
               size = 1.5, colour = "grey20") +
    geom_text(data = endpoints,
              aes(x = value, label = label_text, hjust = hjust),
              size = 6) +
    geom_vline(xintercept = baseline,
               linetype = "dashed",
               linewidth = 0.7) +
    labs(x = xlab, y = NULL) +
    coord_cartesian(xlim = xlims, clip = "off") +
    theme_minimal(base_size = 20) +
    theme(
      # remove all gridlines
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      # y-axis labels with a bit of space
      axis.text.y = element_text(
        size = 14,
        margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 5)
      ),
      # larger margins so labels outside panel are visible
      plot.margin = ggplot2::margin(t = 10, r = 30, b = 10, l = 20)
    )
  
  # --------- conditional axis-label formatting ---------------------
  if (grepl("Utility", xlab, ignore.case = TRUE)) {
    p <- p + scale_x_continuous(
      limits = xlims,
      labels = scales::comma,
      expand = expansion(mult = 0)
    )
  } else {
    p <- p + scale_x_continuous(
      limits = xlims,
      labels = scales::dollar_format(prefix = "$", big.mark = ","),
      expand = expansion(mult = 0)
    )
  }
  
  p
}




plot_new_labels <- read_csv("../tornado_plot_new_label.csv")
param_labels <- c(
  "mced_compliance_rate"         = "MCD Compliance (0.75)",
  "overdiagnosis_factor"         = "Overdiagnosis Factor (0.05)",
  "mced_test_cost"               = "MCD Test Cost ($500)",
  "fp_track_num"                 = "Number of MCD FP Follow-ups (3)",
  "invasive_multi_mult"          = "Invasive Procedure After Two FPs (0.1)",
  "fp_disutility_val"            = "FP Disutility (0.0125)",
  "d_c"                          = "Discount Rate (Cost, 0.03)",
  "d_e"                          = "Discount Rate (Utility, 0.03)",
  "dwell_scenario"               = "Dwell Time Scenario (Average)",
  "screen_cost_pct"              = "UC Screening Cost (Physician Fee Schedule)",
  "mced_fp_compliance"           = "Compliance to UC After MCD (1)",
  "stage_cost_I"                 = "Post-Dx Cost (Stage I) (Reddy 2022)",
  "stage_cost_II"                = "Post-Dx Cost (Stage II) (Reddy 2022)",
  "stage_cost_III"               = "Post-Dx Cost (Stage III) (Reddy 2022)",
  "stage_cost_IV"                = "Post-Dx Cost (Stage IV) (Reddy 2022)",
  "fp_incidental_cost_pct"       = "Follow-up & Incidental Imaging (Physician Fee Schedule)",
  "invasive_cost_pct"            = "UC Invasive Biopsy Cost (Physician Fee Schedule)",
  "survival_mult"                = "Post-Dx Survival (SEER 17 2021)",
  "cost_per_year_mult"           = "Post-Dx Cost (All Stages) (Reddy 2022)",
  "all_cause_mortality_mult"     = "All-Cause Mortality (US Life Table 2021)",
  "total_cancer_mortality_mult"  = "Cancer-Specific Mortality (SEER 17 2021)",
  "all_sens"                     = "Sensitivity (All) (eFigure 2)",
  "all_spec"                     = "Specificity (All) (0.98)",
  "util_mult"                    = "Post-Dx Disutility (Tafazzoli 2021)",
  "invasive_rate_breast"            = "Invasive Rate (Breast) (0.033)",
  "invasive_rate_cervix"            = "Invasive Rate (Cervix) (0.1)",
  "invasive_rate_lung and bronchus" = "Invasive Rate (Lung and Bronchus) (0.03)",
  "invasive_rate_prostate"          = "Invasive Rate (Prostate) (0.2)",
  "invasive_rate_bladder"                          = "Invasive Rate (Bladder) (0.1)",
  "invasive_rate_head and neck"                    = "Invasive Rate (Head and Neck) (0.1)",
  "invasive_rate_liver and intrahepatic bile duct" = "Invasive Rate (Liver and Intrahepatic Bile Duct) (0.1)",
  "invasive_rate_lymphoma"                         = "Invasive Rate (Lymphoma) (0.1)",
  "invasive_rate_pancreas"                         = "Invasive Rate (Pancreas) (0.1)",
  "invasive_rate_ovarian"                          = "Invasive Rate (Ovarian) (0.1)",
  "invasive_rate_uterus"                           = "Invasive Rate (Uterus) (0.1)",
  "invasive_rate_kidney and renal pelvis" = "Invasive Rate (Kidney) (0.1)"
)

clean_param_labels <- c(
  "mced_compliance_rate"         = "MCD Compliance",
  "overdiagnosis_factor"         = "Overdiagnosis Factor",
  "mced_test_cost"               = "MCD Test Cost",
  "fp_track_num"                 = "Number of MCD FP Follow-ups",
  "invasive_multi_mult"          = "Invasive Procedure After Two FPs",
  "fp_disutility_val"            = "FP Disutility",
  "d_c"                          = "Discount Rate (Cost)",
  "d_e"                          = "Discount Rate (Utility)",
  "dwell_scenario"               = "Dwell Time Scenario",
  "screen_cost_pct"              = "UC Screening Cost",
  "mced_fp_compliance"           = "Compliance to UC After MCD",
  "stage_cost_I"                 = "Post-Dx Cost (Stage I)",
  "stage_cost_II"                = "Post-Dx Cost (Stage II)",
  "stage_cost_III"               = "Post-Dx Cost (Stage III)",
  "stage_cost_IV"                = "Post-Dx Cost (Stage IV)",
  "fp_incidental_cost_pct"       = "Follow-up & Incidental Imaging",
  "invasive_cost_pct"            = "UC Invasive Biopsy Cost",
  "survival_mult"                = "Post-Dx Survival",
  "all_cause_mortality_mult"     = "All-Cause Mortality",
  "total_cancer_mortality_mult"  = "Cancer-Specific Mortality",
  "all_sens"                     = "Sensitivity (All)",
  "all_spec"                     = "Specificity (All)",
  "util_mult"                    = "Post-Dx Disutility",
  "invasive_rate_breast"            = "Invasive Rate (Breast)",
  "invasive_rate_cervix"            = "Invasive Rate (Cervix)",
  "invasive_rate_lung and bronchus" = "Invasive Rate (Lung and Bronchus)",
  "invasive_rate_prostate"          = "Invasive Rate (Prostate)",
  "invasive_rate_bladder"                          = "Invasive Rate (Bladder)",
  "invasive_rate_head and neck"                    = "Invasive Rate (Head and Neck)",
  "invasive_rate_liver and intrahepatic bile duct" = "Invasive Rate (Liver and Intrahepatic Bile Duct)",
  "invasive_rate_lymphoma"                         = "Invasive Rate (Lymphoma)",
  "invasive_rate_pancreas"                         = "Invasive Rate (Pancreas)",
  "invasive_rate_ovarian"                          = "Invasive Rate (Ovarian)",
  "invasive_rate_uterus"                           = "Invasive Rate (Uterus)"
)


results_updated <- results %>%
  # 1) bring in old_results’ new_value and plot_new_value by the two keys
  left_join(
    plot_new_labels %>% 
      select(param_name, change_type, 
             new_value_old = new_value, 
             plot_new_value),
    by = c("param_name", "change_type")
  ) %>%
  # 2) overwrite new_value when there’s an old one, and keep plot_new_value
  mutate(
    new_value = coalesce(new_value_old, new_value)
  ) %>%
  # 3) drop the helper column
  select(-new_value_old)

results <- results_updated

dsa_results <- results %>%
  mutate(plot_label = param_labels[param_name]) %>% 
  mutate(`Incremental_Cost_Mil` = Incremental_Cost / 1e6)


ordered <- dsa_results %>% 
  group_by(plot_label) %>% 
  summarise(span = max(Incremental_Utility) - min(Incremental_Utility), .groups = "drop") %>% 
  arrange(desc(span))


baseline_icer_results <- baseline$final_result_list$Overall$final_icer_result

baseline_icer <- baseline_icer_results$Value[baseline_icer_results$Variable == "ICER"]
baseline_cost <- baseline_icer_results$Value[baseline_icer_results$Variable == "Incremental_Cost"]
baseline_utility <- baseline_icer_results$Value[baseline_icer_results$Variable == "Incremental_Utility"]

# ICER plot
p_icer <- make_tornado(
  df        = dsa_results,
  value_col = "ICER",
  baseline  = baseline_icer,
  xlab      = "ICER"
)
print(p_icer)

# Incremental Cost plot
p_cost <- make_tornado(
  df        = dsa_results,
  value_col = "Incremental_Cost_Mil",
  baseline  = (baseline_cost/1e6),
  xlab      = "Incremental Cost (Millions, $USD)"
)
print(p_cost)

# Incremental Utility plot
p_util <- make_tornado(
  df        = dsa_results,
  value_col = "Incremental_Utility",
  baseline  = baseline_utility,
  xlab      = "Incremental Utility"
)
print(p_util)


