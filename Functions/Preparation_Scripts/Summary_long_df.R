generate_full_long_df <- function(final_result_list, cancer_list, m_cancer_deaths, uc_cancer_deaths, final_screen_df) {
  # --- Deaths ---
  m_deaths <- aggregate(Cancer_dead_cases ~ Cancer, data = m_cancer_deaths, sum) %>%
    mutate(Arm = "m", Variable = "Cancer_dead_cases") %>% rename(Value = Cancer_dead_cases)
  
  uc_deaths <- aggregate(Cancer_dead_cases ~ Cancer, data = uc_cancer_deaths, sum) %>%
    mutate(Arm = "uc", Variable = "Cancer_dead_cases") %>% rename(Value = Cancer_dead_cases)
  
  deaths_averted <- merge(uc_deaths, m_deaths, by = "Cancer", suffixes = c("_uc", "_m")) %>%
    mutate(Value = Value_uc - Value_m, Arm = "m_v_uc", Variable = "Cancer_deaths_averted") %>%
    select(Cancer, Arm, Variable, Value)
  
  # --- Overdiagnosis ---
  overdx_result <- lapply(cancer_list, function(ca) {
    df <- final_result_list[["MCED"]][[ca]][["MCED_Cancer_cases"]]
    tibble(
      Cancer = ca,
      overdx_I = sum(df$overdiagnosed_I, na.rm = TRUE),
      overdx_II = sum(df$overdiagnosed_II, na.rm = TRUE),
      overdx_III = sum(df$overdiagnosed_III, na.rm = TRUE),
      overdx_IV = sum(df$overdiagnosed_IV, na.rm = TRUE),
      overdx_Total = sum(df$Total_overdiagnosed, na.rm = TRUE)
    )
  }) %>% bind_rows() %>%
    pivot_longer(-Cancer, names_to = "Variable", values_to = "Value") %>%
    mutate(Arm = "m")
 
  
  # --- FP Cases ---
  fp_case_m <- lapply(cancer_list, function(cancer) {
    mced_data <- final_result_list[["MCED"]][[cancer]][["mced_TN_FP"]]
    m_FP_uc_vars <- intersect(c("FP_i", "FP_FP", "FP_FP_i"), names(mced_data))
    m_uc_vars <- intersect(c("TN_c_FP", "TN_c_FP_i", "nc_FP", "nc_Invasive"), names(mced_data))
    tibble(
      Cancer = cancer,
      Variable = c("FP", "FP_followups", "UC_followups"),
      Value = c(
        sum(mced_data$FP_fu, na.rm = TRUE),
        if (length(m_FP_uc_vars) > 0) sum(unlist(mced_data[m_FP_uc_vars]), na.rm = TRUE) else 0,
        if (length(m_uc_vars) > 0) sum(unlist(mced_data[m_uc_vars]), na.rm = TRUE) else 0
      ),
      Arm = "m"
    )
  }) %>% bind_rows()
  

  library(tidyverse)
  
  # ---------------------------------------------------------------------------
  # 1.  Define the columns that feed each new Variable -------------------------
  # ---------------------------------------------------------------------------
  mced_formulae <- list(
    FP        = "FP_fu",
    FPx2      = c("FP_FP", "FP_i"),
    FPx3      = "FP_FP_i",
    UC_FPx1   = c("TN_c_FP",  "nc_FP"),
    UC_FPx2   = c("nc_Invasive", "TN_c_FP_i")
  )
  
  uc_formulae <- list(
    UC_FPx1     = "FP",
    UC_FPx2   = "Invasive"
  )
  
  # helper to turn a named list of vectors → one tibble row --------------------
  sum_by_formula <- function(data_tbl, formulae) {
    tibble(
      Variable = names(formulae),
      Value    = vapply(formulae, function(cols) {
        present <- intersect(cols, names(data_tbl))
        if (length(present)) sum(data_tbl[present], na.rm = TRUE) else 0
      }, numeric(1))
    )
  }
  
  # ---------------------------------------------------------------------------
  # 2.  Build MCED frame -------------------------------------------------------
  # ---------------------------------------------------------------------------
  fp_case_m <- purrr::map_dfr(cancer_list, function(ca) {
    mced_data <- final_result_list[["MCED"]][[ca]][["mced_TN_FP"]]
    
    sum_by_formula(mced_data, mced_formulae) %>%
      mutate(Cancer = ca, Arm = "m", .before = 1)
  })
  
  # ---------------------------------------------------------------------------
  # 3.  Build Usual-care frame -------------------------------------------------
  # ---------------------------------------------------------------------------
  fp_case_uc <- purrr::map_dfr(cancer_list, function(ca) {
    uc_data <- final_result_list[["usual"]][[ca]][["FP_TN_Screening_Summary"]]
    
    sum_by_formula(uc_data, uc_formulae) %>%
      mutate(Cancer = ca, Arm = "uc", .before = 1)
  })
  

  # --- FP Costs ---
  fp_cost_m <- lapply(cancer_list, function(cancer) {
    mced_data <- final_result_list[["MCED"]][[cancer]][["mced_tn_fp_cost"]]
    tibble(
      Cancer = cancer,
      Variable = c("h_Imaging_cost", "h_Diagnostic_cost", "h_Invasive_cost"),
      Value = c(sum(mced_data$Imaging_cost, na.rm = TRUE), sum(mced_data$Diagnostic_cost, na.rm = TRUE), sum(mced_data$Invasive_cost, na.rm = TRUE)),
      Arm = "m"
    )
  }) %>% bind_rows()

  fp_cost_uc <- lapply(cancer_list, function(cancer) {
    uc_data <- final_result_list[["usual"]][[cancer]][["FP_TN_Screening_Summary"]]
    tibble(
      Cancer = cancer,
      Variable = c("h_Imaging_cost", "h_Diagnostic_cost", "h_Invasive_cost"),
      Value = c(sum(uc_data$Imaging_cost, na.rm = TRUE), sum(uc_data$Diagnostic_cost, na.rm = TRUE), sum(uc_data$Invasive_cost, na.rm = TRUE)),
      Arm = "uc"
    )
  }) %>% bind_rows()
  #####  

  # --- Cancer Costs ---
  cancer_costs_m <- lapply(cancer_list, function(ca) {
    mced_data <- final_result_list[["MCED"]][[ca]][["MCED_Cancer_costs"]]
    tibble(
      Cancer = ca,
      Variable = c("c_Imaging_cost", "c_Diagnostic_cost", "c_Invasive_cost", "mced_dx", "screen_dx", "incidental_dx", "total_pre"),
      Value = c(
        sum(mced_data$total_imaging_per_cycle, na.rm = TRUE),
        sum(mced_data$total_diagnostic_per_cycle, na.rm = TRUE),
        sum(mced_data$total_invasive_per_cycle, na.rm = TRUE),
        sum(mced_data$total_m_dx_per_cycle, na.rm = TRUE),
        sum(mced_data$total_s_dx_per_cycle, na.rm = TRUE),
        sum(mced_data$total_i_dx_per_cycle, na.rm = TRUE),
        sum(mced_data$total_cost_per_cycle, na.rm = TRUE)
      ),
      Arm = "m"
    )
  }) %>% bind_rows()

  cancer_costs_uc <- lapply(cancer_list, function(ca) {
    uc_data <- final_result_list[["usual"]][[ca]][["Cancer_costs"]]
    tibble(
      Cancer = ca,
      Variable = c("c_Imaging_cost", "c_Diagnostic_cost", "c_Invasive_cost", "screen_dx", "incidental_dx", "total_pre"),
      Value = c(
        sum(uc_data$total_imaging_per_cycle, na.rm = TRUE),
        sum(uc_data$total_diagnostic_per_cycle, na.rm = TRUE),
        sum(uc_data$total_invasive_per_cycle, na.rm = TRUE),
        sum(uc_data$total_s_dx_per_cycle, na.rm = TRUE),
        sum(uc_data$total_i_dx_per_cycle, na.rm = TRUE),
        sum(uc_data$total_cost_per_cycle, na.rm = TRUE)
      ),
      Arm = "uc"
    )
  }) %>% bind_rows()

  # --- Cancer Cases (m only) ---
  cancer_m_cases_long <- final_screen_df %>%
    select(Cancer, s_cases_FN, i_cases_FN, s_cases_NC, i_cases_NC, s_cases_TP, i_cases_TP) %>%
    filter(Cancer != "Overall") %>%
    pivot_longer(-Cancer, names_to = "Variable", values_to = "Value") %>%
    mutate(Arm = "m")
  

  # --- Cancer Cases Costs (m) ---
  cancer_m_costs <- lapply(cancer_list, function(ca) {
    df <- final_result_list[["MCED"]][[ca]][["MCED_Cancer_cases"]]
    stages <- c("I", "II", "III", "IV")

    stage_sums <- lapply(stages, function(stage) {
      cols <- paste0(stage, "_post_cost_", c("TP", "FN", "NA"))
      sum(rowSums(df[, cols, drop = FALSE], na.rm = TRUE))
    })

    names(stage_sums) <- paste0(stages, "_post_cost")
    Total_post_cost <- sum(unlist(stage_sums), na.rm = TRUE)
    Total_post_utility <- sum(df$Total_post_utility_TP + df$Total_post_utility_FN + df$Total_post_utility_NA, na.rm = TRUE)

    tibble(
      Cancer = ca,
      Variable = c(names(stage_sums), "Total_post_cost", "Total_post_utility"),
      Value = c(unlist(stage_sums), Total_post_cost, Total_post_utility),
      Arm = "m"
    )
  }) %>% bind_rows()

  # --- Cancer Cases (uc) ---
  cancer_uc_cases <- lapply(cancer_list, function(ca) {
    df <- final_result_list[["usual"]][[ca]][["Cancer_cases_screening"]]
    tibble(
      Cancer = ca,
      Variable = c("s_cases_uc", "i_cases_uc", "I_post_cost", "II_post_cost", "III_post_cost", "IV_post_cost", "Total_post_cost", "Total_post_utility"),
      Value = c(
        sum(df$s_dx_total, na.rm = TRUE),
        sum(df$i_dx_total, na.rm = TRUE),
        sum(df$I_post_cost, na.rm = TRUE),
        sum(df$II_post_cost, na.rm = TRUE),
        sum(df$III_post_cost, na.rm = TRUE),
        sum(df$IV_post_cost, na.rm = TRUE),
        sum(df$Total_post_cost, na.rm = TRUE),
        sum(df$Total_post_utility, na.rm = TRUE)
      ),
      Arm = "uc"
    )
  }) %>% bind_rows()
  
  ##### Combine all #####
  final_long_df <- bind_rows(
    m_deaths, uc_deaths, deaths_averted,
    overdx_result,
    fp_case_m, fp_case_uc,
    fp_cost_m, fp_cost_uc,
    cancer_costs_m, cancer_costs_uc,
    cancer_m_cases_long, cancer_m_costs,
    cancer_uc_cases
  )
  
  return(final_long_df)
}
