# Shiny Plots
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

sum_by_stage_cancer <- function(incident_matrix, stage_levels) {
  
  # Convert to data frame (if it's a matrix)
  df <- as.data.frame(incident_matrix)
  
  # Summation over rows (the cycles), so we get total incidence for each column
  summed <- df %>%
    summarise(across(everything(), ~sum(.x, na.rm = TRUE))) %>%
    pivot_longer(cols = everything(), names_to = "StageCancer", values_to = "Cases") %>%
    # We assume column names are e.g. "I_lung and bronchus", "III_breast", etc.
    separate(StageCancer, into = c("Stage", "Cancer"), sep = "_", extra = "merge") %>%
    mutate(
      Stage = factor(Stage, levels = stage_levels),  # ensure correct stage order if needed
      Cancer = trimws(Cancer)                        # in case there's leftover whitespace
    )
  
  # Return the long data frame
  summed
}

collect_cancer_summaries <- function(final_result_list, selected_cancers) {
  
  # These columns remain the same:
  cancer_cols_mced  <- c("m_dx_total", "s_dx_total", "i_dx_total")
  healthy_cols_mced <- c("FP_fu", "FP_fu_compliant", "FP_FP", "FP_FP_i", 
                         "TN_c_FP", "TN_c_FP_i")
  
  cancer_cols_usual  <- c("m_dx_total", "s_dx_total", "i_dx_total")
  healthy_cols_usual <- c("FP", "Invasive")
  
  results_list <- lapply(seq_along(selected_cancers), function(i) {
    cancer_index <- selected_cancers[i] 
    cancer_name  <- selected_cancers[i]
    
    #------------------------------------
    # 1) Identify the usual screening range
    #------------------------------------
    if (cancer_name %in% names(usual_screening_ranges)) {
      usual_range <- usual_screening_ranges[[cancer_name]]  # e.g. c(40,74)
    } else {
      # The user said "the rest don't have screening" => maybe no coverage:
      usual_range <- c(0,0)
      # or if you want to adopt MCED range for "no standard screening," do:
      # usual_range <- get_mced_range(cancer_name)
    }
    
    #------------------------------------
    # 2) Identify the MCED screening range
    #------------------------------------
    mced_range <- get_mced_range(cancer_name)  # e.g. c(45,80) or c(50,80)
    
    #------------------------------------
    # 3) Extract MCED data frames
    #------------------------------------
    # final_result_list$MCED[[cancer_index]]$MCED_Cancer_cases
    # final_result_list$MCED[[cancer_index]]$mced_TN_FP
    
    mced_cases_df <- final_result_list$MCED[[ cancer_index ]]$MCED_Cancer_cases
    mced_tn_fp_df <- final_result_list$MCED[[ cancer_index ]]$mced_TN_FP
    
    #------------------------------------
    # 4) Extract Usual data frames
    #------------------------------------
    # final_result_list$usual[[cancer_index]]$Cancer_cases_screening
    # final_result_list$usual[[cancer_index]]$FP_TN_Screening_Summary
    
    usual_cases_df <- final_result_list$usual[[ cancer_index ]]$Cancer_cases_screening
    usual_tn_fp_df <- final_result_list$usual[[ cancer_index ]]$FP_TN_Screening_Summary
    
    #------------------------------------
    # 5) Filter by age range
    #------------------------------------
    # We'll assume each data frame has one row per "cycle" or "row index".
    # If row 1 => age 40, that implies row_i => age = i + 39.
    # We'll add an Age column and subset.
    
    # For MCED cancer metrics:
    if (!is.null(mced_cases_df) && nrow(mced_cases_df) > 0) {
      mced_cases_df$Age <- as.numeric(rownames(mced_cases_df)) + 39
      # Keep only rows in [mced_range[1], mced_range[2]]
      mced_cases_df <- subset(mced_cases_df, 
                              Age >= mced_range[1] & Age <= mced_range[2])
    }
    # For MCED healthy metrics:
    if (!is.null(mced_tn_fp_df) && nrow(mced_tn_fp_df) > 0) {
      mced_tn_fp_df$Age <- as.numeric(rownames(mced_tn_fp_df)) + 39
      mced_tn_fp_df <- subset(mced_tn_fp_df, 
                              Age >= mced_range[1] & Age <= mced_range[2])
    }
    
    # For Usual cancer metrics:
    if (!is.null(usual_cases_df) && nrow(usual_cases_df) > 0) {
      usual_cases_df$Age <- as.numeric(rownames(usual_cases_df)) + 39
      usual_cases_df <- subset(usual_cases_df, 
                               Age >= usual_range[1] & Age <= usual_range[2])
    }
    # For Usual healthy metrics:
    if (!is.null(usual_tn_fp_df) && nrow(usual_tn_fp_df) > 0) {
      usual_tn_fp_df$Age <- as.numeric(rownames(usual_tn_fp_df)) + 39
      usual_tn_fp_df <- subset(usual_tn_fp_df,
                               Age >= usual_range[1] & Age <= usual_range[2])
    }
    
    #------------------------------------
    # 6) Summation or normalization
    #------------------------------------
    # Suppose you want "average per year of screening." 
    # => divide the sums by (range_width).
    
    # Range widths
    usual_years <- max(0, usual_range[2] - usual_range[1] + 1)
    mced_years  <- max(0, mced_range[2]  - mced_range[1]  + 1)
    
    # ---- MCED: Cancer metrics (m_dx_total, s_dx_total, i_dx_total) ----
    if (!is.null(mced_cases_df) && nrow(mced_cases_df) > 0) {
      # Keep only columns that exist
      sub_cols <- intersect(names(mced_cases_df), cancer_cols_mced)
      sums_mced_cancer <- colSums(mced_cases_df[sub_cols], na.rm=TRUE)
      # Normalize by number of MCED screening years if desired:
      sums_mced_cancer <- sums_mced_cancer / ifelse(mced_years>0, mced_years, 1)
      
      df_mced_cancer <- data.frame(
        Cancer   = cancer_name,
        Arm      = "MCED",
        Category = "Cancer",
        Metric   = names(sums_mced_cancer),
        Value    = unname(sums_mced_cancer),
        row.names = NULL
      )
    } else {
      df_mced_cancer <- NULL
    }
    
    # ---- MCED: Healthy metrics (FP/TN) ----
    if (!is.null(mced_tn_fp_df) && nrow(mced_tn_fp_df) > 0) {
      sub_cols <- intersect(names(mced_tn_fp_df), healthy_cols_mced)
      sums_mced_healthy <- colSums(mced_tn_fp_df[sub_cols], na.rm=TRUE)
      sums_mced_healthy <- sums_mced_healthy / ifelse(mced_years>0, mced_years, 1)
      
      df_mced_healthy <- data.frame(
        Cancer   = cancer_name,
        Arm      = "MCED",
        Category = "Healthy",
        Metric   = names(sums_mced_healthy),
        Value    = unname(sums_mced_healthy),
        row.names = NULL
      )
    } else {
      df_mced_healthy <- NULL
    }
    
    # ---- Usual: Cancer metrics ----
    if (!is.null(usual_cases_df) && nrow(usual_cases_df) > 0) {
      sub_cols <- intersect(names(usual_cases_df), cancer_cols_usual)
      sums_usual_cancer <- colSums(usual_cases_df[sub_cols], na.rm=TRUE)
      sums_usual_cancer <- sums_usual_cancer / ifelse(usual_years>0, usual_years, 1)
      
      df_usual_cancer <- data.frame(
        Cancer   = cancer_name,
        Arm      = "Usual",
        Category = "Cancer",
        Metric   = names(sums_usual_cancer),
        Value    = unname(sums_usual_cancer),
        row.names = NULL
      )
    } else {
      df_usual_cancer <- NULL
    }
    
    # ---- Usual: Healthy metrics ----
    if (!is.null(usual_tn_fp_df) && nrow(usual_tn_fp_df) > 0) {
      sub_cols <- intersect(names(usual_tn_fp_df), healthy_cols_usual)
      sums_usual_healthy <- colSums(usual_tn_fp_df[sub_cols], na.rm=TRUE)
      sums_usual_healthy <- sums_usual_healthy / ifelse(usual_years>0, usual_years, 1)
      
      df_usual_healthy <- data.frame(
        Cancer   = cancer_name,
        Arm      = "Usual",
        Category = "Healthy",
        Metric   = names(sums_usual_healthy),
        Value    = unname(sums_usual_healthy),
        row.names = NULL
      )
    } else {
      df_usual_healthy <- NULL
    }
    
    # Combine for this cancer
    do.call(rbind, list(df_mced_cancer, df_mced_healthy,
                        df_usual_cancer, df_usual_healthy))
  })
  
  # Combine across all selected cancers
  out <- do.call(rbind, results_list)
  out
}

# Example function to gather the overall screening/FP/MCED data
collect_overall_summaries <- function(final_result_list) {
  # 1) We'll assume final_result_list$Overall is a list with these data frames:
  #    - FP_TN_Screening_Summary  (for "Usual", "Healthy" => columns: "TN", "FP")
  #    - Cancer_cases_screening   (for "Usual", "Cancer" => columns: "s_dx_total")
  #    - MCED_Cancer_cases_screening (for "MCED", "Cancer" => columns: "m_dx_total","s_dx_total","mced_test_cost")
  #    - mced_TN_FP (for "MCED", "Healthy" => columns: "Total_FP_fu","True_TN","mced_test_cost")
  #
  # We'll sum across rows for each relevant set of columns.
  
  # A helper that sums the columns of interest, then returns a data frame:
  sum_columns <- function(df, cols, arm, category) {
    # Subset only the columns that exist
    existing_cols <- intersect(cols, names(df))
    if (length(existing_cols) == 0) {
      # If none exist, return NULL
      return(NULL)
    }
    # Summation
    sums <- colSums(df[existing_cols], na.rm = TRUE)
    # Build a tidy output
    out <- data.frame(
      Arm      = arm,
      Category = category,
      Metric   = names(sums),
      Value    = unname(sums),
      row.names = NULL
    )
    out
  }
  
  # Grab the "Overall" part
  overall_data <- final_result_list$Overall
  
  # 2) Usual - Healthy
  #    columns: TN, FP
  df_usual_healthy <- sum_columns(
    df       = overall_data$FP_TN_Screening_Summary,
    cols     = c("TN", "FP"),
    arm      = "Usual",
    category = "Healthy"
  )
  
  # 3) Usual - Cancer
  #    columns: s_dx_total
  df_usual_cancer <- sum_columns(
    df       = overall_data$Cancer_cases_screening,
    cols     = "s_dx_total",
    arm      = "Usual",
    category = "Cancer"
  )
  
  # 4) MCED - Cancer
  #    columns: m_dx_total, s_dx_total, mced_test_cost
  df_mced_cancer <- sum_columns(
    df       = overall_data$MCED_Cancer_cases_screening_summary,
    cols     = c("m_dx_total", "s_dx_total", "mced_test_cost", "FN", "TP", "Total_overdiagnosed"),
    arm      = "MCED",
    category = "Cancer"
  )
  
  # 5) MCED - Healthy
  df_mced_healthy <- sum_columns(
    overall_data$mced_TN_FP,
    cols = c("fp_counts", "FP_step1","FP_step2","TN_step1","TN_step2",
             "Total_FP_fu","True_TN","mced_test_cost","fp_counts"),
    arm  = "MCED",
    category = "Healthy"
  )
  
  # 6) Combine them
  final_df <- do.call(rbind, list(
    df_usual_healthy,
    df_usual_cancer,
    df_mced_cancer,
    df_mced_healthy
  ))
  
  final_df
}

get_value <- function(metric, arm = NULL, type = NULL, category = NULL) {
  df <- collect_overall_summaries(run_data$final_result_list) |>
    dplyr::filter(Arm == "MCED")                         # only MCED rows
  
  rows <- df$Metric %in% metric
  if (!is.null(category)) rows <- rows & df$Category == category
  sum(df$Value[rows], na.rm = TRUE)
}

build_fp_levels_df <- function(overall_summary,
                               run_info,
                               def_param_list,
                               cancer_colors) {
  
  level_map <- tribble(
    ~Variable,        ~Arm,        ~Level,  ~Path,
    #  --------------- ---------- ------ ----------------------
    "FP",              "m",       1,  "MCED path",
    "FPx2",            "m",       2,  "MCED path",
    "FPx3",            "m",       3,  "MCED path",
    
    "UC_FPx1",         "m",       2,  "Usual-care follow-up",
    "UC_FPx2",         "m",       3,  "Usual-care follow-up",
    
    "UC_FPx1",         "uc",      1,  "Usual-care path",
    "UC_FPx2",         "uc",      2,  "Usual-care path"
  ) %>% 
    mutate(Level_lbl = factor(Level,
                              levels = 0:3,
                              labels = c("Level 0", "Level 1",
                                         "Level 2", "Level 3")))
  
  
  # ---- slice out the FP variables we care about ---------------------------
  wanted_vars <- c("FP", "FPx2", "FPx3", "UC_FPx1", "UC_FPx2")
  df_core <- overall_summary %>%
    filter(Variable %in% wanted_vars,
           Arm %in% c("m", "uc"))          # ensure arm labels match
  
  # ---- merge with level/path lookup ---------------------------------------
  df_core <- df_core %>%
    left_join(level_map,
              by = c("Variable", "Arm")) %>% 
    mutate(Screening = Cancer %in% def_param_list$screening_cancers) %>% 
    # >>> NEW: cumulative running total within each Cancer × Path  <<<
    group_by(Cancer, Path) %>%
    arrange(Level) %>%                   # Level is numeric (0,1,2,3)
    mutate(Value = cumsum(Value)) %>%    # overwrite with cumulative
    ungroup()
  
  
  # ---- bind everything, add tooltip & colour ------------------------------
  df_all <- df_core %>%
    mutate(
      tooltip = glue::glue("{Cancer}\n{Path}\n{Variable}: {scales::comma(round(Value, 0))}"),
      colour  = cancer_colors[Cancer]
    )
  
  df_all
}

darken_hex <- function(hex_color, factor = 0.8) {
  # factor < 1 => make color darker; factor=0.7 => ~30% darker
  rgb_in <- grDevices::col2rgb(hex_color)
  # multiply each channel
  rgb_out <- pmax(0, pmin(255, rgb_in * factor))
  grDevices::rgb(rgb_out[1], rgb_out[2], rgb_out[3], maxColorValue = 255)
}

lighten_hex <- function(hex_color, factor = 0.2) {
  # factor in [0..1] => how much to move channel toward white
  rgb_in <- grDevices::col2rgb(hex_color)
  rgb_out <- rgb_in + (255 - rgb_in) * factor
  rgb_out <- pmax(0, pmin(255, rgb_out))
  grDevices::rgb(rgb_out[1], rgb_out[2], rgb_out[3], maxColorValue = 255)
}

make_monochrome_palette <- function(base_hex) {
  # --------------------------------------------------
  # 1) Check input
  # --------------------------------------------------
  if (length(base_hex) != 1) {
    stop("`make_monochrome_palette()` expects exactly one hex string (e.g. \"#7D3C98\").")
  }
  if (!grepl("^#([A-Fa-f0-9]{6})$", base_hex)) {
    stop("Please supply a valid 6‐digit hex code, e.g. \"#7D3C98\".")
  }
  
  # --------------------------------------------------
  # 2) Convert that single hex into R, G, B (each 0–1)
  # --------------------------------------------------
  # col2rgb(...) returns a 3×N matrix, with rownames "red","green","blue".
  rgb_mat <- grDevices::col2rgb(base_hex) / 255
  # We know col2rgb(base_hex) has exactly one column, so:
  r <- rgb_mat["red",   1]
  g <- rgb_mat["green", 1]
  b <- rgb_mat["blue",  1]
  
  # --------------------------------------------------
  # 3) Convert (r,g,b) to HSV
  # --------------------------------------------------
  # rgb2hsv expects three numeric vectors in [0,1]; here each is length‐1.
  hsv_vals <- grDevices::rgb2hsv(r, g, b)
  h <- hsv_vals["h"]  # hue ∈ [0,1]
  s <- hsv_vals["s"]  # saturation ∈ [0,1]
  v <- hsv_vals["v"]  # value/brightness ∈ [0,1]
  
  # --------------------------------------------------
  # 4) Compute “lighter” and “darker” values
  # --------------------------------------------------
  # Two lighter shades (capped at 1.0)
  v_lighter1 <- pmin(v + 0.30, 1.0)  # ~30% brighter
  v_lighter2 <- pmin(v + 0.15, 1.0)  # ~15% brighter
  # One darker shade (floored at 0.0)
  v_darker   <- pmax(v - 0.15, 0.0)  # ~15% darker
  
  # --------------------------------------------------
  # 5) Convert each (h,s,v_*) back to hex
  # --------------------------------------------------
  hex_lighter1 <- grDevices::hsv(h = h, s = s, v = v_lighter1)
  hex_lighter2 <- grDevices::hsv(h = h, s = s, v = v_lighter2)
  hex_base     <- base_hex
  hex_darker   <- grDevices::hsv(h = h, s = s, v = v_darker)
  
  # --------------------------------------------------
  # 6) Return a length‐4 character vector:
  #    [1] lighter (~v+30%), [2] lighter (~v+15%), [3] base, [4] darker (~v−15%)
  # --------------------------------------------------
  c(hex_lighter1, hex_lighter2, hex_base, hex_darker)
}



build_mced_sunburst <- function(final_result_list, selected_cancers, final_screen_df_stage, mced_incident_matrix) {
  
  getv <- function(tbl, cols) {
    present <- intersect(cols, names(tbl))
    if (length(present) == 0) return(0)
    sum(sapply(present, function(cl) tbl[[cl]]), na.rm = TRUE)
  }
  
  get_tbl <- function(ca) final_result_list[["MCED"]][[ca]]$MCED_Cancer_cases
  
  lvl1_value <- sum(mced_incident_matrix)
  
  id     <- c("All")
  label  <- c("All Cancers")
  parent <- c("")
  value  <- c(lvl1_value)
  
  stages <- c("I", "II", "III", "IV")
  
  for (ca in selected_cancers) {
    tbl <- get_tbl(ca)
    
    tot_TP_over <- sum(getv(tbl, paste0("TP_o_", stages)))
    tot_FN      <- sum(getv(tbl, paste0("FN_",   stages)))
    tot_noncomp <- sum(getv(tbl, paste0("mced_non_compliant_", stages)))
    tot_cancer  <- tot_TP_over + tot_FN + tot_noncomp
    
    # Cancer node
    id     <- c(id, ca)
    label  <- c(label, ca)
    parent <- c(parent, "All")
    value  <- c(value, tot_cancer)
    
    comp_id  <- paste0(ca, "_comp")
    ncomp_id <- paste0(ca, "_ncomp")
    
    id     <- c(id, comp_id, ncomp_id)
    label  <- c(label, "MCED Compliant", "MCED Non‑compliant")
    parent <- c(parent, ca, ca)
    value  <- c(value, tot_TP_over + tot_FN, tot_noncomp)
    
    # TP/FN branches
    tp_id <- paste0(comp_id, "_TP")
    fn_id <- paste0(comp_id, "_FN")
    
    id     <- c(id, tp_id, fn_id)
    label  <- c(label, "TP", "FN")
    parent <- c(parent, comp_id, comp_id)
    value  <- c(value, tot_TP_over, tot_FN)
    
    fs_tbl <- final_screen_df_stage %>% dplyr::filter(Cancer == ca)
    
    # TP → UC compliant/non-compliant
    uc_tp_id  <- paste0(tp_id, "_UC_comp")
    unc_tp_id <- paste0(tp_id, "_UC_ncomp")
    s_tp <- sum(fs_tbl$s_cases_TP, na.rm = TRUE)
    i_tp <- sum(fs_tbl$i_cases_TP, na.rm = TRUE)
    
    id     <- c(id, uc_tp_id, unc_tp_id)
    label  <- c(label, "UC Compliant", "UC Non-compliant")
    parent <- c(parent, tp_id, tp_id)
    value  <- c(value, s_tp, i_tp)
    
    # TP → UC Compliant → MCED screen-dx → Stage
    mced_dx_tp_id <- paste0(uc_tp_id, "_mced_screen")
    id     <- c(id, mced_dx_tp_id)
    label  <- c(label, "MCED screen-dx")
    parent <- c(parent, uc_tp_id)
    value  <- c(value, s_tp)
    
    for (stg in unique(fs_tbl$Stage)) {
      stg_id <- paste0(mced_dx_tp_id, "_", stg)
      val    <- fs_tbl$s_cases_TP[fs_tbl$Stage == stg]
      id     <- c(id, stg_id)
      label  <- c(label, stg)
      parent <- c(parent, mced_dx_tp_id)
      value  <- c(value, val)
    }
    
    # TP → UC Non-compliant → Incidental dx → Stage
    inc_dx_tp_id <- paste0(unc_tp_id, "_incident")
    id     <- c(id, inc_dx_tp_id)
    label  <- c(label, "Incidental dx")
    parent <- c(parent, unc_tp_id)
    value  <- c(value, i_tp)
    
    for (stg in unique(fs_tbl$Stage)) {
      stg_id <- paste0(inc_dx_tp_id, "_", stg)
      val    <- fs_tbl$i_cases_TP[fs_tbl$Stage == stg]
      id     <- c(id, stg_id)
      label  <- c(label, stg)
      parent <- c(parent, inc_dx_tp_id)
      value  <- c(value, val)
    }
    
    # FN
    uc_fn_id  <- paste0(fn_id, "_UC_comp")
    unc_fn_id <- paste0(fn_id, "_UC_ncomp")
    s_fn <- sum(fs_tbl$s_cases_FN, na.rm = TRUE)
    i_fn <- sum(fs_tbl$i_cases_FN, na.rm = TRUE)
    
    id     <- c(id, uc_fn_id, unc_fn_id)
    label  <- c(label, "UC Compliant", "UC Non-compliant")
    parent <- c(parent, fn_id, fn_id)
    value  <- c(value, s_fn, i_fn)
    
    uc_screen_fn_id <- paste0(uc_fn_id, "_uc_screen")
    id     <- c(id, uc_screen_fn_id)
    label  <- c(label, "UC screen-dx")
    parent <- c(parent, uc_fn_id)
    value  <- c(value, s_fn)
    
    for (stg in unique(fs_tbl$Stage)) {
      stg_id <- paste0(uc_screen_fn_id, "_", stg)
      val    <- fs_tbl$s_cases_FN[fs_tbl$Stage == stg]
      id     <- c(id, stg_id)
      label  <- c(label, stg)
      parent <- c(parent, uc_screen_fn_id)
      value  <- c(value, val)
    }
    
    inc_dx_fn_id <- paste0(unc_fn_id, "_incident")
    id     <- c(id, inc_dx_fn_id)
    label  <- c(label, "Incidental dx")
    parent <- c(parent, unc_fn_id)
    value  <- c(value, i_fn)
    
    for (stg in unique(fs_tbl$Stage)) {
      stg_id <- paste0(inc_dx_fn_id, "_", stg)
      val    <- fs_tbl$i_cases_FN[fs_tbl$Stage == stg]
      id     <- c(id, stg_id)
      label  <- c(label, stg)
      parent <- c(parent, inc_dx_fn_id)
      value  <- c(value, val)
    }
    
    # NC (non-compliant to MCED)
    uc_nc_id  <- paste0(ncomp_id, "_UC_comp")
    unc_nc_id <- paste0(ncomp_id, "_UC_ncomp")
    s_nc <- sum(fs_tbl$s_cases_NC, na.rm = TRUE)
    i_nc <- sum(fs_tbl$i_cases_NC, na.rm = TRUE)
    
    id     <- c(id, uc_nc_id, unc_nc_id)
    label  <- c(label, "UC Compliant", "UC Non-compliant")
    parent <- c(parent, ncomp_id, ncomp_id)
    value  <- c(value, s_nc, i_nc)
    
    uc_screen_nc_id <- paste0(uc_nc_id, "_uc_screen")
    id     <- c(id, uc_screen_nc_id)
    label  <- c(label, "UC screen-dx")
    parent <- c(parent, uc_nc_id)
    value  <- c(value, s_nc)
    
    for (stg in unique(fs_tbl$Stage)) {
      stg_id <- paste0(uc_screen_nc_id, "_", stg)
      val    <- fs_tbl$s_cases_NC[fs_tbl$Stage == stg]
      id     <- c(id, stg_id)
      label  <- c(label, stg)
      parent <- c(parent, uc_screen_nc_id)
      value  <- c(value, val)
    }
    
    inc_dx_nc_id <- paste0(unc_nc_id, "_incident")
    id     <- c(id, inc_dx_nc_id)
    label  <- c(label, "Incidental dx")
    parent <- c(parent, unc_nc_id)
    value  <- c(value, i_nc)
    
    for (stg in unique(fs_tbl$Stage)) {
      stg_id <- paste0(inc_dx_nc_id, "_", stg)
      val    <- fs_tbl$i_cases_NC[fs_tbl$Stage == stg]
      id     <- c(id, stg_id)
      label  <- c(label, stg)
      parent <- c(parent, inc_dx_nc_id)
      value  <- c(value, val)
    }
  }
  
  stage_levels <- c("I", "II", "III", "IV")
  
  tibble::tibble(id, label, parent, value) %>%
    mutate(
      label = ifelse(label %in% stage_levels,
                     factor(label, levels = stage_levels, ordered = TRUE),
                     label),
      label = as.character(label)
    )
}


make_bar_df_v3 <- function(overall_summary, cancer_type = "All Cancers") {
  
  df <- overall_summary %>%
    dplyr::filter(
      (Arm == "m"  & Variable %in% c("s_cases_FN", "i_cases_FN",
                                     "s_cases_NC", "i_cases_NC",
                                     "s_cases_TP", "i_cases_TP")) |
        (Arm == "uc" & Variable %in% c("s_cases_uc", "i_cases_uc"))
    ) %>%
    mutate(
      XCategory = dplyr::case_when(
        Arm == "m"  & stringr::str_detect(Variable, "^s_") ~ "Screened",
        Arm == "m"  & stringr::str_detect(Variable, "^i_") ~ "Incidental",
        Arm == "uc" & stringr::str_detect(Variable, "^s_") ~ "UC Screened",
        Arm == "uc" & stringr::str_detect(Variable, "^i_") ~ "UC Incidental",
        TRUE                                               ~ NA_character_
      ),
      CaseType = dplyr::case_when(
        Arm == "m" & stringr::str_detect(Variable, "TP") ~ "MCED TP",
        Arm == "m" & stringr::str_detect(Variable, "FN") ~ "MCED FN",
        Arm == "m" & stringr::str_detect(Variable, "NC") ~ "MCED NC",
        Arm == "uc"                                      ~ "Usual Care",
        TRUE                                             ~ NA_character_
      )
    )
  
  if (cancer_type != "All Cancers") {
    df <- df %>% dplyr::filter(Cancer == cancer_type)
  }
  
  df_sum <- df %>%
    group_by(XCategory, CaseType) %>%
    summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(Cancer = cancer_type)
  
  df_sum <- df_sum %>%
    mutate(CaseType = dplyr::recode(CaseType,
                                    "MCED TP" = "MCD TP",
                                    "MCED FN" = "MCD FN",
                                    "MCED NC" = "Not compliant to MCD"))
  
  return(df_sum)
}

case_labels <- c(
  "MCED TP"    = "MCD TP",
  "MCED FN"    = "MCD FN",
  "MCED NC"    = "Not compliant to MCD",
  "Usual Care" = "Usual Care"
)

build_case_plot_v3 <- function(df) {
  
  case_labels <- c(
    "MCED TP"    = "MCD TP",
    "MCED FN"    = "MCD FN",
    "MCED NC"    = "Not compliant to MCD",
    "Usual Care" = "Usual Care"
  )
  
  fill_colors <- c(
    "MCD TP"               = "#FFB9AA",
    "MCD FN"               = "#802815",
    "Not compliant to MCD" = "#fdd49e",
    "Usual Care"           = "#999999"
  )
  
  x_levels <- c("Screened", "Incidental", "UC Screened", "UC Incidental")
  
  df <- df %>%
    mutate(
      CaseTypeRaw = CaseType,
      CaseType    = dplyr::recode(CaseType, !!!case_labels),
      XCategory   = factor(XCategory, levels = x_levels),
      CaseType    = factor(CaseType, levels = names(fill_colors))
    )
  
  gg <- ggplot(df, aes(
    x = XCategory,
    y = Value,
    fill = CaseType,
    text = paste0(
      "XCategory: ", XCategory,
      "<br>Case Type: ", CaseType,
      "<br>Cases: ", scales::comma(round(Value, 0))
    )
  )) +
    geom_col(position = "stack") +
    scale_fill_manual(values = fill_colors) +
    scale_y_continuous(labels = scales::comma) +   # 👈 ADD THIS
    labs(x = NULL, y = "Number of Cases", fill = "Case Type") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )
  
  
  plotly::ggplotly(gg, tooltip = "text")
}



build_overdx_deaths_plot_v2 <- function(overall_summary,
                                        selected_cancer = "All Cancers",
                                        cancer_colors = NULL) {
  
  assign_color <- function(id) {
    matches <- names(cancer_colors)[
      sapply(names(cancer_colors), function(x) str_detect(tolower(id), fixed(tolower(x))))
    ]
    if (length(matches) > 0) {
      return(cancer_colors[[matches[1]]])
    } else {
      return("#cccccc")  # fallback
    }
  }
  
  cancer_color <- assign_color(selected_cancer)
  
  
  overdx_df <- overall_summary %>%
    filter(Arm == "m", Variable == "overdx_Total") %>%
    transmute(Cancer, Value, Metric = "Over-diagnosed")
  
  deaths_df <- overall_summary %>%
    filter(Arm == "m_v_uc", Variable == "Cancer_deaths_averted") %>%
    transmute(Cancer, Value = pmax(Value, 0), Metric = "Deaths averted")
  
  df_long <- bind_rows(overdx_df, deaths_df)
  
  # Subset or aggregate
  if (selected_cancer != "All Cancers") {
    df_long <- filter(df_long, Cancer == selected_cancer)
    
    # Color by metric
    
    
    fill_values <- c("Over-diagnosed" = darken_hex(cancer_color), "Deaths averted" = lighten_hex(cancer_color))
    
    
    gg <- ggplot(df_long, aes(
      x = Metric, y = Value, fill = Metric,
      text = paste0("Metric: ", Metric,
                    "<br>Value: ", scales::comma(Value))
    )) +
      geom_col(width = 0.6) +
      scale_fill_manual(values = fill_values)
    
  } else {
    df_long <- df_long %>%
      group_by(Metric) %>%
      summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop")
    
    gg <- ggplot(df_long, aes(
      x = Metric, y = Value, fill = Metric,
      text = paste0("Metric: ", Metric,
                    "<br>Value: ", scales::comma(Value))
    )) +
      geom_col(width = 0.6) +
      scale_fill_manual(values = c("Over-diagnosed" = "#550F00",
                                   "Deaths averted" = "#FFB9AA"))
  }
  
  gg <- gg +
    scale_y_continuous(labels = scales::comma) +   # 👈 ADD
    labs(
      x = NULL,
      y = "Number of Cases"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
  
  plotly::ggplotly(gg, tooltip = "text")
}

build_pathway_df <- function(result_list,
                             pathway_type = c("MCED FP", "MCED TN", "Non-Compliant", "ALL"),
                             cancer = "Overall",
                             screening_cancers = c("lung and bronchus", "breast", "colon and rectum", 
                                                   "cervical", "prostate"),
                             for_sunburst = FALSE) {
  pathway_type <- match.arg(pathway_type)
  
  # Always pull shared base steps from Overall
  overall_df <- result_list[["Overall"]][["mced_TN_FP"]]
  base_totals <- colSums(overall_df, na.rm = TRUE)
  
  # Use cancer-specific data only if cancer != "Overall"
  if (cancer == "Overall") {
    cancer_df <- NULL
    cancer_totals <- rep(0, 0)
    is_screening <- TRUE
  } else {
    cancer_df <- result_list[["MCED"]][[cancer]][["mced_TN_FP"]]
    cancer_totals <- colSums(cancer_df, na.rm = TRUE)
    is_screening <- cancer %in% screening_cancers
  }
  
  if (pathway_type == "ALL") {
    df_fp <- build_pathway_df(result_list, "MCED FP", cancer, screening_cancers, for_sunburst)
    df_tn <- build_pathway_df(result_list, "MCED TN", cancer, screening_cancers, for_sunburst)
    df_nc <- build_pathway_df(result_list, "Non-Compliant", cancer, screening_cancers, for_sunburst)
    
    if (!for_sunburst) {
      # Adjust TN and NC to avoid stacking base steps (only for bar plot)
      df_tn$Count[df_tn$Step %in% c("Healthy", "MCED Compliant")] <- 0
      df_nc$Count[df_nc$Step %in% c("Healthy", "MCED Result")]     <- 0
    } else {
      # Don't recode compliance labels either
    }
    
    df_all <- dplyr::bind_rows(
      dplyr::mutate(df_fp, SubPathway = "MCED FP"),
      dplyr::mutate(df_tn, SubPathway = "MCED TN"),
      dplyr::mutate(df_nc, SubPathway = "Non-Compliant")
    )
    return(df_all)
  }
  
  if (pathway_type == "MCED FP") {
    steps <- c(
      "Healthy", "MCED Compliant", "MCED Result",
      "MCED Follow-up (FPx1)", "UC Compliant",
      "UC FP (FPx2)", "UC Invasive (FPx3)"
    )
    
    uc_compliance <- if (!is_screening) 0 else if (!is.null(cancer_df)) cancer_totals["FP_fu"] else base_totals["Total_FP_fu"]
    uc_fp <- if (!is_screening) {
      if (!is.null(cancer_df)) cancer_totals["FP_fu"] else base_totals["Total_FP_fu"]
    } else {
      if (!is.null(cancer_df)) cancer_totals["FP_FP"] else base_totals["FP_FP"]
    }
    
    counts <- c(
      base_totals["H"],
      base_totals["compliant_pop"],
      if (!is.null(cancer_df)) cancer_totals["FP"] else base_totals["fp_counts"],
      if (!is.null(cancer_df)) cancer_totals["FP_fu"] else base_totals["Total_FP_fu"],
      uc_compliance,
      uc_fp,
      if (!is.null(cancer_df)) cancer_totals["FP_FP_i"] + cancer_totals["FP_i"] else base_totals["FP_FP_i"] + base_totals["FP_i"]
    )
    
  } else if (pathway_type == "MCED TN") {
    steps <- c(
      "Healthy", "MCED Compliant", "MCED Result",
      "MCED Follow-up (FPx1)", "UC Compliant",
      "UC FP (FPx2)", "UC Invasive (FPx3)"
    )
    
    counts <- c(
      base_totals["H"],
      base_totals["compliant_pop"],
      base_totals["True_TN"],
      0,
      if (!is.null(cancer_df)) cancer_totals["TN_c"] else base_totals["TN_c"],
      if (!is.null(cancer_df)) cancer_totals["TN_c_FP"] else base_totals["TN_c_FP"],
      if (!is.null(cancer_df)) cancer_totals["TN_c_FP_i"] else base_totals["TN_c_FP_i"]
    )
    
  } else if (pathway_type == "Non-Compliant") {
    steps <- c(
      "Healthy", "MCED Non-Compliant", "MCED Result",
      "MCED Follow-up (FPx1)", "UC Compliant",
      "UC FP (FPx2)", "UC Invasive (FPx3)"
    )
    
    counts <- c(
      base_totals["H"],
      base_totals["non_compliant_pop"],
      0,
      0,
      if (!is.null(cancer_df)) cancer_totals["nc_compliant"] else base_totals["nc_compliant"],
      if (!is.null(cancer_df)) cancer_totals["nc_FP"] else base_totals["nc_FP"],
      if (!is.null(cancer_df)) cancer_totals["nc_Invasive"] else base_totals["nc_Invasive"]
    )
  }
  
  data.frame(
    Cancer = cancer,
    Pathway = pathway_type,
    Step = factor(steps, levels = steps),
    Count = counts
  )
}


build_cost_breakdown_df <- function(final_icer_result, overall_summary, cancer = "All Cancers") {
  get_val <- function(df, varname) {
    val <- df$Value[df$Variable == varname]
    if (length(val) == 0) return(0)
    return(val[1])
  }
  
  get_sum_across_cancers <- function(df, varname, arm, exclude_all = TRUE) {
    filtered <- df %>% filter(Variable == varname, Arm == arm)
    if (exclude_all) {
      filtered <- filtered %>% filter(Cancer != "All Cancers")
    }
    sum(filtered$Value, na.rm = TRUE)
  }
  
  # Fixed test costs (same for all cancer types)
  test_cost_cancer  <- get_val(final_icer_result, "cancer_test_cost")
  test_cost_healthy <- get_val(final_icer_result, "tn_fp_mced_test_cost")
  total_test_cost   <- test_cost_cancer + test_cost_healthy
  
  bucket_names <- c("Test Cost", "Imaging Cost", "Diagnostic Cost", "Invasive Cost", "Post Diagnosis Costs")
  
  if (cancer == "All Cancers") {
    df <- tibble::tibble(
      Bucket = bucket_names,
      Cancer = c(
        test_cost_cancer,
        get_sum_across_cancers(overall_summary, "c_Imaging_cost", "m"),
        get_sum_across_cancers(overall_summary, "c_Diagnostic_cost", "m"),
        get_sum_across_cancers(overall_summary, "c_Invasive_cost", "m"),
        get_sum_across_cancers(overall_summary, "I_post_cost", "m") +
          get_sum_across_cancers(overall_summary, "II_post_cost", "m") +
          get_sum_across_cancers(overall_summary, "III_post_cost", "m") +
          get_sum_across_cancers(overall_summary, "IV_post_cost", "m")
      ),
      Healthy = c(
        test_cost_healthy,
        get_sum_across_cancers(overall_summary, "h_Imaging_cost", "m"),
        get_sum_across_cancers(overall_summary, "h_Diagnostic_cost", "m"),
        get_sum_across_cancers(overall_summary, "h_Invasive_cost", "m"),
        0
      )
    )
  } else {
    os <- overall_summary %>% filter(Cancer == cancer, Arm == "m")
    
    df <- tibble::tibble(
      Bucket = bucket_names,
      Cancer = c(
        test_cost_cancer,
        get_val(os, "c_Imaging_cost"),
        get_val(os, "c_Diagnostic_cost"),
        get_val(os, "c_Invasive_cost"),
        sum(os$Value[os$Variable %in% c("I_post_cost", "II_post_cost", "III_post_cost", "IV_post_cost")])
      ),
      Healthy = c(
        test_cost_healthy,
        get_val(os, "h_Imaging_cost"),
        get_val(os, "h_Diagnostic_cost"),
        get_val(os, "h_Invasive_cost"),
        0
      )
    )
  }
  
  # Do NOT include test costs in "Overall" column
  df <- df %>%
    mutate(Overall = dplyr::if_else(Bucket == "Test Cost", total_test_cost, Cancer + Healthy))
  
  df_long <- df %>%
    tidyr::pivot_longer(cols = c("Cancer", "Healthy", "Overall"), names_to = "Group", values_to = "Cost")
  
  return(df_long)
}
