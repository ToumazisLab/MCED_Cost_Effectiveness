count_tp_incidentals <- function(matrix, screening_data, screening_cancers, 
                                 final_result_list, arm = "usual") {
  # Initialize matrices
  s_ne_i_dx_matrix <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))  
  s_e_i_dx_matrix <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))   
  s_e_s_dx_matrix <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))   
  
  colnames(s_ne_i_dx_matrix) <- colnames(matrix)
  colnames(s_e_i_dx_matrix) <- colnames(matrix)
  colnames(s_e_s_dx_matrix) <- colnames(matrix)
  
  # Track which cancers we have processed to add totals at the end
  processed_cancers <- character()
  
  for (col_name in colnames(matrix)) {
    # Extract stage and cancer_type
    parts <- strsplit(col_name, "_")[[1]]
    stage <- parts[1]
    cancer_type <- tolower(paste(parts[-1], collapse = "_")) 
    
    # Initialize the cancer data frame if needed
    if (!("Cancer_cases_screening" %in% names(final_result_list[[arm]][[cancer_type]]))) {
      final_result_list[[arm]][[cancer_type]]$Cancer_cases_screening <- create_cancer_df(nrow(matrix), arm = arm)
    }
    
    # Classify each cycle
    for (cycle in 1:nrow(matrix)) {
      total_cases <- matrix[cycle, col_name]
      result <- classify_cases(total_cases, cancer_type, cycle, screening_data, screening_cancers)
      
      # Store results in matrices
      s_ne_i_dx_matrix[cycle, col_name] <- result$s_ne_i
      s_e_i_dx_matrix[cycle, col_name] <- result$s_e_i
      s_e_s_dx_matrix[cycle, col_name] <- result$s_e_s
      
      # Update final_result_list data frame
      df <- final_result_list[[arm]][[cancer_type]]$Cancer_cases_screening
      df[cycle, stage] <- total_cases
      df[cycle, paste0("s_e_s_", stage)] <- result$s_e_s
      df[cycle, paste0("s_e_i_", stage)] <- result$s_e_i
      df[cycle, paste0("s_ne_i_", stage)] <- result$s_ne_i
      
      final_result_list[[arm]][[cancer_type]]$Cancer_cases_screening <- df
    }
    
    # Mark this cancer as processed for later summary
    processed_cancers <- unique(c(processed_cancers, cancer_type))
  }
  
  # After processing all columns, add totals by cancer
  for (cancer_type in processed_cancers) {
    df <- final_result_list[[arm]][[cancer_type]]$Cancer_cases_screening
    
    # Calculate totals by stage and across stages
    df$s_e_s_total <- df$s_e_s_I + df$s_e_s_II + df$s_e_s_III + df$s_e_s_IV
    df$s_e_i_total <- df$s_e_i_I + df$s_e_i_II + df$s_e_i_III + df$s_e_i_IV
    df$s_ne_i_total <- df$s_ne_i_I + df$s_ne_i_II + df$s_ne_i_III + df$s_ne_i_IV
    
    # s_dx = s_e_s_total, i_dx = s_e_i_total + s_ne_i_total
    df$s_dx_I <- df$s_e_s_I
    df$s_dx_II <- df$s_e_s_II
    df$s_dx_III <- df$s_e_s_III
    df$s_dx_IV <- df$s_e_s_IV
    df$i_dx_I <- df$s_e_i_I + df$s_ne_i_I
    df$i_dx_II <- df$s_e_i_II + df$s_ne_i_II
    df$i_dx_III <- df$s_e_i_III + df$s_ne_i_III
    df$i_dx_IV <- df$s_e_i_IV + df$s_ne_i_IV
    
    df$s_dx_total <- df$s_e_s_total
    df$i_dx_total <- df$s_e_i_total + df$s_ne_i_total
    df$Total_cases <- df$I + df$II + df$III + df$IV
    final_result_list[[arm]][[cancer_type]]$Cancer_cases_screening <- df
  }
  
  # Summarize results across all cancers
  sum_stage <- function(mat, st) rowSums(mat[, grepl(paste0("^", st, "_"), colnames(mat)), drop=FALSE])
  
  overall_summary <- data.frame(
    Cycle = 1:nrow(matrix),
    I = sum_stage(matrix, "I"),
    II = sum_stage(matrix, "II"),
    III = sum_stage(matrix, "III"),
    IV = sum_stage(matrix, "IV"),
    s_e_s_I = sum_stage(s_e_s_dx_matrix, "I"),
    s_e_s_II = sum_stage(s_e_s_dx_matrix, "II"),
    s_e_s_III = sum_stage(s_e_s_dx_matrix, "III"),
    s_e_s_IV = sum_stage(s_e_s_dx_matrix, "IV"),
    s_e_i_I = sum_stage(s_e_i_dx_matrix, "I"),
    s_e_i_II = sum_stage(s_e_i_dx_matrix, "II"),
    s_e_i_III = sum_stage(s_e_i_dx_matrix, "III"),
    s_e_i_IV = sum_stage(s_e_i_dx_matrix, "IV"),
    s_ne_i_I = sum_stage(s_ne_i_dx_matrix, "I"),
    s_ne_i_II = sum_stage(s_ne_i_dx_matrix, "II"),
    s_ne_i_III = sum_stage(s_ne_i_dx_matrix, "III"),
    s_ne_i_IV = sum_stage(s_ne_i_dx_matrix, "IV")
  )
  
  # Calculate total case summaries across all stages
  overall_summary$Total_cases <- rowSums(overall_summary[, c("I", "II", "III", "IV")])
  overall_summary$s_e_s_total <- rowSums(overall_summary[, c("s_e_s_I", "s_e_s_II", "s_e_s_III", "s_e_s_IV")])
  overall_summary$s_e_i_total <- rowSums(overall_summary[, c("s_e_i_I", "s_e_i_II", "s_e_i_III", "s_e_i_IV")])
  overall_summary$s_ne_i_total <- rowSums(overall_summary[, c("s_ne_i_I", "s_ne_i_II", "s_ne_i_III", "s_ne_i_IV")])
  
  # s_dx and i_dx in the overall
  # s_dx by stage = s_e_s_ by stage
  overall_summary$s_dx_I <- overall_summary$s_e_s_I
  overall_summary$s_dx_II <- overall_summary$s_e_s_II
  overall_summary$s_dx_III <- overall_summary$s_e_s_III
  overall_summary$s_dx_IV <- overall_summary$s_e_s_IV
  
  # i_dx by stage = s_e_i_ + s_ne_i_ by stage
  overall_summary$i_dx_I <- overall_summary$s_e_i_I + overall_summary$s_ne_i_I
  overall_summary$i_dx_II <- overall_summary$s_e_i_II + overall_summary$s_ne_i_II
  overall_summary$i_dx_III <- overall_summary$s_e_i_III + overall_summary$s_ne_i_III
  overall_summary$i_dx_IV <- overall_summary$s_e_i_IV + overall_summary$s_ne_i_IV
  
  overall_summary$s_dx_total <- overall_summary$s_e_s_total
  overall_summary$i_dx_total <- overall_summary$s_e_i_total + overall_summary$s_ne_i_total
  
  final_result_list[["Overall"]][["Cancer_cases_screening_summary"]] <- overall_summary
  
  return(list(
    final_result_list = final_result_list,
    s_ne_i_dx_matrix = s_ne_i_dx_matrix,
    s_e_i_dx_matrix = s_e_i_dx_matrix,
    s_e_s_dx_matrix = s_e_s_dx_matrix
  ))
}

# Helpers
sum_by_cancer <- function(mat, all_cols) {
  cancer_names <- sapply(all_cols, function(x)
    tolower(paste(strsplit(x, "_")[[1]][-1], collapse = "_")))
  tapply(colSums(mat), cancer_names, sum)
}

sum_by_cancer_stage <- function(mat) {
  parts      <- strsplit(colnames(mat), "_")
  stage_vec  <- sapply(parts, `[`, 1)                                   # I, II, …
  cancer_vec <- sapply(parts, function(x) tolower(paste(x[-1], collapse = "_")))
  key        <- paste(stage_vec, cancer_vec, sep = "|")
  totals     <- tapply(colSums(mat), key, sum)
  
  ## return a tidy data frame
  df <- data.frame(
    Stage  = sub("\\|.*", "", names(totals)),
    Cancer = sub(".*\\|",   "", names(totals)),
    Total  = as.numeric(totals),
    row.names = NULL
  )
  df
}


mced_cancer_categories <- function(matrix_TP,
                                   matrix_FN,
                                   matrix_NC,
                                   screening_data,
                                   screening_cancers,
                                   final_result_list,
                                   arm = "MCED",
                                   reduced_comply) {
  
  # --------------------
  # 1. Gather columns and create empty bucket matrices
  # --------------------
  all_cols <- unique(c(colnames(matrix_TP),
                       colnames(matrix_NC),
                       colnames(matrix_FN))) 
  ncycles <- nrow(matrix_TP)  # Assumes matrix_NC has the same number of rows
  
  s_e_m_matrix  <- matrix(0, nrow = ncycles, ncol = length(all_cols), dimnames = list(NULL, all_cols))
  s_e_s_matrix  <- matrix(0, nrow = ncycles, ncol = length(all_cols), dimnames = list(NULL, all_cols))
  s_e_i_matrix  <- matrix(0, nrow = ncycles, ncol = length(all_cols), dimnames = list(NULL, all_cols))
  s_ne_i_matrix <- matrix(0, nrow = ncycles, ncol = length(all_cols), dimnames = list(NULL, all_cols))
  
  s_e_m_matrix_TP <- matrix(0, nrow = ncycles, ncol = length(all_cols), dimnames = list(NULL, all_cols))
  
  s_e_s_matrix_FN  <- matrix(0, nrow = ncycles, ncol = length(all_cols), dimnames = list(NULL, all_cols))
  s_e_i_matrix_FN  <- matrix(0, nrow = ncycles, ncol = length(all_cols), dimnames = list(NULL, all_cols))
  s_ne_i_matrix_FN <- matrix(0, nrow = ncycles, ncol = length(all_cols), dimnames = list(NULL, all_cols))
  
  s_e_s_matrix_NC  <- matrix(0, nrow = ncycles, ncol = length(all_cols), dimnames = list(NULL, all_cols))
  s_e_i_matrix_NC  <- matrix(0, nrow = ncycles, ncol = length(all_cols), dimnames = list(NULL, all_cols))
  s_ne_i_matrix_NC <- matrix(0, nrow = ncycles, ncol = length(all_cols), dimnames = list(NULL, all_cols))
  
  # sum_by_cancer <- function(mat, all_cols) {
  #   cancer_names <- sapply(all_cols, function(x) tolower(paste(strsplit(x, "_")[[1]][-1], collapse = "_")))
  #   summed <- tapply(colSums(mat), cancer_names, sum)
  #   return(summed)
  # }
  
  # --------------------
  # 2. Fill TP buckets
  # --------------------
  if (!is.null(matrix_TP)) {
    
    s_e_m_matrix[, colnames(matrix_TP)] <- matrix_TP
    
    s_e_m_matrix_TP[, colnames(matrix_TP)] <- matrix_TP
    
    # existing code …
    s_e_m_totals         <- sum_by_cancer(s_e_m_matrix, colnames(matrix_TP))
    
    ## --- stage-specific version --------------------------------------------¶ NEW
    s_e_m_totals_stage   <- sum_by_cancer_stage(s_e_m_matrix)
    m_summary_df_stage   <- dplyr::rename(s_e_m_totals_stage,
                                          s_e_m_TP = Total)
    
    m_summary_df <- data.frame(
      Cancer  = names(s_e_m_totals),
      s_e_m_TP     = as.numeric(s_e_m_totals),
      row.names    = NULL
    )
    
    overall_M_totals <- data.frame(
      s_e_m_TP_total = sum(s_e_m_totals)
    )
  } else {
    stop("No TP Cases")
  }
  
  fn_summary_df_stage <- data.frame(Cancer = character(),
                                    Stage  = character())
  # --------------------
  # 3a. Handle matrix_FN via classify_cases() logic
  # --------------------
  if (!is.null(matrix_FN)) {
    for (col_name in colnames(matrix_FN)) {
      parts <- strsplit(col_name, "_")[[1]]
      cancer_type <- tolower(paste(parts[-1], collapse = "_"))
      
      for (cycle in seq_len(nrow(matrix_FN))) {
        total_cases <- matrix_FN[cycle, col_name]
        result <- classify_cases(total_cases, cancer_type, cycle, screening_data, screening_cancers, case_type = "FN", reduced_comply = reduced_comply)
        
        s_e_s_matrix[cycle, col_name]  <- s_e_s_matrix[cycle, col_name]  + result$s_e_s
        s_e_i_matrix[cycle, col_name]  <- s_e_i_matrix[cycle, col_name]  + result$s_e_i
        s_ne_i_matrix[cycle, col_name] <- s_ne_i_matrix[cycle, col_name] + result$s_ne_i
        
        s_e_s_matrix_FN[cycle, col_name]  <- s_e_s_matrix_FN[cycle, col_name]  + result$s_e_s
        s_e_i_matrix_FN[cycle, col_name]  <- s_e_i_matrix_FN[cycle, col_name]  + result$s_e_i
        s_ne_i_matrix_FN[cycle, col_name] <- s_ne_i_matrix_FN[cycle, col_name] + result$s_ne_i
      }
    }
    
    s_e_s_FN_totals  <- sum_by_cancer(s_e_s_matrix_FN, colnames(matrix_FN))
    s_e_i_FN_totals  <- sum_by_cancer(s_e_i_matrix_FN, colnames(matrix_FN))
    s_ne_i_FN_totals <- sum_by_cancer(s_ne_i_matrix_FN, colnames(matrix_FN))
    
    fn_summary_df <- data.frame(
      Cancer   = names(s_e_s_FN_totals),
      s_e_s_FN      = as.numeric(s_e_s_FN_totals),
      s_e_i_FN      = as.numeric(s_e_i_FN_totals),
      s_ne_i_FN     = as.numeric(s_ne_i_FN_totals),
      row.names     = NULL
    )
    
    # after the loop & cancer-totals you already build …
    s_e_s_FN_totals_stage  <- sum_by_cancer_stage(s_e_s_matrix_FN)
    s_e_i_FN_totals_stage  <- sum_by_cancer_stage(s_e_i_matrix_FN)
    s_ne_i_FN_totals_stage <- sum_by_cancer_stage(s_ne_i_matrix_FN)
    
    fn_summary_df_stage <- reduce(
      list(
        dplyr::rename(s_e_s_FN_totals_stage,  s_e_s_FN  = Total),
        dplyr::rename(s_e_i_FN_totals_stage,  s_e_i_FN  = Total),
        dplyr::rename(s_ne_i_FN_totals_stage, s_ne_i_FN = Total)
      ),
      dplyr::full_join, by = c("Cancer", "Stage")
    )
    fn_summary_df_stage$s_cases_FN <- fn_summary_df_stage$s_e_s_FN
    fn_summary_df_stage$i_cases_FN <- fn_summary_df_stage$s_e_i_FN +
      fn_summary_df_stage$s_ne_i_FN
    
    
    overall_FN_totals <- data.frame(
      s_e_s_FN_total  = sum(s_e_s_FN_totals),
      s_e_i_FN_total  = sum(s_e_i_FN_totals),
      s_ne_i_FN_total = sum(s_ne_i_FN_totals)
    )
  }
  
  nc_summary_df_stage <- data.frame(Cancer = character(),
                                    Stage  = character()) # ← 3
  
  # --------------------
  # 3b. Handle matrix_NC via classify_cases() logic
  # --------------------
  if (!is.null(matrix_NC)) {
    for (col_name in colnames(matrix_NC)) {
      parts <- strsplit(col_name, "_")[[1]]
      cancer_type <- tolower(paste(parts[-1], collapse = "_"))
      
      for (cycle in seq_len(nrow(matrix_NC))) {
        total_cases <- matrix_NC[cycle, col_name]
        result <- classify_cases(total_cases, cancer_type, cycle, screening_data, screening_cancers, case_type = "NA", reduced_comply = NULL)
        
        s_e_s_matrix[cycle, col_name]  <- s_e_s_matrix[cycle, col_name]  + result$s_e_s
        s_e_i_matrix[cycle, col_name]  <- s_e_i_matrix[cycle, col_name]  + result$s_e_i
        s_ne_i_matrix[cycle, col_name] <- s_ne_i_matrix[cycle, col_name] + result$s_ne_i
        
        s_e_s_matrix_NC[cycle, col_name]  <- s_e_s_matrix_NC[cycle, col_name]  + result$s_e_s
        s_e_i_matrix_NC[cycle, col_name]  <- s_e_i_matrix_NC[cycle, col_name]  + result$s_e_i
        s_ne_i_matrix_NC[cycle, col_name] <- s_ne_i_matrix_NC[cycle, col_name] + result$s_ne_i
      }
    }
    
    s_e_s_NC_totals  <- sum_by_cancer(s_e_s_matrix_NC, colnames(matrix_NC))
    s_e_i_NC_totals  <- sum_by_cancer(s_e_i_matrix_NC, colnames(matrix_NC))
    s_ne_i_NC_totals <- sum_by_cancer(s_ne_i_matrix_NC, colnames(matrix_NC))
    
    nc_summary_df <- data.frame(
      Cancer   = names(s_e_s_NC_totals),
      s_e_s_NC      = as.numeric(s_e_s_NC_totals),
      s_e_i_NC      = as.numeric(s_e_i_NC_totals),
      s_ne_i_NC     = as.numeric(s_ne_i_NC_totals),
      row.names     = NULL
    )
    
    s_e_s_NC_totals_stage  <- sum_by_cancer_stage(s_e_s_matrix_NC)
    s_e_i_NC_totals_stage  <- sum_by_cancer_stage(s_e_i_matrix_NC)
    s_ne_i_NC_totals_stage <- sum_by_cancer_stage(s_ne_i_matrix_NC)
    
    nc_summary_df_stage <- reduce(
      list(
        dplyr::rename(s_e_s_NC_totals_stage,  s_e_s_NC  = Total),
        dplyr::rename(s_e_i_NC_totals_stage,  s_e_i_NC  = Total),
        dplyr::rename(s_ne_i_NC_totals_stage, s_ne_i_NC = Total)
      ),
      dplyr::full_join, by = c("Cancer", "Stage")
    )
    nc_summary_df_stage$s_cases_NC <- nc_summary_df_stage$s_e_s_NC
    nc_summary_df_stage$i_cases_NC <- nc_summary_df_stage$s_e_i_NC +
      nc_summary_df_stage$s_ne_i_NC
    
    
    overall_NC_totals <- data.frame(
      s_e_s_NC_total  = sum(s_e_s_NC_totals),
      s_e_i_NC_total  = sum(s_e_i_NC_totals),
      s_ne_i_NC_total = sum(s_ne_i_NC_totals)
    )
  }
  
  # --------------------
  # Summarize FN and NC totals by cancer
  # --------------------
  fn_summary_df$s_cases_FN <- fn_summary_df$s_e_s_FN
  fn_summary_df$i_cases_FN <- fn_summary_df$s_e_i_FN + fn_summary_df$s_ne_i_FN
  
  nc_summary_df$s_cases_NC <- nc_summary_df$s_e_s_NC
  nc_summary_df$i_cases_NC <- nc_summary_df$s_e_i_NC + nc_summary_df$s_ne_i_NC
  
  m_summary_df$s_cases_TP <- m_summary_df$s_e_m_TP
  m_summary_df$i_cases_TP <- 0
  
  m_summary_df_stage$s_cases_TP <- m_summary_df_stage$s_e_m_TP
  m_summary_df_stage$i_cases_TP <- 0
  
  final_screen_df_stage <- fn_summary_df_stage %>%
    dplyr::full_join(nc_summary_df_stage, by = c("Cancer", "Stage")) %>%
    dplyr::full_join(m_summary_df_stage,  by = c("Cancer", "Stage")) %>%
    dplyr::arrange(Cancer, Stage)
  
  final_screen_df_stage <- final_screen_df_stage %>%        # ← 4
    dplyr::mutate(Stage = factor(Stage, levels = c("I","II","III","IV"))) %>%
    dplyr::arrange(Cancer, Stage)
  
  # Join all by-cancer summaries
  temp_df <- fn_summary_df %>%
    dplyr::full_join(nc_summary_df, by = "Cancer") %>%
    dplyr::full_join(m_summary_df, by = "Cancer")
  
  # Build overall row
  overall_screen_df <- data.frame(
    Cancer    = "Overall",
    s_e_s_FN       = overall_FN_totals$s_e_s_FN_total,
    s_e_i_FN       = overall_FN_totals$s_e_i_FN_total,
    s_ne_i_FN      = overall_FN_totals$s_ne_i_FN_total,
    s_cases_FN     = overall_FN_totals$s_e_s_FN_total,
    i_cases_FN     = overall_FN_totals$s_e_i_FN_total + overall_FN_totals$s_ne_i_FN_total,
    s_e_s_NC       = overall_NC_totals$s_e_s_NC_total,
    s_e_i_NC       = overall_NC_totals$s_e_i_NC_total,
    s_ne_i_NC      = overall_NC_totals$s_ne_i_NC_total,
    s_cases_NC     = overall_NC_totals$s_e_s_NC_total,
    i_cases_NC     = overall_NC_totals$s_e_i_NC_total + overall_NC_totals$s_ne_i_NC_total,
    s_e_m_TP       = overall_M_totals$s_e_m_TP_total,
    s_cases_TP     = overall_M_totals$s_e_m_TP_total,
    i_cases_TP     = 0
  )
  
  # Final combined dataframe
  final_screen_df <- dplyr::bind_rows(temp_df, overall_screen_df)
  
  # --------------------
  # 4. Update final_result_list for each cancer type & stage
  # --------------------
  unique_cancers <- unique(
    tolower(sapply(all_cols, function(x) paste(strsplit(x, "_")[[1]][-1], collapse = "_")))
  )
  
  stages <- c("I", "II", "III", "IV")
  
  # A helper to ensure we don't overwrite existing columns but DO add missing ones
  add_missing_columns <- function(df, template_df) {
    # If a column from template_df is missing in df, add it as 0
    for (col_name in colnames(template_df)) {
      if (!col_name %in% colnames(df)) {
        df[[col_name]] <- 0
      }
    }
    # Return df with all the needed columns added, plus any original columns
    return(df)
  }
  
  for (col_name in all_cols) {
    parts <- strsplit(col_name, "_")[[1]]
    stage <- parts[1]
    cancer_type <- tolower(paste(parts[-1], collapse = "_"))
    
    # 4a. Get or create the data frame for this arm/cancer_type
    df <- final_result_list[[arm]][[cancer_type]]$MCED_Cancer_cases
    if (is.null(df)) {
      # If there's no existing df at all, start with a fresh template
      df <- create_cancer_df(ncycles)
    } else {
      # Otherwise, add any missing columns from the template (don't overwrite)
      df <- add_missing_columns(df, create_cancer_df(ncycles))
      # Also ensure row counts match; if nrow(df) < ncycles, you might need to expand
      if (nrow(df) < ncycles) {
        # Example: if you need to enlarge df, carefully rbind more rows
        extra_rows <- ncycles - nrow(df)
        df2 <- create_cancer_df(extra_rows)  # brand new rows
        # Copy your existing columns from df into the new structure:
        df2 <- add_missing_columns(df2, df)  # ensure same columns
        df  <- rbind(df, df2)
      }
    }
    
    # 4b. Accumulate the relevant stage/bucket counts
    for (cycle in seq_len(ncycles)) {
      total_cases <- s_e_m_matrix[cycle, col_name] +
        s_e_s_matrix[cycle, col_name] +
        s_e_i_matrix[cycle, col_name] +
        s_ne_i_matrix[cycle, col_name]
      
      df[cycle, paste0(stage, "_calculated")] <- total_cases # df[cycle, stage] + 
      df[cycle, paste0("s_e_m_", stage)]   <- df[cycle, paste0("s_e_m_", stage)]   + s_e_m_matrix[cycle, col_name]
      df[cycle, paste0("s_e_s_", stage)]   <- df[cycle, paste0("s_e_s_", stage)]   + s_e_s_matrix[cycle, col_name]
      df[cycle, paste0("s_e_i_", stage)]   <- df[cycle, paste0("s_e_i_", stage)]   + s_e_i_matrix[cycle, col_name]
      df[cycle, paste0("s_ne_i_", stage)]  <- df[cycle, paste0("s_ne_i_", stage)]  + s_ne_i_matrix[cycle, col_name]
    }
    
    final_result_list[[arm]][[cancer_type]]$MCED_Cancer_cases <- df
  }
  
  # --------------------
  # 5. Compute totals/derived variables within each cancer df
  # --------------------
  for (cancer_type in unique_cancers) {
    df <- final_result_list[[arm]][[cancer_type]]$MCED_Cancer_cases
    if (!is.null(df)) {
      # Again, add any missing columns in case some steps changed
      df <- add_missing_columns(df, create_cancer_df(ncycles))
      
      # Totals by detection category
      df$s_e_m_total  <- df$s_e_m_I + df$s_e_m_II + df$s_e_m_III + df$s_e_m_IV
      df$s_e_s_total  <- df$s_e_s_I + df$s_e_s_II + df$s_e_s_III + df$s_e_s_IV
      df$s_e_i_total  <- df$s_e_i_I + df$s_e_i_II + df$s_e_i_III + df$s_e_i_IV
      df$s_ne_i_total <- df$s_ne_i_I + df$s_ne_i_II + df$s_ne_i_III + df$s_ne_i_IV
      
      # Overall total
      df$Total_cases <- df$I_calculated + df$II_calculated + df$III_calculated + df$IV_calculated
      
      # Derived columns
      for (stg in stages) {
        df[[paste0("m_dx_", stg)]] <- df[[paste0("s_e_m_", stg)]]
        df[[paste0("s_dx_", stg)]] <- df[[paste0("s_e_s_", stg)]]
        df[[paste0("i_dx_", stg)]] <- df[[paste0("s_e_i_", stg)]] + df[[paste0("s_ne_i_", stg)]]
      }
      
      df$m_dx_total <- df$s_e_m_total
      df$s_dx_total <- df$s_e_s_total
      df$i_dx_total <- df$s_e_i_total + df$s_ne_i_total
      
      final_result_list[[arm]][[cancer_type]]$MCED_Cancer_cases <- df
    }
  }
  
  # --------------------
  # 6. Build / Update the overall DF with the same columns
  # --------------------
  overall_summary_name <- if (arm == "usual") {
    "Cancer_cases_screening_summary"
  } else {
    "MCED_Cancer_cases_screening_summary"
  }
  
  overall_df <- final_result_list[["Overall"]][[overall_summary_name]]
  if (is.null(overall_df)) {
    overall_df <- create_cancer_df(ncycles)
  } else {
    # Add any missing columns from template
    overall_df <- add_missing_columns(overall_df, create_cancer_df(ncycles))
    # If we need more rows, handle that similarly to above
    if (nrow(overall_df) < ncycles) {
      extra_rows <- ncycles - nrow(overall_df)
      df2 <- create_cancer_df(extra_rows)
      # Merge columns
      df2 <- add_missing_columns(df2, overall_df)
      overall_df <- rbind(overall_df, df2)
    }
  }
  
  # We'll sum up all the cancers for each column. Because we need to do row sums
  # of the final "MCED_Cancer_cases" DFs, let's gather them:
  
  # First, zero out the columns in overall_df so we can 'add onto' them
  # (We don’t want to double-count if columns already have values)
  for (coln in setdiff(colnames(overall_df), "Cycle")) {
    overall_df[[coln]] <- 0
  }
  
  # For each cancer type, add its per-cycle columns
  # We assume same row structure, i.e. "Cycle" 1..ncycles
  for (cancer_type in unique_cancers) {
    df <- final_result_list[[arm]][[cancer_type]]$MCED_Cancer_cases
    if (!is.null(df)) {
      # Ensure df has the same shape & columns
      df <- add_missing_columns(df, overall_df)
      
      # Now add each column's values to overall_df (except "Cycle")
      # i.e. overall_df[col] += df[col]
      for (coln in setdiff(colnames(overall_df), "Cycle")) {
        overall_df[[coln]] <- overall_df[[coln]] + df[[coln]]
      }
    }
  }
  
  # Save updated overall_df
  final_result_list[["Overall"]][[overall_summary_name]] <- overall_df
  
  # --------------------
  # 7. Return everything
  # --------------------
  return(list(
    final_result_list      = final_result_list,
    s_e_m_matrix           = s_e_m_matrix,
    s_e_s_matrix           = s_e_s_matrix,
    s_e_i_matrix           = s_e_i_matrix,
    s_ne_i_matrix          = s_ne_i_matrix,
    final_screen_df        = final_screen_df,
    final_screen_df_stage  = final_screen_df_stage         
  ))
}


