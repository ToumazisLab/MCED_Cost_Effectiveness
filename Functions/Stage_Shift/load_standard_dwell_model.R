load_seer_dwell_model <- function(name_mapping, cancer_list) {
  # Define stages
  xStage <- c("I", "II", "III", "IV")
  
  # Construct file paths dynamically
  dwell_model_group_file <- file.path("Data/raw/20200728_dwell_time_groups.tsv")
  dwell_model_timing_file <- file.path("Data/raw/20200728_dwell_group_timing.tsv")
  
  # Read the input files
  dwell_model_group_df <- read_tsv(dwell_model_group_file)
  dwell_model_timing_df <- read_tsv(dwell_model_timing_file)
  
  # Process and join the datasets
  dwell_model_all_df <- dwell_model_timing_named %>%
    rename(dwell_group = group) %>%
    full_join(dwell_model_timing_df %>% rename(dwell_group = group), relationship = "many-to-many") %>%
    mutate(number_stage = match(Stage, xStage))
  
  # Apply mapping and filter based on cancer list
  dwell_model_all_df <- apply_mapping_and_filter(dwell_model_all_df, name_mapping, cancer_list)
  
  # Return the processed data frame
  return(dwell_model_all_df)
}
