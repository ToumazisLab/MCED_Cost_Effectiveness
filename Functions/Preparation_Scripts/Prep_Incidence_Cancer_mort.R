library(dplyr)
library(stringr)
library(tidyr)

incidence_seer <- read.csv("Data/NEW_SEER/incidence_feb17_2025.csv")
cancer_mort_seer <- read.csv("Data/NEW_SEER/cancer_mortality_feb17_2025.csv")

incidence_seer <- incidence_seer %>% 
  rename("Cancer" = MCED_Cancers_2, "Age" = X40_85, "Stage" = AJCC_I_IV, "Rate" = Age.Adjusted.Rate) %>% 
  mutate(Cancer = tolower(Cancer), Rate = (Count/Population))

cancer_mort_seer <- cancer_mort_seer %>% 
  rename("Cancer" = MCED_Cancers, "Age" = X40_85, "Rate" = Age.Adjusted.Rate) %>% 
  mutate(Cancer = tolower(Cancer), Rate = (Count/Population))


expand_age_range <- function(age_range, max_age = 100) {
  # Check if the age range contains a "+"
  if (str_detect(age_range, "\\+")) {
    # Extract the starting age (e.g., "85" from "85+ years")
    start <- as.numeric(str_extract(age_range, "\\d+"))
    # Create a sequence from the start to max_age
    return(seq(start, max_age))
  } else {
    # Extract the two numbers (e.g., "40" and "44" from "40-44 years")
    numbers <- str_extract_all(age_range, "\\d+")[[1]]
    start <- as.numeric(numbers[1])
    end <- as.numeric(numbers[2])
    return(seq(start, end))
  }
}

incidence_seer <- incidence_seer %>%
  mutate(Age_seq = lapply(Age, expand_age_range)) %>%  
  unnest(Age_seq) %>% 
  select(-Age) %>% 
  rename(Age = Age_seq)

df_bins <- incidence_seer %>%
  mutate(rate    = Rate,
         bin_lo  = 5 * floor(Age/5),
         bin_hi  = bin_lo + 5) %>%
  group_by(Cancer, Stage, bin_lo, bin_hi) %>%
  summarize(rate = mean(rate), .groups="drop") %>%
  # turn each bin into two “knot” points, at lo & hi
  rowwise() %>%
  do(data.frame(
    Cancer = .$Cancer,
    Stage  = .$Stage,
    Age    = c(.$bin_lo, .$bin_hi),
    rate   = .$rate
  )) %>%
  arrange(Cancer, Stage, Age)

ages <- sort(unique(incidence_seer$Age))

# 2) fit an FMM spline (cubic) through those boundary knots
smooth_bound <- df_bins %>%
  group_by(Cancer, Stage) %>%
  do({
    f <- splinefun(.$Age, .$rate, method = "fmm")
    tibble(Age = ages,
           rate_smooth = f(ages))
  }) %>%
  ungroup()

ggplot() +
  # original “flat” rates (lower alpha so you can see the smoother on top)
  geom_line(data = df_bins,
            aes(Age, rate, color = Stage),
            alpha = 0.4) +
  # smoothed spline
  geom_line(data = smooth_bound,
            aes(Age, rate_smooth, color = Stage),
            size  = 1) +
  facet_wrap(~ Cancer, scales = "free_y") +
  labs(x = "Age",
       y = "Incidence rate per 1e6",
       color = "Stage") +
  theme_minimal()

smooth_bound <- smooth_bound %>% 
  mutate(Rate = rate_smooth)

cancer_mort_seer <- cancer_mort_seer %>%
  mutate(Age_seq = lapply(Age, expand_age_range), "probability" = Rate) %>%  
  unnest(Age_seq) %>% 
  select(-Age) %>% 
  rename(Age = Age_seq)


# write.csv(smooth_bound, "/Users/rsnargund/Documents/MCED_Cohort_Model/FINAL_MCED_Model/Processed_Data/overall/incidence_final.csv", row.names = FALSE)
# write.csv(cancer_mort_seer, "/Users/rsnargund/Documents/MCED_Cohort_Model/FINAL_MCED_Model/Processed_Data/raw/cancer_mort.csv", row.names = FALSE)

processCancerData <- function(cancer_mort_seer, all_cause, cancer_list) {
  # Step 1: Filter rows where the 'Cancer' column is in the provided cancer_list
  filtered_cancer <- subset(cancer_mort_seer, Cancer %in% cancer_list)
  
  # Step 2: Aggregate the cancer mortality probabilities by Age
  total_cancer_mortality_by_age <- aggregate(probability ~ Age, 
                                             data = filtered_cancer, 
                                             sum)
  colnames(total_cancer_mortality_by_age)[2] <- "Total_cancer_probability"
  
  # Step 3: Merge the aggregated cancer mortality data with the all-cause mortality data
  all_cancer_mort <- merge(all_cause, total_cancer_mortality_by_age, 
                           by = "Age", all.x = TRUE)
  
  # Handle ages with no cancer data by replacing NA values with 0
  all_cancer_mort$Total_cancer_probability[is.na(all_cancer_mort$Total_cancer_probability)] <- 0
  
  # Calculate the non-cancer all-cause mortality probability
  all_cancer_mort$Non_cancer_all_cause_probability <- 
    all_cancer_mort$All_cause_mortality - all_cancer_mort$Total_cancer_probability
  
  # Return the resulting data frame
  return(all_cancer_mort)
}
all_cancer_mort_default9 <- processCancerData(cancer_mort_seer, all_cause, cancer_types_MDA)
all_cancer_mort_overall <- processCancerData(cancer_mort_seer, all_cause, cancer_types)

# write.csv(all_cancer_mort_overall, "/Users/rsnargund/Documents/MCED_Cohort_Model/FINAL_MCED_Model/Processed_Data/overall/all_cancer_mort_overall.csv", row.names = FALSE)
# write.csv(all_cancer_mort_default9, "/Users/rsnargund/Documents/MCED_Cohort_Model/FINAL_MCED_Model/Processed_Data/default9/all_cancer_mort_default9.csv", row.names = FALSE)




survival_raw <- read_csv("Data/NEW_SEER/Suvival_feb18.csv")
# Update survival_raw with a new column 'CancerGroup'
survival_raw <- survival_raw %>%
  mutate(CancerGroup = case_when(
    # head and neck
    Cancer %in% c("Nose, Nasal Cavity and Middle Ear", "Larynx", "Lip", "Tongue",
                  "Salivary Gland", "Floor of Mouth", "Gum and Other Mouth", "Nasopharynx",
                  "Tonsil", "Oropharynx", "Hypopharynx", "Other Oral Cavity and Pharynx") ~ "head and neck",
    
    # colon and rectum
    Cancer %in% c("Cecum", "Appendix", "Ascending Colon", "Hepatic Flexure", "Transverse Colon",
                  "Splenic Flexure", "Descending Colon", "Sigmoid Colon", "Large Intestine, NOS",
                  "Rectosigmoid Junction", "Rectum") ~ "colon and rectum",
    
    # anus
    Cancer == "Anus" ~ "anus",
    
    # liver and intrahepatic bile duct
    Cancer %in% c("Liver", "Intrahepatic Bile Duct") ~ "liver and intrahepatic bile duct",
    
    # lung and bronchus
    Cancer == "Lung and Bronchus" ~ "lung and bronchus",
    
    # breast
    Cancer == "Breast" ~ "breast",
    
    # cervix
    Cancer == "Cervix" ~ "cervix",
    
    # uterus
    Cancer %in% c("Corpus Uteri", "Uterus, NOS") ~ "uterus",
    
    # ovarian
    Cancer == "Ovarian" ~ "ovarian",
    
    # prostate
    Cancer == "Prostate" ~ "prostate",
    
    # pancreas
    Cancer == "Pancreas" ~ "pancreas",
    
    # kidney and renal pelvis
    Cancer == "Kidney and Renal Pelvis" ~ "kidney and renal pelvis",
    
    # lymphoma
    Cancer %in% c("Hodgkin - Nodal", "Hodgkin - Extranodal", "NHL - Nodal", "NHL - Extranodal") ~ "lymphoma",
    
    # bladder
    Cancer == "Bladder" ~ "bladder",
    
    # esophagus
    Cancer == "Esophagus" ~ "esophagus",
    
    # stomach
    Cancer == "Stomach" ~ "stomach",
    
    # default case
    TRUE ~ "Other"
  ))

aggregated_survival <- survival_raw %>%
  mutate(survival_year = Survival_months / 12) %>% 
  group_by(CancerGroup, Stage, Age) %>%
  summarise(
    avg_survival_yr    = mean(survival_year, na.rm = TRUE),
    median_survival_yr = median(survival_year, na.rm = TRUE),
    avg_survival_mth   = mean(Survival_months, na.rm = TRUE),
    median_survival_mth= median(Survival_months, na.rm = TRUE),
    n                  = n(),
    .groups            = "drop"
  )


weighted_survival <- aggregated_survival %>%
  mutate(Stage_simple = case_when(
    grepl("^IV", Stage) ~ "IV",
    grepl("^III", Stage) ~ "III",
    grepl("^II", Stage) ~ "II",
    grepl("^I", Stage) ~ "I",
    TRUE ~ NA_character_
  )) %>%
  group_by(CancerGroup, Age, Stage_simple) %>%
  summarise(
    weighted_avg_survival_mth = weighted.mean(avg_survival_mth, w = n, na.rm = TRUE),
    total_count             = sum(n),
    .groups                 = "drop"
  )

final_survival <- weighted_survival %>% 
  rename(Stage = Stage_simple, Cancer = CancerGroup, expected_survival = weighted_avg_survival_mth)  %>%
  mutate(expected_survival = expected_survival/12) %>% 
  select(Age, Cancer, Stage, expected_survival)


# Nest data by Age and Cancer
nested_data <- final_survival_cum %>%
  group_by(Age, Cancer) %>%
  nest()

# Fit the exponential model
fitted_models <- nested_data %>%
  mutate(model = map(data, ~ tryCatch(
    nls(cumulative_expected_survival ~ S0 * exp(-lambda * Stage_num),
        start = list(S0 = max(.x$cumulative_expected_survival, na.rm = TRUE), lambda = 0.1),
        data = .x),
    error = function(e) NULL  # Return NULL if the model fails
  )))

# Generate predictions if the model was successfully fit
fitted_data <- fitted_models %>%
  mutate(predictions = map2(data, model, ~ if (!is.null(.y)) {
    .x %>%
      mutate(predicted = predict(.y, newdata = .x))
  } else {
    .x %>% mutate(predicted = NA_real_) 
  })) %>%
  unnest(predictions)


survival_plot <- fitted_data %>% 
  filter(Cancer == "lung and bronchus", Age %in% c("50-54 years", "55-59 years"))



ggplot(survival_plot, aes(x = Stage_num)) +
  geom_point(aes(y = cumulative_expected_survival, color = "Actual"), size = 2) +
  geom_line(aes(y = predicted, color = "Fitted"), size = 1) +
  ylim(0,20) +
  facet_wrap(~ Age + Cancer, scales = "free_y") + 
  labs(title = "Exponential Decay Fit: Actual vs. Predicted",
       x = "Cancer Stage", 
       y = "Cumulative Expected Survival") +
  scale_color_manual(values = c("Actual" = "red", "Fitted" = "blue")) +
  theme_minimal()



