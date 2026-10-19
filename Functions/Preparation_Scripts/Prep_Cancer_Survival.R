survival_raw <- read_csv("Data/NEW_SEER/Survival_31925.csv")
#### Update survival_raw with a new column 'CancerGroup' ####
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
    
    # Acute Myeloid Leukemia
    Cancer == "Acute Myeloid Leukemia" ~ "AML",
    
    # default case
    TRUE ~ "Other"
  ))

aggregated_survival <- survival_raw %>%
  mutate(Survival_months = as.numeric(Survival_months)) %>% 
  mutate(survival_year =   Survival_months / 12) %>% 
  group_by(CancerGroup, Stage, Age_group) %>%
  summarise(
    avg_survival_yr    = mean(survival_year, na.rm = TRUE),
    median_survival_yr = median(survival_year, na.rm = TRUE),
    avg_survival_mth   = mean(Survival_months, na.rm = TRUE),
    median_survival_mth= median(Survival_months, na.rm = TRUE),
    n                  = n(),
    .groups            = "drop"
  )

aggregated_survival <- aggregated_survival %>%
  mutate(Stage_simple = case_when(
    grepl("^IV", Stage) ~ "IV",
    grepl("^III", Stage) ~ "III",
    grepl("^II", Stage) ~ "II",
    grepl("^I", Stage) ~ "I",
    TRUE ~ NA_character_
  )) %>%
  group_by(CancerGroup, Age, Stage_simple) %>%
  summarise(
    avg_survival_mth = mean(avg_survival_mth, na.rm = TRUE),
    total_count             = sum(n),
    .groups                 = "drop"
  ) %>% 
  mutate(
    avg_survival_yr = avg_survival_mth / 12
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



#### AML ####
survival_raw <- read_csv("Data/NEW_SEER/survival_AML.csv")
#### Update survival_raw with a new column 'CancerGroup' ####
survival_raw <- survival_raw %>%
  mutate(COD_group = case_when(
    # head and neck
    COD == "Alive" ~ "Alive",
    
    COD == "Acute Myeloid Leukemia" ~ "cancer",
    
    # default case
    TRUE ~ "Other"
  ))

aggregated_survival <- survival_raw %>%
  mutate(Survival_months = as.numeric(Survival_months)) %>% 
  mutate(survival_year =   Survival_months / 12) %>% 
  group_by(Age, COD_group) %>%
  summarise(
    avg_survival_yr    = mean(survival_year, na.rm = TRUE),
    median_survival_yr = median(survival_year, na.rm = TRUE),
    avg_survival_mth   = mean(Survival_months, na.rm = TRUE),
    median_survival_mth= median(Survival_months, na.rm = TRUE),
    n                  = n(),
    .groups            = "drop"
  )
