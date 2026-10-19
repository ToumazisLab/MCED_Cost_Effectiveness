library(isotone)
library(dplyr)

estimate_isotonic <- function(data) {
  # Ensure data is ordered by Age
  data <- data %>% arrange(Age)
  
  # Apply isotonic regression with decreasing constraint
  fit <- gpava(data$Age, data$Survival_years, weights = rep(1, nrow(data)))
  
  # Extract fitted values ensuring they follow an isotonic decay
  data$Survival_years <- fit$x
  
  # Predict values for Age 40 and 45
  age_40 <- approx(data$Age, data$Survival_years, xout = 40, rule = 2)$y
  age_45 <- approx(data$Age, data$Survival_years, xout = 45, rule = 2)$y
  
  # Ensure monotonicity across stages (Stage I ≥ Stage II ≥ Stage III ≥ Stage IV)
  data <- data %>%
    bind_rows(data.frame(Age = c(40, 45), Survival_years = c(age_40, age_45))) %>%
    arrange(Age, Stage_num)
  
  return(data)
}

df <- data_list$survival_expanded %>% 
  rename("Survival_years" = Count) %>% 
  mutate(Stage_num = case_when(
    Stage == "I"  ~ 1,
    Stage == "II" ~ 2,
    Stage == "III" ~ 3,
    Stage == "IV" ~ 4,
    TRUE ~ NA_real_
  ))  %>% 
  mutate(Age_interval = case_when(
    Age %in% c(50:54) ~ "50-54",
    Age %in% c(55:59) ~ "55-59",
    Age %in% c(60:64) ~ "60-64",
    Age %in% c(65:69) ~ "65-69",
    Age %in% c(70:74) ~ "70-74",
    Age %in% c(75:79) ~ "75-79",
    Age %in% c(80:84) ~ "80-84",
    Age %in% c(85:100) ~ "85+"
  ))


df_unique <- df %>%
  group_by(Cancer, Stage, Stage_num, Age_interval) %>%
  summarise(
    Age_val = min(Age),
    Survival_years = first(Survival_years),
    .groups = "drop"
  ) %>%
  arrange(Cancer, Stage_num, Age_val)

df_unique
# (b) For each Cancer–Stage group, fit an exponential model and predict for ages 40 and 45.
exp_predictions <- df_unique %>%
  group_by(Cancer, Stage, Stage_num) %>%
  group_modify(~ {
    if(nrow(.x) < 2) {
      # If there is only one age group, fallback to the same value (though this loses the trend)
      tibble(
        Age_val = c(40, 45),
        Survival_years = rep(.x$Survival_years[1], 2),
        Age_interval = c("40-44", "45-49")
      )
    } else {
      # Fit a linear model to log(Survival_years)
      mod <- lm(log(Survival_years) ~ Age_val, data = .x)
      newAges <- data.frame(Age_val = c(40, 45))
      # Predict on the log-scale, then exponentiate.
      preds <- exp(predict(mod, newdata = newAges))
      tibble(
        Age_val = c(40, 45),
        Survival_years = preds,
        Age_interval = c("40-44", "45-49")
      )
    }
  }) %>%
  ungroup()

# (c) Combine with your original (collapsed) data if desired:
combined_exp <- bind_rows(df_unique, exp_predictions) %>% 
  arrange(Cancer, Stage_num, Age_val)

survival_int <- combined_exp %>%
  mutate(Age_seq = lapply(Age_interval, expand_age_range)) %>%  
  unnest(Age_seq) %>% 
  select(-Age_val) %>% 
  rename(Age = Age_seq)

write.csv(survival_int, "FINAL_MCED_Model/Processed_Data/raw/Survival_dummy.csv",row.names = FALSE)

transposed_survival_int <- survival_int %>%
  mutate(Stage_Cancer = paste(Stage, Cancer, sep = "_")) %>%
  select(Age, Stage_Cancer, Survival_years) %>%
  pivot_wider(names_from = Stage_Cancer, values_from = Survival_years) %>%
  arrange(Age) %>% 
  select(-Age)


write.csv(transposed_survival_int, "FINAL_MCED_Model/Processed_Data/raw/transposed_survival_dummy.csv",row.names = FALSE)


