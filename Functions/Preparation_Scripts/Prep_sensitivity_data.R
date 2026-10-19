df <- df_sensitivity
eps <- .Machine$double.eps
margin <- 1e-3   # 0.1% away from the boundary

df_adj <- df %>%
  rowwise() %>%
  mutate(
    Sens_base = mean(c(Galleri, Thunder, thromboSeq), na.rm = TRUE)
    # you could use median(...) instead if you want robustness to outliers
  ) %>%
  ungroup() %>%
  
  # 2) overwrite your Sensitivity column with that baseline
  mutate(Sensitivity = Sens_base) %>%
  
  # 3) turn Stage into an ordered integer (I→1, II→2, …)
  mutate(Stage_num = case_when(
    Stage == "I"   ~ 1L,
    Stage == "II"  ~ 2L,
    Stage == "III" ~ 3L,
    Stage == "IV"  ~ 4L,
    TRUE           ~ NA_integer_
  )) %>%
  group_by(Cancer) %>%
  arrange(Stage_num, .by_group = TRUE) %>%
  mutate(
    sens1 = pava(Sensitivity, decreasing = FALSE),
    # clamp into [margin, 1-margin]
    sens2 = pmin(pmax(sens1, margin), 1 - margin),
    # cap I/II, etc…
    sens3 = if_else(Stage_num <= 2,
                    pmin(sens2, 0.65 - margin),
                    sens2),
    Sensitivity_adj = pava(sens3, decreasing = FALSE)
  ) %>%
  ungroup() %>%
  
  # drop intermediate columns if you like
  select(-Sens_base, -Stage_num, -sens1, -sens2, -sens3) %>% 
  mutate(Sensitivity = Sensitivity_adj)
