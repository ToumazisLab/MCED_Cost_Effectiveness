##
## 1. MCED Cost
##
mced_cost_order <- c(
  # Red
  "FP_Imaging_cost", "FP_Diagnostic_cost", "FP_Invasive_cost",
  # Green
  "TN_Diagnostic_cost", "TN_Imaging_cost", "TN_Invasive_cost",
  # Black 
  "cancer_nc_cost",
  # Blue
  "tn_fp_mced_test_cost", "cancer_test_cost",
  # Pink
  "cancer_diagnosed_cost",
  "cancer_incidentally_diagnosed_cost",
  "cancer_screen_diagnosed_cost",
  # Orange shade 1
  "I_pd_cost_TP", "I_pd_cost_NA", "I_pd_cost_FN",
  # Orange shade 2
  "II_pd_cost_TP", "II_pd_cost_NA", "II_pd_cost_FN",
  # Orange shade 3
  "III_pd_cost_TP", "III_pd_cost_NA", "III_pd_cost_FN",
  # Orange shade 4
  "IV_pd_cost_TP", "IV_pd_cost_NA", "IV_pd_cost_FN"
)

##
## 2. MCED Utility
##
# - "diagnostic_disutility" -> light green
# - I_pd_utility for TP/FN/NA -> orange shade 1
# - II_pd_utility for TP/FN/NA -> orange shade 2
# - III & IV_pd_utility for TP/FN/NA -> orange shade 3 & 4
mced_utility_order <- c(
  "adjusted_utility",
  # I
  "I_pd_utility_TP", "I_pd_utility_NA", "I_pd_utility_FN",
  # II
  "II_pd_utility_TP", "II_pd_utility_NA", "II_pd_utility_FN",
  # III
  "III_pd_utility_TP", "III_pd_utility_NA", "III_pd_utility_FN",
  # IV
  "IV_pd_utility_TP", "IV_pd_utility_NA", "IV_pd_utility_FN"
)

##
## 3. Usual Utility
##
# - "adjusted_utility" -> light green
# - I_pd_utility -> orange shade 1
# - II_pd_utility -> orange shade 2
# - III & IV_pd_utility -> orange shade 3 & 4
#
usual_utility_order <- c(
  "adjusted_utility",
  "I_pd_utility",
  "II_pd_utility",
  "III_pd_utility",
  "IV_pd_utility"
)

##
## 4. Usual Cost
##
# - FP_Imaging_Cost, FP_Invasive_cost, FP_Diagnostic_cost -> Red
# - TN_Imaging_cost, TN_Diagnostic_cost -> Green
# - cancer_incidentally_diagnosed, cancer_screen_diagnosed -> Pink
# - I–IV_pd_cost (with TP/NA/FN) -> shades of orange
usual_cost_order <- c(
  # Red
  "FP_Imaging_cost", "FP_Diagnostic_cost", "FP_Invasive_cost",
  # Green
  "TN_Imaging_cost", "TN_Diagnostic_cost", 
  # Pink
  "cancer_incidentally_diagnosed_cost", "cancer_screen_diagnosed_cost",
  # Orange (stage I–IV, each with TP, NA, FN)
  "I_pd_cost", "II_pd_cost", "III_pd_cost", "IV_pd_cost"
)


var_colors <- c(
  ################################
  ## MCED COST
  ################################
  # Red
  "FP_Imaging_cost"    = "#efc0c2",
  "FP_Diagnostic_cost" = "#dc7175",
  "FP_Invasive_cost"   = "#be2e33",
  
  # Green
  "TN_Imaging_cost"    = "#db94ba",
  "TN_Diagnostic_cost" = "#C04185",
  "TN_Invasive_cost"   = "#7a2954",
  
  # Black 
  "cancer_nc_cost" = "black", 
  
  # Blue
  "tn_fp_mced_test_cost" = "#B0B096",
  "cancer_test_cost"     = "#747458",
  
  # Pink
  "cancer_diagnosed_cost"              = "#bdd8ef",
  "cancer_incidentally_diagnosed_cost" = "#4893d3",
  "cancer_screen_diagnosed_cost"       = "#1c4a73",
  
  # Stage I
  "I_pd_cost_TP" = "#EFEDBD",
  "I_pd_cost_NA" = "#D3c348",
  "I_pd_cost_FN" = "#838020",
  
  # Stage II
  "II_pd_cost_TP" = "#EFE1bd",
  "II_pd_cost_NA" = "#d3ab48",
  "II_pd_cost_FN" = "#735a1c",
  
  # Stage III
  "III_pd_cost_TP" = "#EBCAAD",
  "III_pd_cost_NA" = "#D38948",
  "III_pd_cost_FN" = "#A46228",
  
  # Stage IV
  "IV_pd_cost_TP" = "#945a24",
  "IV_pd_cost_NA" = "#73461c",
  "IV_pd_cost_FN" = "#422810",
  
  ################################
  ## MCED UTILITY
  ################################
  "adjusted_utility" = "#48D389",
  
  # Stage I
  "I_pd_utility_TP" = "#EFEDBD",
  "I_pd_utility_NA" = "#D3c348",
  "I_pd_utility_FN" = "#838020",
  
  # Stage II
  "II_pd_utility_TP" = "#EFE1bd",
  "II_pd_utility_NA" = "#d3ab48",
  "II_pd_utility_FN" = "#735a1c",
  
  # Stage III
  "III_pd_utility_TP" = "#EBCAAD",
  "III_pd_utility_NA" = "#D38948",
  "III_pd_utility_FN" = "#A46228",
  
  # Stage IV
  "IV_pd_utility_TP" = "#945a24",
  "IV_pd_utility_NA" = "#73461c",
  "IV_pd_utility_FN" = "#422810",
  
  ################################
  ## USUAL UTILITY
  ################################
  "adjusted_utility" = "#48D389",
  
  # (Same color logic)
  "I_pd_utility"   = "#D3c348",
  "II_pd_utility"  = "#d3ab48",
  "III_pd_utility" = "#D38948",
  "IV_pd_utility"  = "#73461c",
  
  ################################
  ## USUAL COST
  ################################
  "FP_Imaging_Cost"      = "#efc0c2",
  "FP_Diagnostic_cost"    = "#dc7175",
  "FP_Invasive_cost"     = "#be2e33",
  
  "TN_Imaging_cost"      = "#db94ba",
  "TN_Diagnostic_cost"    = "#C04185",
  
  "cancer_incidentally_diagnosed_cost" = "#4893d3",
  "cancer_screen_diagnosed_cost"       = "#1c4a73",
  
  # Stages I–IV again in orange
  "I_pd_cost"  = "#D3c348",
  "II_pd_cost" = "#d3ab48",
  "III_pd_cost" = "#D38948",
  "IV_pd_cost"  = "#73461c"
)


