def_param_list <- list(
  cycle_length = 1,
  n_age_init = 39,
  n_age_max = 100,
  screening_cancers = c("colon and rectum", 
                        "lung and bronchus", 
                        "breast", "cervix", "prostate"),
  num_people = 1e6,
  female_prop = 0.5101351, # https://www.census.gov/data/tables/2017/demo/age-and-sex/2017-age-sex-composition.html Table 1
  male_prop = 0.4898649,
  female_cancers = c("breast", 
                     "cervix", 
                     "ovarian",
                     "uterus"),
  male_cancers = c("prostate"),
  v_names_states = c("H", "I", "II", "III", "IV", "D"),
  xStage = c("I","II","III","IV")
)

n_states = length(def_param_list$v_names_states)

def_param_list[["n_states"]] = c(def_param_list[["n_states"]], n_states)

cancer_types <- cancer_types_overall <- c(
  "anus",
  "bladder",
  "breast",
  "cervix",
  "colon and rectum",
  "esophagus",
  "head and neck",
  "kidney and renal pelvis",
  "liver and intrahepatic bile duct",
  "lung and bronchus",
  "lymphoma",
  "ovarian",
  "pancreas",
  "prostate",
  "stomach",
  "uterus"
)

cancer_types_MDA <- c(
  "colon and rectum",
  "lung and bronchus",
  "breast",
  "pancreas",
  "ovarian",
  "esophagus",
  "liver and intrahepatic bile duct",
  "stomach",
  "prostate"    
)
