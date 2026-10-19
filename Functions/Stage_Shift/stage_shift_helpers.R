isotone_fix <- function(sens, Stage, num) {
  out_val <- sens
  ndx <- match(Stage, Stage)  # Match stages
  good_ndx <- ndx[!is.na(ndx)]  # Keep valid indices
  
  if (length(good_ndx) > 1) {
    y <- sens[good_ndx]
    w <- num[good_ndx]
    val <- pava(y, w)  # Apply isotonic regression
    out_val[good_ndx] <- val  # Update sensitivity with isotonic regression values
  }
  
  out_val
}

name_mapping <- c(
  "Anus" = "anus",
  "Bladder" = "bladder",
  "Breast" = "breast",
  "Cervix" = "cervix",
  "Colon/Rectum" = "colon and rectum",
  "Esophagus" = "esophagus",
  "Head and Neck" = "head and neck",
  "Kidney" = "kidney and renal pelvis",
  "Liver/Bile-duct" = "liver and intrahepatic bile duct",
  "Lung" = "lung and bronchus",
  "Lymphoma" = "lymphoma",
  "Ovary" = "ovarian",
  "Pancreas" = "pancreas",
  "Prostate" = "prostate",
  "Stomach" = "stomach",
  "Uterus" = "uterus",
  "Urothelial Tract" = "urothelial"
)

# MISSING:
# Gallbladder, Melanoma, Sarcoma, Thyroid, Urothelial Tract, 
# Lymphoid Leukemia, Myeloid Neoplasm, Plasma Cell Neoplasm 

apply_mapping_and_filter <- function(df, name_mapping, cancer_list) {
  df$Cancer <- name_mapping[df$Cancer]
  df <- df[df$Cancer %in% cancer_list, ]
  print("mapping")
  print(paste(unique(df$Cancer)))
  return(df)
}

# dwell_model_timing_named <- apply_mapping_and_filter(dwell_model_group_df, name_mapping, cancer_types)
