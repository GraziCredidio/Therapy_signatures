# EZE cohort: Therapy Signatures
  # Linear Mixed Model
  # Prednisolone vs No systemic therapy
  # Inactive and all patients models
  # Author: Graziella Credidio

rm(list = ls())

folders <- c("Output_files/Maaslin2/pred_vs_noSyst/inactive", "Output_files/Maaslin2/pred_vs_noSyst/allPatients")
for (i in folders){
  if (!dir.exists(i)) {
    dir.create(i)}
}

# Loading packages ----
library(tidyverse)
library(Maaslin2)

# Loading data ----
# Inactive disease patients
coldata_pred_noSyst <- read.table("Cleaned_tables/models/pred/EZECohort_coldata_pred.txt", sep = "\t")
vst_pred_noSyst <- read.table("Cleaned_tables/models/pred/vst_counts_pred.txt", sep = "\t")

# All patients
all_coldata_pred_noSyst <- read.table("Cleaned_tables/models/pred/EZECohort_coldata_pred_allPatients.txt", sep = "\t")
all_vst_pred_noSyst <- read.table("Cleaned_tables/models/pred/vst_counts_pred_allPatients.txt", sep = "\t")

# Coldata preprocessing ----
preprocessing_coldata <- function(coldata){
  # Setting reference values
  coldata$pred_vs_noSyst_mod = coldata$pred_vs_noSyst
  coldata$pred_vs_noSyst_mod[coldata$pred_vs_noSyst_mod == "no_syst"] = "a_noSyst"
  coldata$pred_vs_noSyst_mod[coldata$pred_vs_noSyst_mod == "pred"] = "z_pred"
  
  coldata <- coldata %>% 
    mutate(across(c(diagnosis_class, age_group,  
                    sex, bmi_class, pred_vs_noSyst_mod), as.factor))
  
  return(coldata)
}

coldata_pred_noSyst <- preprocessing_coldata(coldata_pred_noSyst)  
all_coldata_pred_noSyst <- preprocessing_coldata(all_coldata_aza_noSyst) 

# Maaslin2 model ----
fit_model_pred <- function(coldata, vst_counts, outputPath, resultsPath){
  model = Maaslin2( 
    input_data = vst_counts, 
    input_metadata = coldata, 
    analysis_method = "LM",
    normalization = "NONE",
    transform = "NONE",
    output = outputPath,
    fixed_effects = c("pred_vs_noSyst_mod", "crp_log", "sex", "biologics"),
    random_effects = c("diagnosis_class", "age_group", "bmi_class"))
  
  # Generating tables
  maaslin2_all_results_pred_noSyst <- model$results
  maaslin2_results_pred_noSyst <- maaslin2_all_results_pred_noSyst %>% filter(metadata == 'pred_vs_noSyst_mod') 
  maaslin2_results_pred_noSyst$qval <- p.adjust(maaslin2_results_pred_noSyst$pval, method = 'BH')
  
  write.table(maaslin2_results_pred_noSyst, resultsPath, sep = "\t",  quote = FALSE)
}


fit_model_pred(coldata = coldata_pred_noSyst, vst_counts = vst_pred_noSyst, 
              outputPath = "Output_files/Maaslin2/pred_vs_noSyst/inactive",
              resultsPath = "Output_files/Maaslin2/results/inactive/maaslin2_results_pred_vs_noSyst.txt")

fit_model_pred(coldata = all_coldata_pred_noSyst, vst_counts = all_vst_pred_noSyst, 
              outputPath = "Output_files/Maaslin2/pred_vs_noSyst/allPatients",
              resultsPath = "Output_files/Maaslin2/results/allPatients/maaslin2_results_pred_vs_noSyst_allPatients.txt")