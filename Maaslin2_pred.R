# EZE cohort: Therapy Signatures
  # Linear Mixed Model
  # Prednisolone vs No Systemic Therapies
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(tidyverse)
library(Maaslin2)

# Loading data ----
coldata_pred_noSyst <- read.table("Cleaned_tables/models/pred/EZECohort_coldata_pred.txt", sep = "\t")
vst_pred_noSyst <- read.table("Cleaned_tables/models/pred/vst_counts_pred.txt", sep = "\t")

# Setting no_syst as the reference value ----
coldata_pred_noSyst$pred_vs_noSyst_mod = coldata_pred_noSyst$pred_vs_noSyst
coldata_pred_noSyst$pred_vs_noSyst_mod[coldata_pred_noSyst$pred_vs_noSyst_mod == "no_syst"] =
  "a_noSyst"
coldata_pred_noSyst$pred_vs_noSyst_mod[coldata_pred_noSyst$pred_vs_noSyst_mod == "pred"] =
  "z_pred"

coldata_pred_noSyst <- coldata_pred_noSyst %>% 
  mutate(across(c(diagnosis_class, age_group,  
                  sex, bmi_class, pred_vs_noSyst_mod), as.factor))

dir.create("Output_files/Maaslin2/pred_vs_noSyst")

fit_data_pred_noSyst = Maaslin2( 
  input_data = vst_pred_noSyst, 
  input_metadata = coldata_pred_noSyst, 
  analysis_method = "LM",
  normalization = "NONE",
  transform = "NONE",
  output = "Output_files/Maaslin2/pred_vs_noSyst", 
  fixed_effects = c("pred_vs_noSyst_mod", "crp_log", "sex", "biologics"),
  random_effects = c("diagnosis_class", "age_group", "bmi_class"))

# Generating tables
maaslin2_all_results_pred_noSyst <- fit_data_pred_noSyst$results
maaslin2_results_pred_noSyst <- maaslin2_all_results_pred_noSyst %>% filter(metadata == 'pred_vs_noSyst_mod') 
maaslin2_results_pred_noSyst$qval <- p.adjust(maaslin2_results_pred_noSyst$pval, method = 'BH')

write.table(maaslin2_results_pred_noSyst, "Output_files/Maaslin2/maaslin2_results_pred_vs_noSyst.txt", sep = "\t",  quote = FALSE)