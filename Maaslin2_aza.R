# EZE cohort: Therapy Signatures
  # Linear Mixed Model
  # Azathioprine vs No Systemic Therapies
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(tidyverse)
library(Maaslin2)

# Loading data ----
coldata_aza_noSyst <- read.table("Cleaned_tables/models/aza/EZECohort_coldata_aza.txt", sep = "\t")
vst_aza_noSyst <- read.table("Cleaned_tables/models/aza/vst_counts_aza.txt", sep = "\t")

# Setting no_syst as the reference value ----
coldata_aza_noSyst$aza_vs_noSyst_mod = coldata_aza_noSyst$aza_vs_noSyst
coldata_aza_noSyst$aza_vs_noSyst_mod[coldata_aza_noSyst$aza_vs_noSyst_mod == "no_syst"] =
  "a_noSyst"
coldata_aza_noSyst$aza_vs_noSyst_mod[coldata_aza_noSyst$aza_vs_noSyst_mod == "aza"] =
  "z_aza"

coldata_aza_noSyst <- coldata_aza_noSyst %>% 
  mutate(across(c(diagnosis_class, age_group,
                  sex, bmi_class, aza_vs_noSyst_mod), as.factor))


dir.create("Output_files/Maaslin2/aza_vs_noSyst")

fit_data_aza_noSyst = Maaslin2( 
  input_data = vst_aza_noSyst, 
  input_metadata = coldata_aza_noSyst, 
  analysis_method = "LM",
  normalization = "NONE",
  transform = "NONE",
  output ="Output_files/Maaslin2/aza_vs_noSyst",
  fixed_effects = c("aza_vs_noSyst_mod", "crp_log", "sex", "biologics", "diagnosis_class"),
  random_effects = c("age_group", "bmi_class"))

# Generating tables
maaslin2_all_results_aza_noSyst <- fit_data_aza_noSyst$results
maaslin2_results_aza_noSyst <- maaslin2_all_results_aza_noSyst %>% filter(metadata == 'aza_vs_noSyst_mod') 
maaslin2_results_aza_noSyst$qval <- p.adjust(maaslin2_results_aza_noSyst$pval, method = 'BH')

write.table(maaslin2_results_aza_noSyst, "Output_files/Maaslin2/aza_vs_noSyst/maaslin2_results_aza_vs_noSyst.txt", sep = "\t",  quote = FALSE)
