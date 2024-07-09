# EZE cohort: Therapy Signatures
  # Linear Mixed Model
  # AntiTNF vs No Biologics
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(tidyverse)
library(Maaslin2)

# Loading data ----
coldata_noBio <- read.table("Cleaned_tables/models/antiTNF/EZECohort_coldata_antiTNF.txt", sep = "\t")
vst_noBio <- read.table("Cleaned_tables/models/antiTNF/vst_counts_antiTNF.txt", sep = "\t")

# Setting no_biologics as the reference value ----
coldata_noBio$antiTNF_vs_noBiologics_mod = coldata_noBio$antiTNF_vs_noBiologics
coldata_noBio$antiTNF_vs_noBiologics_mod[coldata_noBio$antiTNF_vs_noBiologics_mod == "no_biologics"] =
  "a_noBiologics"
coldata_noBio$antiTNF_vs_noBiologics_mod[coldata_noBio$antiTNF_vs_noBiologics_mod == "anti_tnf"] =
  "z_antiTnf"

coldata_noBio <- coldata_noBio %>% 
  mutate(across(c(diagnosis_class, age_group, sex, bmi_class, Prednisolon, antiTNF_vs_noBiologics_mod), as.factor))

dir.create("Output_files/Maaslin2/antiTNF_vs_noBiologics")

fit_data_noBio = Maaslin2( 
  input_data = vst_noBio, 
  input_metadata = coldata_noBio, 
  analysis_method = "LM",
  normalization = "NONE",
  transform = "NONE",
  output = "Output_files/Maaslin2/antiTNF_vs_noBiologics", 
  fixed_effects = c("antiTNF_vs_noBiologics_mod", "crp_log", "sex", "Prednisolon"),
  random_effects = c("diagnosis_class", "age_group", "bmi_class"),
  reference = c("antiTNF_vs_noBiologics,no_biologics"))

maaslin2_all_results_noBio <- fit_data_noBio$results
maaslin2_results_noBio <- maaslin2_all_results_noBio %>% filter(metadata == 'antiTNF_vs_noBiologics_mod')
maaslin2_results_noBio$qval <- p.adjust(maaslin2_results_noBio$pval, method = 'BH')

write.table(maaslin2_results_noBio, "Output_files/Maaslin2/maaslin2_results_antiTNF_vs_noBio.txt", sep = "\t",  quote = FALSE)