# EZE cohort: Therapy Signatures
  # Linear Mixed Model
  # Azathioprine vs No Systemic Therapies
  # Inactive and all patients models
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(tidyverse)
library(Maaslin2)

# Loading data ----
# Inactive disease patients
coldata_aza_noSyst <- read.table("Cleaned_tables/models/aza/EZECohort_coldata_aza.txt", sep = "\t")
vst_aza_noSyst <- read.table("Cleaned_tables/models/aza/vst_counts_aza.txt", sep = "\t")

# All patients
all_coldata_aza_noSyst <- read.table("Cleaned_tables/models/aza/EZECohort_coldata_aza_allPatients.txt", sep = "\t")
all_vst_aza_noSyst <- read.table("Cleaned_tables/models/aza/vst_counts_aza_allPatients.txt", sep = "\t")

# Coldata preprocessing ----
preprocessing_coldata <- function(coldata){
  # Setting reference values
  coldata$aza_vs_noSyst_mod = coldata$aza_vs_noSyst
  coldata$aza_vs_noSyst_mod[coldata$aza_vs_noSyst_mod == "no_syst"] = "a_noSyst"
  coldata$aza_vs_noSyst_mod[coldata$aza_vs_noSyst_mod == "aza"] = "z_aza"
  
  coldata <- coldata %>% 
    mutate(across(c(diagnosis_class, age_group,
                    sex, bmi_class, aza_vs_noSyst_mod), as.factor))
  
  return(coldata)
}

coldata_aza_noSyst <- preprocessing_coldata(coldata_aza_noSyst)  
all_coldata_aza_noSyst <- preprocessing_coldata(all_coldata_aza_noSyst) 

# Creating output directories ----
dir.create("Output_files/Maaslin2/aza_vs_noSyst/inactive")
dir.create("Output_files/Maaslin2/aza_vs_noSyst/allPatients")

# Maaslin2 model ----
fit_model_aza <- function(coldata, vst_counts, outputPath, resultsPath){
  model = Maaslin2( 
    input_data = vst_counts, 
    input_metadata = coldata, 
    analysis_method = "LM",
    normalization = "NONE",
    transform = "NONE",
    output = "outputPath",
    fixed_effects = c("aza_vs_noSyst_mod", "crp_log", "sex", "biologics", "diagnosis_class"),
    random_effects = c("age_group", "bmi_class"))
  
  # Generating tables
  maaslin2_all_results_aza_noSyst <- model$results
  maaslin2_results_aza_noSyst <- maaslin2_all_results_aza_noSyst %>% filter(metadata == 'aza_vs_noSyst_mod') 
  maaslin2_results_aza_noSyst$qval <- p.adjust(maaslin2_results_aza_noSyst$pval, method = 'BH')
  
  write.table(maaslin2_results_aza_noSyst, resultsPath, sep = "\t",  quote = FALSE)
}

fit_model_aza(coldata = coldata_aza_noSyst, vst_counts = vst_aza_noSyst, 
              outputPath = "Output_files/Maaslin2/aza_vs_noSyst/inactive",
              resultsPath = "Output_files/Maaslin2/results/inactive/maaslin2_results_aza_vs_noSyst.txt")

fit_model_aza(coldata = all_coldata_aza_noSyst, vst_counts = all_vst_aza_noSyst, 
              outputPath = "Output_files/Maaslin2/aza_vs_noSyst/allPatients",
              resultsPath = "Output_files/Maaslin2/results/allPatients/maaslin2_results_aza_vs_noSyst_allPatients.txt")
