# EZE cohort: Therapy Signatures
  # Linear Mixed Model
  # AntiTNF vs No Biologics
  # Inactive and all patients models
  # Author: Graziella Credidio

rm(list = ls())

folders <- c("Output_files/Maaslin2/antiTNF_vs_noBio/inactive", "Output_files/Maaslin2/antiTNF_vs_noBio/allPatients")
for (i in folders){
  if (!dir.exists(i)) {
    dir.create(i)}
}

# Loading packages ----
library(tidyverse)
library(Maaslin2)

# Loading data ----
# Inactive disease patients
coldata_noBio <- read.table("Cleaned_tables/models/antiTNF/EZECohort_coldata_antiTNF.txt", sep = "\t")
vst_noBio <- read.table("Cleaned_tables/models/antiTNF/vst_counts_antiTNF.txt", sep = "\t")

# All patients
all_coldata_noBio <- read.table("Cleaned_tables/models/antiTNF/EZECohort_coldata_antiTNF_allPatients.txt", sep = "\t")
all_vst_noBio <- read.table("Cleaned_tables/models/antiTNF/vst_counts_antiTNF_allPatients.txt", sep = "\t")

# Coldata preprocessing ----
preprocessing_coldata <- function(coldata){
  # Setting reference values
  coldata$antiTNF_vs_noBiologics_mod = coldata$antiTNF_vs_noBiologics
  coldata$antiTNF_vs_noBiologics_mod[coldata$antiTNF_vs_noBiologics_mod == "no_biologics"] = "a_noBiologics"
  coldata$antiTNF_vs_noBiologics_mod[coldata$antiTNF_vs_noBiologics_mod == "anti_tnf"] = "z_antiTnf"
  
  coldata <- coldata %>% 
    mutate(across(c(diagnosis_class, age_group, sex, 
                    bmi_class, Prednisolon, antiTNF_vs_noBiologics_mod), as.factor))
  
  return(coldata)
}

coldata_noBio <- preprocessing_coldata(coldata_noBio)  
all_oldata_noBio <- preprocessing_coldata(all_coldata_noBio) 

# Maaslin2 model ----
fit_model_antiTNF <- function(coldata, vst_counts, outputPath, resultsPath){
  model = Maaslin2( 
    input_data = vst_counts, 
    input_metadata = coldata, 
    analysis_method = "LM",
    normalization = "NONE",
    transform = "NONE",
    output = outputPath, 
    fixed_effects = c("antiTNF_vs_noBiologics_mod", "crp_log", "sex", "Prednisolon"),
    random_effects = c("diagnosis_class", "age_group", "bmi_class"),
    reference = c("antiTNF_vs_noBiologics,no_biologics"))
  
  # Generating tables
  maaslin2_all_results_noBio <- model$results
  maaslin2_all_results_noBio <- maaslin2_all_results_noBio %>% filter(metadata == 'antiTNF_vs_noBiologics_mod') 
  maaslin2_all_results_noBio$qval <- p.adjust(maaslin2_all_results_noBio$pval, method = 'BH')
  
  write.table(maaslin2_all_results_noBio, resultsPath, sep = "\t",  quote = FALSE)
}

fit_model_antiTNF(coldata = coldata_noBio, vst_counts = vst_noBio, 
              outputPath = "Output_files/Maaslin2/antiTNF_vs_noBiologics/inactive",
              resultsPath = "Output_files/Maaslin2/results/inactive/maaslin2_results_antiTNF_vs_noBio.txt")

fit_model_antiTNF(coldata = all_coldata_noBio, vst_counts = all_vst_noBio, 
              outputPath = "Output_files/Maaslin2/antiTNF_vs_noBiologics/allPatients",
              resultsPath = "Output_files/Maaslin2/results/allPatients/maaslin2_results_antiTNF_vs_noBio_allPatients.txt")
