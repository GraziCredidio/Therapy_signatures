# Modeling EZE cohort
# MaAslin2 (gene expression model): Prednisolon x No systemic therapy
# Inactive disease + crp < 5

# graphics.off()
# rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("D:\\Documentos/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

library(tidyverse)
library(Maaslin2)


coldata_pred_noSyst <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_pred_R_crp_04.05.txt", sep = "\t")
vst_pred_noSyst <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_pred_R_crp_04.05.txt", sep = "\t")

# setting no_syst as the reference value
unique(coldata_pred_noSyst$pred_vs_noSyst)
coldata_pred_noSyst$pred_vs_noSyst_mod = coldata_pred_noSyst$pred_vs_noSyst
coldata_pred_noSyst$pred_vs_noSyst_mod[coldata_pred_noSyst$pred_vs_noSyst_mod == "no_syst"] =
  "a_noSyst"
coldata_pred_noSyst$pred_vs_noSyst_mod[coldata_pred_noSyst$pred_vs_noSyst_mod == "pred"] =
  "z_pred"


coldata_pred_noSyst <- coldata_pred_noSyst %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, pred_vs_noSyst_mod), as.factor))


dir.create("Output_files/Maaslin2/Results/Remission_crp/pred_vs_noSyst_crp")

fit_data_pred_noSyst = Maaslin2( 
  input_data = vst_pred_noSyst, 
  input_metadata = coldata_pred_noSyst, 
  analysis_method = "LM",
  normalization = "NONE",
  transform = "NONE",
  output = "Output_files/Maaslin2/Results/Remission_crp/pred_vs_noSyst_crp", 
  fixed_effects = c("pred_vs_noSyst_mod", "crp_log", "sex", "biologics"),
  random_effects = c("diagnosis_class", "age_group2", "bmi_class"))

# Generating tables
maaslin2_all_results_pred_noSyst <- fit_data_pred_noSyst$results
maaslin2_results_pred_noSyst <- maaslin2_all_results_pred_noSyst %>% filter(metadata == 'pred_vs_noSyst_mod') 
maaslin2_results_pred_noSyst$qval <- p.adjust(maaslin2_results_pred_noSyst$pval, method = 'BH')

write.table(maaslin2_results_pred_noSyst, "Output_files/Maaslin2/Results/Remission_crp/pred_vs_noSyst_crp/maaslin2_results_pred_vs_noSyst_R_crp.txt", sep = "\t",  quote = FALSE)

