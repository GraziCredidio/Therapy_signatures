# Modeling EZE cohort
# MaAslin2 (gene expression model): Azathioprin x No systemic therapy (only IBD patients)
# Inactive disease + crp < 5

graphics.off()
rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("D:\\Documentos/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

library(tidyverse)
library(Maaslin2)


coldata_aza_noSyst <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_aza_R_crp_04.05.txt", sep = "\t")
vst_aza_noSyst <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_aza_R_crp_04.05.txt", sep = "\t")


# setting no_syst as the reference value
unique(coldata_aza_noSyst$aza_vs_noSyst)
coldata_aza_noSyst$aza_vs_noSyst_mod = coldata_aza_noSyst$aza_vs_noSyst
coldata_aza_noSyst$aza_vs_noSyst_mod[coldata_aza_noSyst$aza_vs_noSyst_mod == "no_syst"] =
  "a_noSyst"
coldata_aza_noSyst$aza_vs_noSyst_mod[coldata_aza_noSyst$aza_vs_noSyst_mod == "aza"] =
  "z_aza"

coldata_aza_noSyst <- coldata_aza_noSyst %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, aza_vs_noSyst_mod), as.factor))


dir.create("Output_files/Maaslin2/Results/Remission_crp/aza_vs_noSyst_crp")

fit_data_aza_noSyst = Maaslin2( 
  input_data = vst_aza_noSyst, 
  input_metadata = coldata_aza_noSyst, 
  analysis_method = "LM",
  normalization = "NONE",
  transform = "NONE",
  output ="Output_files/Maaslin2/Results/Remission_crp/aza_vs_noSyst_crp",
  fixed_effects = c("aza_vs_noSyst_mod", "crp_log", "sex", "biologics", "diagnosis_class"),
  random_effects = c("age_group2", "bmi_class"))

# Generating tables
maaslin2_all_results_aza_noSyst <- fit_data_aza_noSyst$results
maaslin2_results_aza_noSyst <- maaslin2_all_results_aza_noSyst %>% filter(metadata == 'aza_vs_noSyst_mod') 
maaslin2_results_aza_noSyst$qval <- p.adjust(maaslin2_results_aza_noSyst$pval, method = 'BH')

write.table(maaslin2_results_aza_noSyst, "Output_files/Maaslin2/Results/Remission_crp/aza_vs_noSyst_crp/maaslin2_results_aza_vs_noSyst_R_crp.txt", sep = "\t",  quote = FALSE)
