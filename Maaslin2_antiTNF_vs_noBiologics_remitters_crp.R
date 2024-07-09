# Modeling EZE cohort
# MaAslin2 (gene expression model): Anti-TNF x No biologics
# Inactive disease + crp < 5

#graphics.off()
#rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("D:\\Documentos/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

library(tidyverse)
library(Maaslin2)


coldata_noBio <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_antiTNF_R_crp_04.05.txt", sep = "\t")
vst_noBio <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_antiTNF_R_crp_04.05.txt", sep = "\t")

# setting no_biologics as the reference value
unique(coldata_noBio$antiTNF_vs_noBiologics)
coldata_noBio$antiTNF_vs_noBiologics_mod = coldata_noBio$antiTNF_vs_noBiologics
coldata_noBio$antiTNF_vs_noBiologics_mod[coldata_noBio$antiTNF_vs_noBiologics_mod == "no_biologics"] =
  "a_noBiologics"
coldata_noBio$antiTNF_vs_noBiologics_mod[coldata_noBio$antiTNF_vs_noBiologics_mod == "anti_tnf"] =
  "z_antiTnf"

coldata_noBio <- coldata_noBio %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, Prednisolon, antiTNF_vs_noBiologics_mod), as.factor))


dir.create("Output_files/Maaslin2/Results/Remission_crp/antiTNF_vs_noBiologics_crp")

fit_data_noBio = Maaslin2( 
  input_data = vst_noBio, 
  input_metadata = coldata_noBio, 
  analysis_method = "LM",
  normalization = "NONE",
  transform = "NONE",
  output = "Output_files/Maaslin2/Results/Remission_crp/antiTNF_vs_noBiologics_crp", 
  fixed_effects = c("antiTNF_vs_noBiologics_mod", "crp_log", "sex", "Prednisolon"),
  random_effects = c("diagnosis_class", "age_group2", "bmi_class"),
  reference = c("antiTNF_vs_noBiologics,no_biologics"))

maaslin2_all_results_noBio <- fit_data_noBio$results
maaslin2_results_noBio <- maaslin2_all_results_noBio %>% filter(metadata == 'antiTNF_vs_noBiologics_mod')
maaslin2_results_noBio$qval <- p.adjust(maaslin2_results_noBio$pval, method = 'BH')

write.table(maaslin2_results_noBio, "Output_files/Maaslin2/Results/Remission_crp/antiTNF_vs_noBiologics_crp/maaslin2_results_antiTNF_vs_noBio_R_crp.txt", sep = "\t",  quote = FALSE)

