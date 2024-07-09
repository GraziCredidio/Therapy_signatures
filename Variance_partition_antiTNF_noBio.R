# EZE cohort 
# Variance Partition: Anti TNF x no Bio
  # Date: 05.05


graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC


# Loading packages ----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("variancePartition")

library(variancePartition)
library(tidyverse)

# Loading files ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_R_noBio <- read.table("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_antiTNF_vs_noBio_diag_R_20.03.txt", sep = "\t")
vst_counts_R_noBio <- read.table("Output_files/DESeq2/Maaslin2/Remission/DESeq2_vst_counts_antiTNF_vs_noBiologics_R_29.03.txt", sep = "\t")
  
sig_genes <- read.table("Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_R_vst_q0.05.txt", sep = "\t")

vst_sig <- vst_counts_R_noBio[sig_genes$feature,]

# Variance Partition ----
coldata_R_noBio <- coldata_R_noBio %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, Prednisolon, antiTNF_vs_noBiologics), as.factor))

form <- ~  crp_log + (1|antiTNF_vs_noBiologics) + (1|Prednisolon) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_noBio)

write.table(varPart, "Output_files/VariancePartition/Remission/sig_genes/noBio/varPart_antiTNF_noBio_R.txt", sep="\t", row.names=TRUE)

# Plots
vp_plot <- sortCols(varPart) 
rownames(vp_plot) <-  make.names(unique_ensg2gene[rownames(vp_plot), ]$hgnc_symbol, unique = TRUE)
plotVarPart(vp_plot) 

vp_noBio <- as.data.frame(vp_plot) %>% 
  arrange(desc(antiTNF_vs_noBiologics)) # plotting with gene names
plotPercentBars(vp_noBio[1:30,])


vp <- sortCols(varPart)
vp_noBio <- as.data.frame(vp) %>% 
  arrange(desc(antiTNF_vs_noBiologics)) # saving with ensembl id 

write.table(vp_noBio, "Output_files/VariancePartition/Remission/sig_genes/noBio/varPart_bynoBio_sigGenes.txt", sep="\t", row.names=TRUE)


# Checking correlation between variables ----
formula_corr <- ~ crp_log + diagnosis_class + age_group2 + bmi_class + antiTNF_vs_noBiologics + sex + Prednisolon
  
C = canCorPairs(formula_corr, coldata_R_noBio)
plotCorrMatrix(C)

############################################################
# Including leuko, thrombo and erythrocytes in formula ----
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_noBio$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)

idx <- match(coldata_R_noBio$study_id, coldata_R_noBio$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_noBio_lab <- left_join(coldata_R_noBio, ord_redcap_patients, by = "study_id")
names(coldata_R_noBio_lab)

# Variance Partition 
coldata_R_noBio_lab <- coldata_R_noBio_lab %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, Prednisolon, antiTNF_vs_noBiologics), as.factor))


form <- ~ crp_log + leucocytes +  erythrocytes + thrombocytes +(1|antiTNF_vs_noBiologics) + (1|Prednisolon) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_noBio_lab)

write.table(varPart, "Output_files/VariancePartition/Remission/sig_genes/noBio/lab_values/varPart_noBio_labValues_sigGenes.txt", sep="\t", row.names=TRUE)
#varPart <- read.table("Output_files/VariancePartition/Remission/sig_genes/varPart_aza_labValues_sigGenes.txt", sep="\t")

# Plots
vp_plot_lab <- sortCols(varPart)
rownames(vp_plot_lab) <-  make.names(unique_ensg2gene[rownames(vp_plot_lab), ]$hgnc_symbol, unique = TRUE)
plotVarPart(vp_plot_lab) 

vp_noBio_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(antiTNF_vs_noBiologics))
plotPercentBars(vp_noBio_lab[1:30,])

vp_lab <- sortCols(varPart)
vp_noBio_lab <- as.data.frame(vp_lab) %>% 
  arrange(desc(antiTNF_vs_noBiologics))

write.table(vp_noBio_lab, "Output_files/VariancePartition/Remission/sig_genes/noBio/lab_values/varPart_bynoBio_sigGenes_labValues.txt", sep="\t", row.names=TRUE)






#################### Remitters + crp ######################
# Loading files ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_R_noBio <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_antiTNF_R_crp_04.05.txt", sep = "\t")
vst_counts_R_noBio <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_antiTNF_R_crp_04.05.txt", sep = "\t")
sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_R_crp.txt", sep = "\t")

vst_sig <- vst_counts_R_noBio[sig_genes$feature,]

# Including leuko, thrombo and erythrocytes in formula ----
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_noBio$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)

idx <- match(coldata_R_noBio$study_id, coldata_R_noBio$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_noBio_lab <- left_join(coldata_R_noBio, ord_redcap_patients, by = "study_id")
names(coldata_R_noBio_lab)

# Variance Partition 
coldata_R_noBio_lab <- coldata_R_noBio_lab %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, Prednisolon, antiTNF_vs_noBiologics), as.factor))


form <- ~ crp_log + leucocytes +  erythrocytes + thrombocytes +(1|antiTNF_vs_noBiologics) + (1|Prednisolon) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_noBio_lab)

write.table(varPart, "Output_files/VariancePartition/Remission_crp/noBio/varPart_noBio_labValues_sigGenes.txt", sep="\t", row.names=TRUE)
varPart <- read.table("Output_files/VariancePartition/Remission_crp/noBio/varPart_noBio_labValues_sigGenes.txt", sep="\t")


# Plots
vp_plot_lab <- sortCols(varPart)
rownames(vp_plot_lab) <-  make.names(unique_ensg2gene[rownames(vp_plot_lab), ]$hgnc_symbol, unique = TRUE)
plotVarPart(vp_plot_lab) 

vp_noBio_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(antiTNF_vs_noBiologics))
plotPercentBars(vp_noBio_lab[1:15,])

vp_bmi <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(bmi_class))
plotPercentBars(vp_noBio_lab[1:30,])

write.table(vp_noBio_lab, "Output_files/VariancePartition/Remission_crp/noBio/varPart_bynoBio_sigGenes_labValues.txt", sep="\t", row.names=TRUE)


# Correlation
formula_corr <- ~ crp_log + diagnosis_class + age_group2 + bmi_class + antiTNF_vs_noBiologics + sex + Prednisolon + leucocytes + erythrocytes + thrombocytes

C = canCorPairs(formula_corr, coldata_R_noBio_lab)
plotCorrMatrix(C)




################# Data integration: methylation ###################
rm(list = ls())

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_R_noBio <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_antiTNF_R_crp_04.05.txt", sep = "\t")
vst_counts_R_noBio <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_antiTNF_R_crp_04.05.txt", sep = "\t")
sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_R_crp.txt", sep = "\t")

cor_genes <- read.table("Output_files/Methylation/antiTNF/antiTNF_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = "\t")
cor_genes <- unique(cor_genes$Gene)

vst_sig <- vst_counts_R_noBio[cor_genes,]

# Including leuko, thrombo and erythrocytes in formula 
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_noBio$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)

idx <- match(coldata_R_noBio$study_id, coldata_R_noBio$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_noBio_lab <- left_join(coldata_R_noBio, ord_redcap_patients, by = "study_id")
names(coldata_R_noBio_lab)

# Variance Partition 
coldata_R_noBio_lab <- coldata_R_noBio_lab %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, Prednisolon, antiTNF_vs_noBiologics), as.factor))


form <- ~ crp_log + leucocytes +  erythrocytes + thrombocytes +(1|antiTNF_vs_noBiologics) + (1|Prednisolon) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_noBio_lab)


write.table(varPart, "Output_files/Methylation/antiTNF/varPart/varPart_antiTNF_corGenes.txt", sep="\t", row.names=TRUE)

# Plots
vp_plot_lab <- sortCols(varPart)
rownames(vp_plot_lab) <-  make.names(unique_ensg2gene[rownames(vp_plot_lab), ]$hgnc_symbol, unique = TRUE)
plotVarPart(vp_plot_lab) 

vp_noBio_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(antiTNF_vs_noBiologics))
plotPercentBars(vp_noBio_lab[1:30,])


# creating varpart > 25%
top25percent_varPart_lab <- vp_noBio_lab %>% 
  filter(antiTNF_vs_noBiologics >= 0.25)

