# EZE cohort 
# Variance Partition: Aza x no Syst (remitters) - significant genes
# Date: 05.05


graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC


# Loading packages ----
if (!require("BiocManager", quietly = TRUE)) BiocManager::install("variancePartition")
library('variancePartition')
library(tidyverse)

# Loading files ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_R_aza <- read.table("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_aza_vs_noSyst_R_17.03.txt", sep = "\t")
vst_counts_R_aza <- read.table("Output_files/DESeq2/Maaslin2/Remission/DESeq2_vst_counts_aza_vs_noSyst_R_29.03.txt", sep = "\t")

sig_genes <- read.table( "Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_vst_q0.05.txt", sep = "\t")

vst_sig <- vst_counts_R_aza[sig_genes$feature,]

# Variance Partition ----
coldata_R_aza <- coldata_R_aza %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, aza_vs_noSyst), as.factor))

form <- ~  crp_log + (1|aza_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)

#varPart <- fitExtractVarPartModel(vst_counts_R_aza, form, coldata_R_aza) #all genes
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_aza) #significant genes
write.table(varPart, "Output_files/VariancePartition/Remission/sig_genes/aza/varPart_aza_noSyst_sigGenes.txt", sep="\t", row.names=TRUE)

# Plots ----
vp_plot <- sortCols(varPart) #sort variables (i.e. columns) by median fraction of variance explained
rownames(vp_plot) <-  make.names(unique_ensg2gene[rownames(vp_plot), ]$hgnc_symbol, unique = TRUE)
plotPercentBars(vp_plot[1:30,])   # Bar plot of variance fractions for the first 10 genes. 
plotVarPart(vp_plot) # violin plot of contribution of each variable to total variance

vp_aza <- as.data.frame(vp_plot) %>% 
  arrange(desc(aza_vs_noSyst)) # plotting with gene names
plotPercentBars(vp_aza[1:30,])


vp <- sortCols(varPart)
vp_aza <- as.data.frame(vp) %>% 
  arrange(desc(aza_vs_noSyst)) # saving with ensembl id 
write.table(vp_aza, "Output_files/VariancePartition/Remission/sig_genes/aza/varPart_byAza_sigGenes.txt", sep="\t", row.names=TRUE)

# vp_res <- as.data.frame(vp) %>% 
#   arrange(desc(Residuals))
# plotPercentBars(vp_res[1:10,])
# 
# vp_crp <- as.data.frame(vp) %>% 
#   arrange(desc(crp_log))
# plotPercentBars(vp_crp[1:10,])

# Checking correlation between variables ----
formula_corr <- ~ crp_log + diagnosis_class + age_group2 + bmi_class + aza_vs_noSyst + sex + biologics

C = canCorPairs(formula_corr, coldata_R_aza)
plotCorrMatrix( C )

# Checking if genes that have most variance explained by aza x no Syst are also present in intersections ----
intersec_genes <- readRDS("Output_files/Maaslin2/significant_results/list_sig_genes_all_remitters_intersect_diff_06.04.rds")
int_R_aza_topvar <- intersec_genes$intersect_R_aza_topVar

top25percent_varPart <- vp_aza %>% 
  filter(aza_vs_noSyst >= 0.25)
top25percent_varPart_genes <- rownames(top25percent_varPart) #201 genes

topvarPart_intersectVenn <- as.data.frame(intersect(top25percent_varPart_genes, int_R_aza_topvar)) # 93 genes
colnames(topvarPart_intersectVenn) <- "Genes"

# adding gene names
topvarPart_intersectVenn$Gene_names <-  make.names(unique_ensg2gene[topvarPart_intersectVenn$Genes, ]$hgnc_symbol, unique = TRUE)
write.table(topvarPart_intersectVenn, "Output_files/VariancePartition/Remission/sig_genes/aza/genes_top25percent_varPart_intersectVenn_aza_R.txt", sep="\t", row.names=TRUE)



##### Including leuko, thrombo and erythrocytes in formula ----
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_aza$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)

idx <- match(coldata_R_aza$study_id, redcap_patients$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_aza_lab <- left_join(coldata_R_aza, ord_redcap_patients, by = "study_id")
names(coldata_R_aza_lab)

# Variance Partition 
coldata_R_aza_lab <- coldata_R_aza_lab %>% 
  mutate(across(c(diagnosis_class,age_group2, 
                  sex, bmi_class, aza_vs_noSyst), as.factor))

form <- ~  crp_log + leucocytes +  erythrocytes + thrombocytes + (1|aza_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_aza_lab)
write.table(varPart, "Output_files/VariancePartition/Remission/sig_genes/aza/varPart_aza_labValues_sigGenes.txt", sep="\t", row.names=TRUE)
varPart <- read.table("Output_files/VariancePartition/Remission/sig_genes/aza/varPart_aza_labValues_sigGenes.txt", sep="\t")

# Plots
vp_plot_lab <- sortCols(varPart)
rownames(vp_plot_lab) <-  make.names(unique_ensg2gene[rownames(vp_plot_lab), ]$hgnc_symbol, unique = TRUE)
plotPercentBars(vp_plot_lab[1:30,]) 
plotVarPart(vp_plot_lab) 

vp_aza_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(aza_vs_noSyst))
plotPercentBars(vp_aza[1:30,])


vp_sex_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(sex))
plotPercentBars(vp_sex_lab[1:30,])

vp_crp_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(crp_log))
plotPercentBars(vp_crp_lab[1:30,])

vp_residuals_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(Residuals))
plotPercentBars(vp_residuals_lab[1:30,])


vp_lab <- sortCols(varPart)
vp_aza_lab <- as.data.frame(vp_lab) %>% 
  arrange(desc(aza_vs_noSyst))

write.table(vp_aza_lab, "Output_files/VariancePartition/Remission/sig_genes/aza/lab_values/varPart_byAza_sigGenes_labValues.txt", sep="\t", row.names=TRUE)

# Checking correlation between variables ----
formula_corr <- ~ crp_log + diagnosis_class + age_group2 + bmi_class + aza_vs_noSyst + sex + biologics + leucocytes + erythrocytes + thrombocytes

C = canCorPairs(formula_corr, coldata_R_aza_lab)
plotCorrMatrix( C )

# Checking if genes that have most variance explained by aza x no Syst with lab values are also present in intersections ----
intersec_genes <- readRDS("Output_files/Maaslin2/significant_results/list_sig_genes_all_remitters_intersect_diff_06.04.rds")
int_R_aza_topvar <- intersec_genes$intersect_R_aza_topVar

top25percent_varPart_lab <- vp_aza_lab %>% 
  filter(aza_vs_noSyst >= 0.25)
top25percent_varPart_lab_genes <- rownames(top25percent_varPart_lab) #151 genes

topvarPart_intersectVenn_lab <- as.data.frame(intersect(top25percent_varPart_lab_genes, int_R_aza_topvar)) # 70 genes
colnames(topvarPart_intersectVenn_lab) <- "Genes"

# adding gene names
topvarPart_intersectVenn_lab$Gene_names <-  make.names(unique_ensg2gene[topvarPart_intersectVenn_lab$Genes, ]$hgnc_symbol, unique = TRUE)
write.table(topvarPart_intersectVenn_lab, "Output_files/VariancePartition/Remission/sig_genes/aza/genes_top25percent_varPart_intersectVenn_aza_R_lab.txt", sep="\t", row.names=TRUE)






################# Remitters + crp ###################
rm(list = ls())

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_R_aza <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_aza_R_crp_04.05.txt", sep = "\t")
vst_counts_R_aza <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_aza_R_crp_04.05.txt", sep = "\t")
sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_crp.txt", sep = "\t")

vst_sig <- vst_counts_R_aza[sig_genes$feature,]

# Including leuko, thrombo and erythrocytes in formula 
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_aza$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)

idx <- match(coldata_R_aza$study_id, redcap_patients$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_aza_lab <- left_join(coldata_R_aza, ord_redcap_patients, by = "study_id")
names(coldata_R_aza_lab)

# Variance Partition 
coldata_R_aza_lab <- coldata_R_aza_lab %>% 
  mutate(across(c(diagnosis_class,age_group2, 
                  sex, bmi_class, aza_vs_noSyst), as.factor))

form <- ~  crp_log + leucocytes +  erythrocytes + thrombocytes + (1|aza_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_aza_lab)

write.table(varPart, "Output_files/VariancePartition/Remission_crp/aza/varPart_aza_labValues_sigGenes.txt", sep="\t", row.names=TRUE)
#varPart <- read.table("Output_files/VariancePartition/Remission_crp/aza/varPart_aza_labValues_sigGenes.txt", sep="\t")
# Plots
vp_plot_lab <- sortCols(varPart)
rownames(vp_plot_lab) <-  make.names(unique_ensg2gene[rownames(vp_plot_lab), ]$hgnc_symbol, unique = TRUE)
plotVarPart(vp_plot_lab) 

vp_aza_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(aza_vs_noSyst))
plotPercentBars(vp_aza_lab[1:15,]) + theme(legend.position="bottom")


vp_sex_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(sex))
plotPercentBars(vp_sex_lab[1:30,])

vp_crp_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(crp_log))
plotPercentBars(vp_crp_lab[1:30,])

vp_age_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(age_group2))
plotPercentBars(vp_age_lab[1:30,])


vp_lab <- sortCols(varPart)
vp_aza_lab <- as.data.frame(vp_lab) %>% 
  arrange(desc(aza_vs_noSyst))

write.table(vp_aza_lab, "Output_files/VariancePartition/Remission_crp/aza/varPart_byAza_sigGenes_labValues.txt", sep="\t", row.names=TRUE)

# Checking correlation between variables ----
formula_corr <- ~ crp_log + diagnosis_class + age_group2 + bmi_class + aza_vs_noSyst + sex + biologics + leucocytes + erythrocytes + thrombocytes

C = canCorPairs(formula_corr, coldata_R_aza_lab)
plotCorrMatrix( C )



# creating varpart > 25%
top25percent_varPart_lab <- vp_aza_lab %>% 
  filter(aza_vs_noSyst >= 0.25)

write.table(top25percent_varPart_lab, "Output_files/VariancePartition/Remission_crp/aza/top25percent_varPart_aza_R_crp.txt", sep="\t", row.names=TRUE)



# 
# ################# Data integration: methylation ###################
# rm(list = ls())
# 
# ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
# unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
# rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id
# 
# coldata_R_aza <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_aza_R_crp_04.05.txt", sep = "\t")
# vst_counts_R_aza <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_aza_R_crp_04.05.txt", sep = "\t")
# sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_crp.txt", sep = "\t")
# cor_genes <- read.table("Output_files/Methylation/aza/aza_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = "\t")
# cor_genes <- unique(cor_genes$Gene)
# 
# vst_sig <- vst_counts_R_aza[cor_genes,]
# 
# # Including leuko, thrombo and erythrocytes in formula 
# redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")
# 
# # Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
# redcap_patients <- redcap[redcap$study_id %in% coldata_R_aza$study_id,]
# redcap_patients <- redcap_patients %>% 
#   dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)
# 
# idx <- match(coldata_R_aza$study_id, redcap_patients$study_id) #matching order of patients
# ord_redcap_patients <- redcap_patients[idx,]
# 
# coldata_R_aza_lab <- left_join(coldata_R_aza, ord_redcap_patients, by = "study_id")
# names(coldata_R_aza_lab)
# 
# # Variance Partition 
# coldata_R_aza_lab <- coldata_R_aza_lab %>% 
#   mutate(across(c(diagnosis_class,age_group2, 
#                   sex, bmi_class, aza_vs_noSyst), as.factor))
# 
# form <- ~  crp_log + leucocytes +  erythrocytes + thrombocytes + (1|aza_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
# varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_aza_lab)
# 
# write.table(varPart, "Output_files/Methylation/aza/varPart/varPart_aza_corGenes.txt", sep="\t", row.names=TRUE)
# 
# # Plots
# vp_plot_lab <- sortCols(varPart)
# rownames(vp_plot_lab) <-  make.names(unique_ensg2gene[rownames(vp_plot_lab), ]$hgnc_symbol, unique = TRUE)
# plotVarPart(vp_plot_lab) 
# 
# vp_aza_lab <- as.data.frame(vp_plot_lab) %>% 
#   arrange(desc(aza_vs_noSyst))
# plotPercentBars(vp_aza_lab[1:30,])
# 
# vp_lab <- sortCols(varPart)
# vp_aza_lab <- as.data.frame(vp_lab) %>% 
#   arrange(desc(aza_vs_noSyst))
# 
# # creating varpart > 25%
# top25percent_varPart_lab <- vp_aza_lab %>% 
#   filter(aza_vs_noSyst >= 0.25)
# 
# write.table(top25percent_varPart_lab, "Output_files/Methylation/aza/varPart/varPart_aza_corGenes_25percent.txt", sep="\t", row.names=TRUE) #99 genes
# 
# 
# 
# 
# ################# Data integration: methylation + topVar ###################
# rm(list = ls())
# 
# ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
# unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
# rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id
# 
# topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
# topVarGenes_names <- rownames(topVarGenes)
# 
# coldata_R_aza <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_aza_R_crp_04.05.txt", sep = "\t")
# vst_counts_R_aza <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_aza_R_crp_04.05.txt", sep = "\t")
# sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_crp.txt", sep = "\t")
# cor_genes <- read.table("Output_files/Methylation/aza/aza_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = "\t")
# cor_genes <- unique(cor_genes$Gene)
# 
# cor_genes_topVar <- intersect(cor_genes, topVarGenes_names)
# vst_sig <- vst_counts_R_aza[cor_genes_topVar,]
# 
# # Including leuko, thrombo and erythrocytes in formula 
# redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")
# 
# # Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
# redcap_patients <- redcap[redcap$study_id %in% coldata_R_aza$study_id,]
# redcap_patients <- redcap_patients %>% 
#   dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)
# 
# idx <- match(coldata_R_aza$study_id, redcap_patients$study_id) #matching order of patients
# ord_redcap_patients <- redcap_patients[idx,]
# 
# coldata_R_aza_lab <- left_join(coldata_R_aza, ord_redcap_patients, by = "study_id")
# names(coldata_R_aza_lab)
# 
# # Variance Partition 
# coldata_R_aza_lab <- coldata_R_aza_lab %>% 
#   mutate(across(c(diagnosis_class,age_group2, 
#                   sex, bmi_class, aza_vs_noSyst), as.factor))
# 
# form <- ~  crp_log + leucocytes +  erythrocytes + thrombocytes + (1|aza_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
# varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_aza_lab)
# 
# write.table(varPart, "Output_files/Methylation/aza/varPart/corGenes + topVar/varPart_aza_corGenes_topVar.txt", sep="\t", row.names=TRUE)
# 
# # Plots
# vp_plot_lab <- sortCols(varPart)
# rownames(vp_plot_lab) <-  make.names(unique_ensg2gene[rownames(vp_plot_lab), ]$hgnc_symbol, unique = TRUE)
# plotVarPart(vp_plot_lab) 
# 
# vp_aza_lab <- as.data.frame(vp_plot_lab)%>% 
# dplyr::arrange(desc(aza_vs_noSyst))
# plotPercentBars(vp_aza_lab[1:30,])
# 
# vp_lab <- sortCols(varPart)
# vp_aza_lab <- as.data.frame(vp_lab) %>% 
#   arrange(desc(aza_vs_noSyst))
# 
# # creating varpart > 25%
# top25percent_varPart_lab <- vp_aza_lab %>% 
#   filter(aza_vs_noSyst >= 0.25)
# 
# write.table(top25percent_varPart_lab, "Output_files/Methylation/aza/varPart/corGenes + topVar/varPart_aza_corGenes_topVar_25percent.txt", sep="\t", row.names=TRUE) #99 genes
