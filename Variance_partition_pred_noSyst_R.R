# EZE cohort 
# Variance Partition: Pred x no Syst (remitters)
  # Date: 05.05


graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC


# Loading packages ----
if (!require("BiocManager", quietly = TRUE))
  BiocManager::install("variancePartition")
library('variancePartition')
library(tidyverse)

# Loading files ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_R_pred <- read.table("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_pred_vs_noSyst_R_17.03.txt", sep = "\t")
vst_counts_R_pred <- read.table("Output_files/DESeq2/Maaslin2/Remission/DESeq2_vst_counts_pred_vs_noSyst_R_29.03.txt", sep = "\t")

sig_genes <- read.table( "Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_vst_q0.05.txt", sep = "\t")

vst_sig <- vst_counts_R_pred[sig_genes$feature,]

# Variance Partition ----
coldata_R_pred <- coldata_R_pred %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, Prednisolon, pred_vs_noSyst), as.factor))

form <- ~  crp_log + (1|pred_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
#varPart <- fitExtractVarPartModel(vst_counts_R_pred, form, coldata_R_pred) #all genes

varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_pred) #sig genes
write.table(varPart, "Output_files/VariancePartition/Remission/sig_genes/pred/varPart_pred_noSyst.txt", sep="\t", row.names=TRUE)

# Plots ----
vp_plot <- sortCols(varPart) 
rownames(vp_plot) <- make.names(unique_ensg2gene[rownames(vp_plot), ]$hgnc_symbol, unique = TRUE)
plotVarPart(vp_plot) 

vp_pred <- as.data.frame(vp_plot) %>% 
  arrange(desc(pred_vs_noSyst)) # plotting with gene names
plotPercentBars(vp_pred[1:30,])


vp <- sortCols(varPart)
vp_pred <- as.data.frame(vp) %>% 
  arrange(desc(pred_vs_noSyst)) # saving with ensembl id 
write.table(vp_pred, "Output_files/VariancePartition/Remission/sig_genes/pred/varPart_byPred_sigGenes.txt", sep="\t", row.names=TRUE)

# Checking correlation between variables ----
formula_corr <- ~ crp_log + diagnosis_class + age_group2 + bmi_class + pred_vs_noSyst + sex + biologics

C = canCorPairs( formula_corr, coldata_R_pred)
plotCorrMatrix( C )


# Checking if genes that have most variance explained by pred x no Syst are also present in intersections
intersec_genes <- readRDS("Output_files/Maaslin2/significant_results/list_sig_genes_all_remitters_intersect_diff_06.04.rds")
int_R_pred_topvar <- intersec_genes$intersect_R_pred_topVar

top25percent_varPart <- vp_pred %>% 
  filter(pred_vs_noSyst >= 0.25)
top25percent_varPart_genes <- rownames(top25percent_varPart) #37 genes

topvarPart_intersectVenn <- as.data.frame(intersect(top25percent_varPart_genes, int_R_pred_topvar)) # 23 genes
colnames(topvarPart_intersectVenn) <- "Genes"

# adding gene names
topvarPart_intersectVenn$Gene_names <-  make.names(unique_ensg2gene[topvarPart_intersectVenn$Genes, ]$hgnc_symbol, unique = TRUE)
write.table(topvarPart_intersectVenn, "Output_files/VariancePartition/Remission/sig_genes/pred/genes_top25percent_varPart_intersectVenn_pred_R.txt", sep="\t", row.names=TRUE)




############################################################
rm(list = ls())

# Including leuko, thrombo and erythrocytes in formula ----
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_pred$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)

idx <- match(coldata_R_pred$study_id, redcap_patients$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_pred_lab <- left_join(coldata_R_pred, ord_redcap_patients, by = "study_id")

# Variance Partition 
coldata_R_pred_lab <- coldata_R_pred_lab %>% 
  mutate(across(c(diagnosis_class,age_group2, 
                  sex, bmi_class, pred_vs_noSyst), as.factor))

form <- ~  crp_log + leucocytes +  erythrocytes + thrombocytes + (1|pred_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_pred_lab)

write.table(varPart, "Output_files/VariancePartition/Remission/sig_genes/pred/lab_values/varPart_pred_labValues_sigGenes.txt", sep="\t", row.names=TRUE)
varPart <- read.table( "Output_files/VariancePartition/Remission/sig_genes/pred/lab_values/varPart_pred_labValues_sigGenes.txt", sep="\t")

# Plots
vp_plot_lab <- sortCols(varPart)
rownames(vp_plot_lab) <- make.names(unique_ensg2gene[rownames(vp_plot_lab), ]$hgnc_symbol, unique = TRUE)
plotVarPart(vp_plot_lab) 

vp_pred_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(pred_vs_noSyst))
plotPercentBars(vp_pred_lab[1:30,])


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
vp_pred_lab <- as.data.frame(vp_lab) %>% 
  arrange(desc(pred_vs_noSyst))

write.table(vp_pred_lab, "Output_files/VariancePartition/Remission/sig_genes/pred/varPart_byPred_sigGenes_labValues.txt", sep="\t", row.names=TRUE)

# Checking correlation between variables ----
formula_corr <- ~ crp_log + diagnosis_class + age_group2 + bmi_class + pred_vs_noSyst + sex + biologics + leucocytes + erythrocytes + thrombocytes

C = canCorPairs(formula_corr, coldata_R_pred_lab)
plotCorrMatrix( C )

# Checking if genes that have most variance explained by aza x no Syst with lab values are also present in intersections ----
intersec_genes <- readRDS("Output_files/Maaslin2/significant_results/list_sig_genes_all_remitters_intersect_diff_06.04.rds")
int_R_pred_topvar <- intersec_genes$intersect_R_pred_topVar

top25percent_varPart_lab <- vp_pred_lab %>% 
  filter(pred_vs_noSyst >= 0.25)
top25percent_varPart_lab_genes <- rownames(top25percent_varPart_lab) # 13 genes

topvarPart_intersectVenn_lab <- as.data.frame(intersect(top25percent_varPart_lab_genes, int_R_pred_topvar)) # 70 genes
colnames(topvarPart_intersectVenn_lab) <- "Genes"

# adding gene names
topvarPart_intersectVenn_lab$Gene_names <-  make.names(unique_ensg2gene[topvarPart_intersectVenn_lab$Genes, ]$hgnc_symbol, unique = TRUE)
write.table(topvarPart_intersectVenn_lab, "Output_files/VariancePartition/Remission/sig_genes/pred/genes_top25percent_varPart_intersectVenn_pred_R_lab.txt", sep="\t", row.names=TRUE)










################# Remitters + crp ###################
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_R_pred <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_pred_R_crp_04.05.txt", sep = "\t")
vst_counts_R_pred <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_pred_R_crp_04.05.txt", sep = "\t")
sig_genes <- read.table( "Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt", sep = "\t")

vst_sig <- vst_counts_R_pred[sig_genes$feature,]

# Including leuko, thrombo and erythrocytes in formula ----
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_pred$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)

idx <- match(coldata_R_pred$study_id, redcap_patients$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_pred_lab <- left_join(coldata_R_pred, ord_redcap_patients, by = "study_id")

# Variance Partition 
coldata_R_pred_lab <- coldata_R_pred_lab %>% 
  mutate(across(c(diagnosis_class,age_group2, 
                  sex, bmi_class, pred_vs_noSyst), as.factor))

form <- ~  crp_log + leucocytes +  erythrocytes + thrombocytes + (1|pred_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_pred_lab)

write.table(varPart, "Output_files/VariancePartition/Remission_crp/pred/varPart_pred_labValues_sigGenes.txt", sep="\t", row.names=TRUE)
varPart <- read.table( "Output_files/VariancePartition/Remission_crp/pred/varPart_pred_labValues_sigGenes.txt", sep="\t")

# Plots
vp_plot_lab <- sortCols(varPart)
rownames(vp_plot_lab) <- make.names(unique_ensg2gene[rownames(vp_plot_lab), ]$hgnc_symbol, unique = TRUE)
plotVarPart(vp_plot_lab) 

vp_pred_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(pred_vs_noSyst))
plotPercentBars(vp_pred_lab[1:30,])


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
vp_pred_lab <- as.data.frame(vp_lab) %>% 
  arrange(desc(pred_vs_noSyst))

write.table(vp_pred_lab, "Output_files/VariancePartition/Remission_crp/pred/varPart_byPred_sigGenes_labValues.txt", sep="\t", row.names=TRUE)

# Checking correlation between variables ----
formula_corr <- ~ crp_log + diagnosis_class + age_group2 + bmi_class + pred_vs_noSyst + sex + biologics + leucocytes + erythrocytes + thrombocytes

C = canCorPairs(formula_corr, coldata_R_pred_lab)
plotCorrMatrix( C )

# creating varpart > 25%
top25percent_varPart_lab <- vp_pred_lab %>% 
  filter(pred_vs_noSyst >= 0.25)

write.table(top25percent_varPart_lab, "Output_files/VariancePartition/Remission_crp/pred/top25percent_varPart_pred_R_crp.txt", sep="\t", row.names=TRUE)






################# Remitters + crp (L2FC 0.5) ###################
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_R_pred <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_pred_R_crp_04.05.txt", sep = "\t")
vst_counts_R_pred <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_pred_R_crp_04.05.txt", sep = "\t")
sig_genes <- read.table( "Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt", sep = "\t")

sig_genes <- sig_genes %>% 
  filter(abs(coef) > 0.5)

vst_sig <- vst_counts_R_pred[sig_genes$feature,]

# Including leuko, thrombo and erythrocytes in formula ----
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_pred$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes, neutrophils)

idx <- match(coldata_R_pred$study_id, redcap_patients$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_pred_lab <- left_join(coldata_R_pred, ord_redcap_patients, by = "study_id")

# Variance Partition 
coldata_R_pred_lab <- coldata_R_pred_lab %>% 
  mutate(across(c(diagnosis_class,age_group2, 
                  sex, bmi_class, pred_vs_noSyst), as.factor))

form <- ~  crp_log + leucocytes +  erythrocytes + thrombocytes + (1|pred_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_pred_lab)


write.table(varPart, "Output_files/VariancePartition/Remission_crp/pred/lfc0.5/varPart_pred_labValues_sigGenes_lfc0.5.txt", sep="\t", row.names=TRUE)

# Plots
vp_plot_lab <- sortCols(varPart)
rownames(vp_plot_lab) <- make.names(unique_ensg2gene[rownames(vp_plot_lab), ]$hgnc_symbol, unique = TRUE)
plotVarPart(vp_plot_lab) 

vp_pred_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(pred_vs_noSyst))
plotPercentBars(vp_pred_lab[1:30,])

vp_lab <- sortCols(varPart)
vp_pred_lab <- as.data.frame(vp_lab) %>% 
  arrange(desc(pred_vs_noSyst))


################# Data integration: methylation ###################
rm(list = ls())

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id


coldata_R_pred <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_pred_R_crp_04.05.txt", sep = "\t")
vst_counts_R_pred <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_pred_R_crp_04.05.txt", sep = "\t")
sig_genes <- read.table( "Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt", sep = "\t")

cor_genes <- read.table("Output_files/Methylation/pred/pred_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = "\t")
cor_genes <- unique(cor_genes$Gene)

vst_sig <- vst_counts_R_pred[cor_genes,]

# Including leuko, thrombo and erythrocytes in formula 
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_pred$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)

idx <- match(coldata_R_pred$study_id, redcap_patients$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_pred_lab <- left_join(coldata_R_pred, ord_redcap_patients, by = "study_id")

# Variance Partition 
coldata_R_pred_lab <- coldata_R_pred_lab %>% 
  mutate(across(c(diagnosis_class,age_group2, 
                  sex, bmi_class, pred_vs_noSyst), as.factor))

form <- ~  crp_log + leucocytes +  erythrocytes + thrombocytes + (1|pred_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_pred_lab)


write.table(varPart, "Output_files/Methylation/pred/varPart/varPart_pred_corGenes.txt", sep="\t", row.names=TRUE)


varPart <- read.table("Output_files/Methylation/pred/varPart/varPart_pred_corGenes.txt", sep="\t")

# Plots
vp_plot_lab <- sortCols(varPart)
rownames(vp_plot_lab) <-  make.names(unique_ensg2gene[rownames(vp_plot_lab), ]$hgnc_symbol, unique = TRUE)
plotVarPart(vp_plot_lab) 

vp_pred_lab <- as.data.frame(vp_plot_lab) %>% 
  arrange(desc(pred_vs_noSyst))
plotPercentBars(vp_pred_lab[1:30,])


# creating varpart > 25%
top25percent_varPart_lab <- vp_pred_lab %>% 
  filter(pred_vs_noSyst >= 0.25)

write.table(top25percent_varPart_lab, "Output_files/Methylation/pred/varPart/varPart_pred_corGenes_25percent.txt", sep="\t", row.names=TRUE) #99 genes







################# Data integration: methylation + topVar ###################
rm(list = ls())

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVarGenes_names <- rownames(topVarGenes)

coldata_R_pred <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_pred_R_crp_04.05.txt", sep = "\t")
vst_counts_R_pred <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_pred_R_crp_04.05.txt", sep = "\t")
sig_genes <- read.table( "Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt", sep = "\t")

cor_genes <- read.table("Output_files/Methylation/pred/pred_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = "\t")
cor_genes <- unique(cor_genes$Gene)

cor_genes_topVar <- intersect(cor_genes, topVarGenes_names)
vst_sig <- vst_counts_R_pred[cor_genes_topVar,]

# Including leuko, thrombo and erythrocytes in formula 
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_pred$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)

idx <- match(coldata_R_pred$study_id, redcap_patients$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_pred_lab <- left_join(coldata_R_pred, ord_redcap_patients, by = "study_id")

# Variance Partition 
coldata_R_pred_lab <- coldata_R_pred_lab %>% 
  mutate(across(c(diagnosis_class,age_group2, 
                  sex, bmi_class, pred_vs_noSyst), as.factor))

form <- ~  crp_log + leucocytes +  erythrocytes + thrombocytes + (1|pred_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group2) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_R_pred_lab)

write.table(varPart, "Output_files/Methylation/pred/varPart/corGenes + topVar/varPart_pred_corGenes_topVar.txt", sep="\t", row.names=TRUE)

# Plots
vp_plot_lab <- sortCols(varPart)
rownames(vp_plot_lab) <-  make.names(unique_ensg2gene[rownames(vp_plot_lab), ]$hgnc_symbol, unique = TRUE)
plotVarPart(vp_plot_lab) 

vp_pred_lab <- as.data.frame(vp_plot_lab) %>% 
  dplyr::arrange(desc(pred_vs_noSyst))
plotPercentBars(vp_pred_lab[1:30,])


# creating varpart > 25%
top25percent_varPart_lab <- vp_pred_lab %>% 
  filter(pred_vs_noSyst >= 0.25)

write.table(top25percent_varPart_lab, "Output_files/Methylation/pred/varPart/corGenes + topVar/varPart_pred_corGenes_toVar_25percent.txt", sep="\t", row.names=TRUE) #99 genes

