# EZE cohort: Therapy Signatures
  # Variance partition: prednisolone x no systemic therapies
  # DEGS with absLFC > 0.5
  # Figures 9A, 9B, 9C  (Master's thesis) 
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
if (!require("BiocManager", quietly = TRUE)) BiocManager::install("variancePartition")
library(variancePartition)
library(tidyverse)

folder <- "Output_files/Variance_partition/pred"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading files ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_pred <- read.table("Cleaned_tables/models/pred/EZECohort_coldata_pred.txt", sep = "\t")
vst_counts <- read.table("Cleaned_tables/models/pred/vst_counts_pred.txt", sep = "\t")
sig_genes <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_pred_vs_noSyst.txt", sep = "\t")

# Filtering counts to significant genes with absLFC > 0.5
sig_genes <- sig_genes %>% 
  filter(abs(coef) > 0.5)

vst_sig <- vst_counts[sig_genes$feature,]

# Variance Partition ----
# Transforming columns as factor
coldata_pred <- coldata_pred %>% 
  mutate(across(c(diagnosis_class, age_group,  
                  sex, bmi_class, Prednisolon, pred_vs_noSyst), as.factor))

# Formula
form <- ~  crp_log + leucocytes +  erythrocytes + thrombocytes + (1|pred_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_pred)

write.table(varPart, "Output_files/Variance_partition/pred/varPart_pred.txt", sep="\t", row.names=TRUE)

# Plots ----
varpart_ord <- sortCols(varPart) 
rownames(varpart_ord) <- make.names(unique_ensg2gene[rownames(varpart_ord), ]$hgnc_symbol, unique = TRUE)
plotVarPart(varpart_ord) #saved as: "Output_files/Variance_partition/pred/plots/varpart_violin_pred"

varpart_ord_genes <- as.data.frame(varpart_ord) %>% 
  arrange(desc(pred_vs_noSyst)) 
plotPercentBars(varpart_ord_genes[1:15,]) #saved as: "Output_files/Variance_partition/pred/plots/varpart_pred_top15"

# Correlations ----
formula_corr <- ~ crp_log + diagnosis_class + age_group + bmi_class + pred_vs_noSyst + sex + biologics + leucocytes + erythrocytes + thrombocytes

C = canCorPairs(formula_corr, coldata_pred)
plotCorrMatrix(C) #saved as: "Output_files/Variance_partition/pred/plots/correlation"
