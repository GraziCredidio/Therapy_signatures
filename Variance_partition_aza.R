# EZE cohort: Therapy Signatures
  # Variance partition: azathioprine x no systemic therapies 
  # Figures 14A, 14B, 14C  (Master's thesis) 
  # Author: Graziella Credidio


rm(list = ls())

# Loading packages ----
if (!require("BiocManager", quietly = TRUE)) BiocManager::install("variancePartition")
library(variancePartition)
library(tidyverse)

folder <- "Output_files/Variance_partition/aza"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading files ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_aza <- read.table("Cleaned_tables/models/aza/EZECohort_coldata_aza.txt", sep = "\t")
vst_counts <- read.table("Cleaned_tables/models/aza/vst_counts_aza.txt", sep = "\t")
sig_genes <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_aza_vs_noSyst.txt", sep = "\t")

# Filtering counts to significant genes
vst_sig <- vst_counts[sig_genes$feature,]

# Variance Partition----
# Transforming columns as factor
coldata_aza <- coldata_aza %>% 
  mutate(across(c(diagnosis_class,age_group, 
                  sex, bmi_class, aza_vs_noSyst), as.factor))

# Formula
form <- ~  crp_log + leucocytes +  erythrocytes + thrombocytes + (1|aza_vs_noSyst) + (1|biologics) + (1|diagnosis_class) + (1|age_group) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_aza)

write.table(varPart, "Output_files/Variance_partition/aza/varPart_aza.txt", sep="\t", row.names=TRUE)

# Plots ----
varPart_ord <- sortCols(varPart)
rownames(varPart_ord) <- make.names(unique_ensg2gene[rownames(varPart_ord), ]$hgnc_symbol, unique = TRUE)
plotVarPart(varPart_ord) #saved as: "Output_files/Variance_partition/aza/plots/varpart_violin_aza"

varPart_ord_genes <- as.data.frame(varPart_ord) %>% 
  arrange(desc(aza_vs_noSyst))
plotPercentBars(varPart_ord_genes[1:15,]) + theme(legend.position="bottom") #saved as: "Output_files/Variance_partition/aza/plots/varpart_aza_top15"

# Correlations ----
formula_corr <- ~ crp_log + diagnosis_class + age_group + bmi_class + aza_vs_noSyst + sex + biologics + leucocytes + erythrocytes + thrombocytes

C = canCorPairs(formula_corr, coldata_aza)
plotCorrMatrix(C) #saved as: "Output_files/Variance_partition/aza/plots/correlation"

# Save variance partition output for genes with aza explaining > 25% of variability ----
top25percent_varPart <- varPart_ord_genes %>% 
  filter(aza_vs_noSyst >= 0.25)

write.table(top25percent_varPart, "Output_files/Variance_partition/aza/top25percent_varPart_aza.txt", sep="\t", row.names=TRUE)