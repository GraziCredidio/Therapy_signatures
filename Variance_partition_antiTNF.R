# EZE cohort: Therapy Signatures
  # Variance partition: antiTNF x no biologics
  # Figures 4A, 4B, 4C  (Master's thesis) 
  # Author: Graziella Credidio

rm(list = ls())

folder <- "Output_files/Variance_partition/antiTNF"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading packages ----
library(variancePartition)
library(tidyverse)

# Loading files ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_antiTNF <- read.table("Cleaned_tables/models/antiTNF/EZECohort_coldata_antiTNF.txt", sep = "\t")
vst_counts <- read.table("Cleaned_tables/models/antiTNF/vst_counts_antiTNF.txt", sep = "\t")
sig_genes <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_antiTNF_vs_noBio.txt", sep = "\t")

# Filtering counts to significant genes
vst_sig <- vst_counts[sig_genes$feature,]

# Variance Partition ----
# Transforming columns as factor
coldata_antiTNF <- coldata_antiTNF %>% 
  mutate(across(c(diagnosis_class, age_group,
                  sex, bmi_class, Prednisolon, antiTNF_vs_noBiologics), as.factor))


form <- ~ crp_log + leucocytes +  erythrocytes + thrombocytes +(1|antiTNF_vs_noBiologics) + (1|Prednisolon) + 
  (1|diagnosis_class) + (1|age_group) + (1|bmi_class)  + (1|sex)
varPart <- fitExtractVarPartModel(vst_sig, form, coldata_antiTNF)

write.table(varPart, "Output_files/Variance_partition/antiTNF/varPart_antiTNF.txt", sep="\t", row.names=TRUE)

# Plots ----
varpart_ord <- sortCols(varPart) 
rownames(varpart_ord) <- make.names(unique_ensg2gene[rownames(varpart_ord), ]$hgnc_symbol, unique = TRUE)
plotVarPart(varpart_ord) #saved as: "Output_files/Variance_partition/antiTNF/plots/varpart_violin_antiTNF"

varpart_ord_genes <- as.data.frame(varpart_ord) %>% 
  arrange(desc(antiTNF_vs_noBiologics)) 
plotPercentBars(varpart_ord_genes[1:15,]) #saved as: "Output_files/Variance_partition/antiTNF/plots/varpart_antiTNF_top15"

# Correlations ----
formula_corr <- ~ crp_log + diagnosis_class + age_group + bmi_class + antiTNF_vs_noBiologics + sex + Prednisolon + leucocytes + erythrocytes + thrombocytes

C = canCorPairs(formula_corr, coldata_antiTNF)
plotCorrMatrix(C) #saved as: "Output_files/Variance_partition/antiTNF/plots/correlation"
