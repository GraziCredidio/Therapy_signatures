# EZE cohort: Therapy Signatures
  # Linear Mixed Model: filtering significant genes
  # Author: Graziella Credidio

rm(list = ls())

library(tidyverse)

# Loading files ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = T, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

maaslin_analysis_R_crp <- c("antiTNF_vs_noBio", "aza_vs_noSyst", "pred_vs_noSyst")


for (i in maaslin_analysis_R_crp) {
  
  comparisons_result <- read.table(file.path(paste("Output_files/Maaslin2/maaslin2_results_", 
                                                   i, ".txt", sep = "")), sep = "\t")
  
  sig_results <- subset(comparisons_result, comparisons_result$qval <= 0.05, drop = FALSE) 
  sig_results <- sig_results[order(sig_results$qval),]
  sig_results$genes <- unique_ensg2gene[sig_results$feature, ]$hgnc_symbol
  
  write.table(sig_results, file = file.path(paste("Output_files/Maaslin2/significant_results/maaslin2_significant_results_",
                                                  i, ".txt", sep = "")), sep = "\t", quote = FALSE)
}
