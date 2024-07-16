# EZE cohort: Therapy Signatures
  # Linear Mixed Models: filtering significant genes
  # Inactive disease patients and all patients models
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(tidyverse)

# Loading files ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = T, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

maaslin_analysis_inactive <- c("antiTNF_vs_noBio", "aza_vs_noSyst", "pred_vs_noSyst")
maaslin_analysis_all <- c("antiTNF_vs_noBio_allPatients", "aza_vs_noSyst_allPatients", "pred_vs_noSyst_allPatients")

# Extraction of significant genes from maaslin2 results
significant_genes <- function(analysis, source_table_path, output_table_path){
  for (i in analysis) {
    comparisons_result <- read.table(file.path(paste(source_table_path, 
                                                     i, ".txt", sep = "")), sep = "\t")
    
    sig_results <- subset(comparisons_result, comparisons_result$qval <= 0.05, drop = FALSE) 
    sig_results <- sig_results[order(sig_results$qval),]
    sig_results$genes <- unique_ensg2gene[sig_results$feature, ]$hgnc_symbol
    
    write.table(sig_results, file = file.path(paste(output_table_path,
                                                    i, ".txt", sep = "")), sep = "\t", quote = FALSE)
  }
}

# inactive disease patients
significant_genes(maaslin_analysis_inactive, 
                  "Output_files/Maaslin2/results/inactive/maaslin2_results_",
                  "Output_files/Maaslin2/significant_results/maaslin2_significant_results_")

# all patients
significant_genes(maaslin_analysis_all, 
                  "Output_files/Maaslin2/results/allPatients/maaslin2_results_",
                  "Output_files/Maaslin2/significant_results/allPatients/maaslin2_significant_results_")
