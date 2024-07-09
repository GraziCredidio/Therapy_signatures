# Modeling EZE cohort
# significant results maaslin: reading and filtering
# date: 20.03

graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

library(tidyverse)

# Loading files ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = T, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

################################################### All patients ###############################################################
maaslin_analysis <- c("antiTNF_vs_noBio_vst",  "pred_vs_noPred_vst", "pred_vs_noSyst_vst", "antiTNF_vs_nonAntiTNF_vst",
                      "aza_vs_noAza_vst", "aza_vs_noSyst_vst", "mtx_vs_noMtx_vst", "mtx_vs_noSyst_vst") 

for (i in maaslin_analysis) {
  
  comparisons_result <- read.table(file.path(paste("Output_files/Maaslin2/Results/maaslin2_results_", 
                                                   i, ".txt", sep = "")), sep = "\t")
 
  sig_results <- subset(comparisons_result, comparisons_result$qval <= 0.05, drop = FALSE) 
  sig_results <- sig_results[order(sig_results$qval),]
  sig_results$genes <- unique_ensg2gene[sig_results$feature, ]$hgnc_symbol

  write.table(sig_results, file = file.path(paste("Output_files/Maaslin2/significant_results/maaslin2_significant_results_",
                                                                i, "_q0.05.txt", sep = "")), sep = "\t", quote = FALSE)
}


## qvalue threshold 0.1 (TNFs and mtx)
result_mtx_noSyst <- read.table("Output_files/Maaslin2/Results/maaslin2_results_mtx_vs_noSyst_vst.txt", sep = "\t")
sig_genes_mtx_noSyst <- subset(result_mtx_noSyst, result_mtx_noSyst$qval <= 0.1, drop = FALSE) #2297 genes
sig_genes_mtx_noSyst <- sig_genes_mtx_noSyst[order(sig_genes_mtx_noSyst$qval),]
sig_genes_mtx_noSyst$genes <- unique_ensg2gene[sig_genes_mtx_noSyst$feature, ]$hgnc_symbol
write.table(sig_genes_mtx_noSyst, file = "Output_files/Maaslin2/significant_results/maaslin2_significant_results_mtx_vs_noSyst_vst_q0.1.txt", sep = "\t")



result_mtx_noMtx <- read.table("Output_files/Maaslin2/Results/maaslin2_results_mtx_vs_noMtx_vst.txt", sep = "\t")
sig_genes_mtx_noMtx <- subset(result_mtx_noMtx, result_mtx_noMtx$qval <= 0.1, drop = FALSE) #57 genes
sig_genes_mtx_noMtx <- sig_genes_mtx_noMtx[order(sig_genes_mtx_noMtx$qval),]
sig_genes_mtx_noMtx$genes <- unique_ensg2gene[sig_genes_mtx_noMtx$feature, ]$hgnc_symbol
write.table(sig_genes_mtx_noMtx, file = "Output_files/Maaslin2/significant_results/maaslin2_significant_results_mtx_vs_noMtx_vst_q0.1.txt", sep = "\t")


result_antiTNF_noBio <- read.table("Output_files/Maaslin2/Results/maaslin2_results_antiTNF_vs_noBio_vst.txt", sep = "\t")
sig_genes_antiTNF_noBio <- subset(result_antiTNF_noBio, result_antiTNF_noBio$qval <= 0.1, drop = FALSE) #7536 genes
sig_genes_antiTNF_noBio <- sig_genes_antiTNF_noBio[order(sig_genes_antiTNF_noBio$qval),]
sig_genes_antiTNF_noBio$genes <- unique_ensg2gene[sig_genes_antiTNF_noBio$feature, ]$hgnc_symbol
write.table(sig_genes_antiTNF_noBio, file = "Output_files/Maaslin2/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_vst_q0.1.txt", sep = "\t")

result_antiTNF_nonAntiTNF <- read.table("Output_files/Maaslin2/Results/maaslin2_results_antiTNF_vs_nonAntiTNF_vst.txt", sep = "\t")
sig_genes_antiTNF_nonAntiTNF <- subset(result_antiTNF_nonAntiTNF, result_antiTNF_nonAntiTNF$qval <= 0.1, drop = FALSE) #4 genes
sig_genes_antiTNF_nonAntiTNF <- sig_genes_antiTNF_nonAntiTNF[order(sig_genes_antiTNF_nonAntiTNF$qval),]
sig_genes_antiTNF_nonAntiTNF$genes <- unique_ensg2gene[sig_genes_antiTNF_nonAntiTNF$feature, ]$hgnc_symbol
write.table(sig_genes_antiTNF_nonAntiTNF, file = "Output_files/Maaslin2/significant_results/maaslin2_significant_results_antiTNF_vs_nonAntiTNF_vst_q0.1.txt", sep = "\t")




################################################### Remitters ###############################################################

maaslin_analysis_remission <- c("antiTNF_vs_noBio_R_vst", "aza_vs_noAza_R_vst", "aza_vs_noSyst_R_vst",
                                "mtx_vs_noMtx_R_vst", "mtx_vs_noSyst_R_vst", "pred_vs_noPred_R_vst",
                                 "pred_vs_noSyst_R_vst")



for (i in maaslin_analysis_remission) {
  
  comparisons_result <- read.table(file.path(paste("Output_files/Maaslin2/Results/Remission/maaslin2_results_", 
                                                   i, ".txt", sep = "")), sep = "\t")
  
  sig_results <- subset(comparisons_result, comparisons_result$qval <= 0.05, drop = FALSE) 
  sig_results <- sig_results[order(sig_results$qval),]
  sig_results$genes <- unique_ensg2gene[sig_results$feature, ]$hgnc_symbol
  
  write.table(sig_results, file = file.path(paste("Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_",
                                                  i, "_q0.05.txt", sep = "")), sep = "\t", quote = FALSE)
}


## qvalue threshold 0.1
R_result_mtx_noSyst <- read.table("Output_files/Maaslin2/Results/Remission/maaslin2_results_mtx_vs_noSyst_R_vst.txt", sep = "\t")
R_sig_genes_mtx_noSyst <- subset(R_result_mtx_noSyst, R_result_mtx_noSyst$qval <= 0.1, drop = FALSE) #0

R_result_mtx_noMtx <- read.table("Output_files/Maaslin2/Results/Remission/maaslin2_results_mtx_vs_noMtx_R_vst.txt", sep = "\t")
R_sig_genes_mtx_noMtx <- subset(R_result_mtx_noMtx, R_result_mtx_noMtx$qval <= 0.1, drop = FALSE) #0


# Loading significant genes files ----
# All patients
setwd("C:/Documents/Masters thesis/EZE_cohort/Output_files/Maaslin2/significant_results") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort/Output_files/Maaslin2/significant_results") #pc

files = list.files(pattern="*.txt")
myfiles = lapply(files, read.delim)

noBio <- myfiles[[1]]
nonTNF <- myfiles[[2]]

aza_noAza <- myfiles[[3]]
aza_noSyst <- myfiles[[4]]

mtx_noMtx <- myfiles[[5]]
mtx_noSyst <- myfiles[[6]]

pred_noPred <- myfiles[[7]]
pred_noSyst <- myfiles[[8]]

# Remitters
setwd("C:/Documents/Masters thesis/EZE_cohort/Output_files/Maaslin2/Results/Remission/significant_results")

files_R = list.files(pattern="*.txt")
myfiles_R = lapply(files_R, read.delim)

R_noBio <- myfiles_R[[1]]

R_aza_noAza <- myfiles_R[[2]]
R_aza_noSyst <- myfiles_R[[3]]

R_mtx_noMtx <- myfiles_R[[4]]
R_mtx_noSyst <- myfiles_R[[5]]

R_pred_noPred <- myfiles_R[[6]]
R_pred_noSyst <- myfiles_R[[7]]


################################################### Remitters + CRP ###############################################################

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = T, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

maaslin_analysis_R_crp <- c("antiTNF_vs_noBio_R_crp", "aza_vs_noSyst_R_crp", "pred_vs_noSyst_R_crp")


for (i in maaslin_analysis_R_crp) {
  
  comparisons_result <- read.table(file.path(paste("Output_files/Maaslin2/Results/Remission_crp/maaslin2_results_", 
                                                   i, ".txt", sep = "")), sep = "\t")
  
  sig_results <- subset(comparisons_result, comparisons_result$qval <= 0.05, drop = FALSE) 
  sig_results <- sig_results[order(sig_results$qval),]
  sig_results$genes <- unique_ensg2gene[sig_results$feature, ]$hgnc_symbol
  
  write.table(sig_results, file = file.path(paste("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_",
                                                  i, ".txt", sep = "")), sep = "\t", quote = FALSE)
}

# Loading significant genes files
setwd("C:/Documents/Masters thesis/EZE_cohort/Output_files/Maaslin2/Results/Remission/significant_results")
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort/Output_files/Maaslin2/Results/Remission_crp/significant_results")

files_R_crp = list.files(pattern="*.txt")
myfiles_R_crp = lapply(files_R_crp, read.delim)

noBio_R_crp <- myfiles_R_crp[[1]]

aza_R_crp <- myfiles_R_crp[[2]]

pred_R_crp <- myfiles_R_crp[[3]]

