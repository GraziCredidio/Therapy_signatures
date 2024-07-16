# EZE cohort: Therapy Signatures
  # DNA methylation-linked genes: Correlation between DNAm intensity of sites located 5000 bp up- and downstream of the TSS DEGs
  # Azathioprine x No Systemic Therapy 
  # Author: Graziella Credidio

rm(list = ls())

folder <- "Output_files/Methylation/aza"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading packages ----
library(org.Hs.eg.db)
library(tibble)
library(tidyverse)

# Loading data ----
# Methylation data
methylation_data <- read.csv("Raw_tables/methylation/m_values_hg38_round.txt", header = TRUE, row.names = NULL, 
                             sep = '\t')
meth_sample_info <- read.csv("Raw_tables/methylation/DZHK_EZE_epicarray_metadata_with_remission_status_090223.csv",
                             header = TRUE, sep = ',')
gene_meth_sites <- read.csv("Raw_tables/methylation/gene_meth_sites_5000bp.txt", sep = '\t')

# Normalized gene counts and coldata
transcriptome_data <- read.table("Cleaned_tables/models/aza/normalized_counts_aza.txt", sep = "\t")
transcriptome_info_file <- read.table("Cleaned_tables/models/aza/EZECohort_coldata_aza.txt", sep = "\t")

# DEGs
degs <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_aza_vs_noSyst.txt", sep = '\t')
deg_gene_list <- degs$feature

# To add gene names
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Methylation data cleaning ----
# Removing NAs
which(is.na(methylation_data$row.names)) #10359 = NA
methylation_data <- methylation_data[-10359,]

# Filter methylation data for sites in meth_sample_info
methylation_data_filtered <- methylation_data[, colnames(methylation_data) %in% paste('X', as.character(meth_sample_info$Sample_ID), sep = '')]
rownames(methylation_data_filtered) <- methylation_data$row.names
names(methylation_data_filtered) <- gsub(pattern= "X", replacement = "", x=names(methylation_data_filtered))

# Matching order of sample info and methylation data filtered
meth_sample_info_filterd <- meth_sample_info[meth_sample_info$Sample_ID %in% colnames(methylation_data_filtered),]
idx <- match(colnames(methylation_data_filtered), meth_sample_info_filterd$Sample_ID) 
ord.meth_sample_info_filtered <- meth_sample_info_filterd[idx, ]

all(colnames(methylation_data_filtered) %in% ord.meth_sample_info_filtered$Sample_ID)
all(colnames(methylation_data_filtered) == ord.meth_sample_info_filtered$Sample_ID)

# Replacing sample_id by Study ID
colnames(methylation_data_filtered) <- ord.meth_sample_info_filtered$Study.ID

# Matching order of samples in transcriptome_info_file and transcriptome_data
#transcriptome_data <- transcriptome_data[, as.character(transcriptome_info_file$sample_id)] #normalized counts filtered by aza patients
idx2 <- match(colnames(transcriptome_data), transcriptome_info_file$sample_id) 
ord.transcriptome_info_file <- transcriptome_info_file[idx2, ]

all(colnames(transcriptome_data) %in% ord.transcriptome_info_file$sample_id)
all(colnames(transcriptome_data) == ord.transcriptome_info_file$sample_id)

colnames(transcriptome_data) <- ord.transcriptome_info_file$study_id

# Matching samples in transcriptome and methylome data
common_samples <- intersect(colnames(transcriptome_data), colnames(methylation_data_filtered))
transcriptome_data <- transcriptome_data[, common_samples]
methylation_data_filtered <- methylation_data_filtered[, common_samples]

# Read and filter gene and linked methylation sites info, remove sites on sex/mythoochondria chromosomes
gene_meth_sites <- gene_meth_sites[,-9] #removing column X
gene_meth_sites <- subset(gene_meth_sites, gene_meth_sites$Chr != 'chrX' & gene_meth_sites$Chr != 'chrY' & 
                            gene_meth_sites$Chr != 'chrM')
rownames(gene_meth_sites) <- gene_meth_sites$Gene_id

deg_gene_list <- intersect(deg_gene_list, rownames(gene_meth_sites))

deg_gene_meth_sites <- gene_meth_sites[as.character(deg_gene_list), ]
deg_gene_meth_sites <- subset(deg_gene_meth_sites, deg_gene_meth_sites$no_of_meth_sites > 0)

# Calculate correlation coefficient between the DEGs gene expression and the methylation intensity of TSS methylation sites----
deg_gene_meth_sites$meth_sites <- as.character(deg_gene_meth_sites$meth_sites)
deg_gene_meth_sites$no_of_meth_sites <- as.numeric(deg_gene_meth_sites$no_of_meth_sites)
gene_meth_site_correlation <- matrix(nrow = sum(deg_gene_meth_sites$no_of_meth_sites), ncol = 5)
rownum = 1

for (i in 1:nrow(deg_gene_meth_sites)) {
  meth_sites <- strsplit(deg_gene_meth_sites[i,8], split = ';')
  print(meth_sites)
  for (j in 1:length(meth_sites[[1]])) {
    site_id_info <- strsplit(meth_sites[[1]][j], split = ':')
    site_id <- site_id_info[[1]][1]
    distance <- site_id_info[[1]][3]
    print(site_id)
    
    if(site_id %in% rownames(methylation_data_filtered)){ #exclude sites that are not in methylation data
      corr_data_frame <- data.frame(meth=t(methylation_data_filtered[site_id,]), 
                                    expr=t(transcriptome_data[as.character(deg_gene_meth_sites[i,1]),]))
    }
    
    corr_data_frame <- corr_data_frame[complete.cases(corr_data_frame),] 
    rho <- cor(corr_data_frame[,1], corr_data_frame[,2], method="spearman")
    fdr <- 0
    for (k in 1:1000) {
      x <- sample(corr_data_frame[,1])
      y <- sample(corr_data_frame[,2])
      rand_rho <- cor(x, y, method="spearman")
      if(abs(rand_rho) >= abs(rho)){
        fdr <- fdr + 1
      }
    }
    fdr <- fdr/1000
    gene_meth_site_correlation[rownum, ] <- c(as.character(deg_gene_meth_sites[i,1]), site_id, distance, rho, fdr)
    rownum <- rownum + 1
  }
}

gene_meth_site_correlation <- as.data.frame(gene_meth_site_correlation)
colnames(gene_meth_site_correlation) <- c("Gene", "Site", "Distance_from_TSS", "Rho", "FDR")
gene_meth_site_correlation$FDR <- as.numeric(as.character(gene_meth_site_correlation$FDR)) 

write.table(gene_meth_site_correlation, file = "Output_files/Methylation/aza/aza_DEG_meth_site_correlation_5000bp_1000rep.txt", 
            quote = FALSE, sep = '\t')

# Extract significant DNAm linked DEGs  ----
sig_gene_meth_site_correlation <- gene_meth_site_correlation %>% 
  filter(FDR < 0.05)

# Save significantly correlated sites file
sig_gene_meth_site_correlation$Gene_name <- unique_ensg2gene[as.character(sig_gene_meth_site_correlation$Gene), "hgnc_symbol"]
write.table(sig_gene_meth_site_correlation, file = "Output_files/Methylation/aza/aza_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", 
            quote = FALSE, sep = '\t')
