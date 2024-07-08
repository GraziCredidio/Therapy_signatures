# EZE cohort: Therapy Signatures
  # Top 2000 Variying Genes (TVG) file creation
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(DESeq2)

# Loading data ----
counts <- read.table("Cleaned_tables/EZECohort_counts_ord.txt", sep = "\t")
coldata <- read.table("Cleaned_tables/EZECohort_coldata_clean_ord.txt", sep = "\t")

# Normalize counts: DESeq2 object ----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = coldata,
                                     design = ~ 1)
dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- as.data.frame(assay(vst))

# Top 2000 varying genes among all patients ----
topVarGenes_idx <- head(order(rowVars(assay(vst)),decreasing=TRUE),2000) 
topVarGenes <- vst_counts_norm[topVarGenes_idx,]
topVarGenes_names <- rownames(topVarGenes) 

# Saving data ----
write.table(topVarGenes, "Cleaned_tables/topVar_2000_genes.txt", sep = "\t", quote = FALSE)
