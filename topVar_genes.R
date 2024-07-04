# EZE cohort
# Top var genes - all and remitters

graphics.off()
rm(list = ls())

library(DESeq2)

# TVG all patients
counts <- read.table("Cleaned_tables/EZECohort_ord.counts_maaslin_16.03.txt", sep = "\t")
coldata <- read.table("Cleaned_tables/EZECohort_ord.coldata_maaslin_16.03.txt", sep = "\t")

# DESeq2 object
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = coldata,
                                     design = ~ 1)
dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- as.data.frame(assay(vst))

# Top variating genes
topVarGenes_idx <- head(order(rowVars(assay(vst)),decreasing=TRUE),2000) #5000 or 2000
topVarGenes <- vst_counts_norm[topVarGenes_idx,]
topVarGenes_names <- rownames(topVarGenes) 

write.table(topVarGenes, "Cleaned_tables/topVar_2000_genes.txt", sep = "\t", quote = FALSE)

# Loading files ----
topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVarGenes_names <- rownames(topVarGenes) 


# TVG remitters
coldata_R <- coldata %>% 
  filter(remission == "R" & crp < 5)
counts_R <- counts[,colnames(counts) %in% coldata_R$sample_id]

all(rownames(coldata_R) %in% colnames(counts_R))
all(colnames(counts_R) == rownames(coldata_R))

dds_counts <- DESeqDataSetFromMatrix(countData = counts_R,
                                     colData = coldata_R,
                                     design = ~ 1)
dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts_R), ]
dds_counts <- estimateSizeFactors(dds_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- as.data.frame(assay(vst))

# Top variating genes
topVarGenes_idx_R <- head(order(rowVars(assay(vst)),decreasing=TRUE),2000) #5000 or 2000
topVarGenes_R <- vst_counts_norm[topVarGenes_idx_R,]
topVarGenes_names_R <- rownames(topVarGenes_R) 

write.table(topVarGenes_R, "Cleaned_tables/topVar_R_2000_genes_28.06.txt", sep = "\t", quote = FALSE)

