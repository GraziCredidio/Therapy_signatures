# EZE cohort: Therapy Signatures
  # Linear Mixed Model: Data preparation
  # AntiTNF vs No Biologics
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(tidyverse)
library(DESeq2)

folder <- "Cleaned_tables/models/antiTNF"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading data ----
coldata <- read.csv("Cleaned_tables/EZECohort_coldata_clean_ord.txt", sep="\t")
counts <-  read.csv("Cleaned_tables/EZECohort_counts_ord.txt", sep="\t")

# Cleaning and formatting coldata ----
# Filtering for anti-TNF and no Biologics patients and removing Arthrosis, SLE and Pso diagnosis
coldata <- coldata %>% 
  filter(biologics_TNF == "anti_tnf" | biologics_TNF == "no_biologics") %>% 
  dplyr::rename(antiTNF_vs_noBiologics = biologics_TNF) %>% 
  filter(!(diagnosis_class == "Arthrosis" | diagnosis_class == "SLE" | diagnosis_class == "Pso"))

# Filtering for patients with CRP < 5
coldata <- coldata %>% 
  filter(crp < 5)

# Transforming categorical columns as factor
coldata <- coldata %>% 
  mutate(across(c(diagnosis_class, age_group, 
                  sex, bmi_class, antiTNF_vs_noBiologics), as.factor))

# Matching coldata and counts samples ----
counts <- counts[,colnames(counts) %in% coldata$sample_id]
idx <- match(colnames(counts), rownames(coldata)) 
ord.coldata <- coldata[idx, ]

all(rownames(ord.coldata) %in% colnames(counts))
all(colnames(counts) == rownames(ord.coldata))

# Saving anti-TNF vs no biologics coldata
write.table(ord.coldata, "Cleaned_tables/models/antiTNF/EZECohort_coldata_antiTNF.txt", sep = "\t", row.names = TRUE)

# Normalizing counts ----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = ord.coldata,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group + antiTNF_vs_noBiologics)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized = TRUE)
df_normalized_counts <- as.data.frame(normalized_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

write.table(df_normalized_counts, "Cleaned_tables/models/antiTNF/normalized_counts_antiTNF.txt", sep = "\t",  quote = FALSE)
write.table(vst_counts_norm, "Cleaned_tables/models/antiTNF/vst_counts_antiTNF.txt", sep = "\t",  quote = FALSE)