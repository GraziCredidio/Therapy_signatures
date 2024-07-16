# EZE cohort: Therapy Signatures
  # Linear Mixed Model: Data preparation
  # Azathioprine vs No Systemic Therapy
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(tidyverse)
library(DESeq2)

folder <- "Cleaned_tables/models/aza"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading data ----
coldata <- read.csv("Cleaned_tables/EZECohort_coldata_clean_ord.txt", sep="\t")
counts <-  read.csv("Cleaned_tables/EZECohort_counts_ord.txt", sep="\t")

# Exploring, cleaning and formatting coldata ----
# Filtering for aza and no systemic IBD patients 
coldata <- coldata %>% 
  filter(diagnosis_class == "CD" | diagnosis_class == "UC") %>% 
  filter(Azathioprin == "1" | No_syst == "1") #214 IBD patients either use aza or no other systemic therapy

coldata <- coldata %>% 
  mutate(aza_vs_noSyst = case_when(
    (Azathioprin == "1")~ "aza",
    (No_syst == "1") ~ "no_syst")) %>% 
  relocate(aza_vs_noSyst, .after = No_syst)

# Transforming categorical columns as factor
coldata <- coldata %>% 
  mutate(across(c(diagnosis_class, age_group,
                  sex, bmi_class, aza_vs_noSyst), as.factor))


# Matching coldata and counts samples ----
# All patients
counts_all <- counts[,colnames(counts) %in% coldata$sample_id]
idx <- match(colnames(counts_all), rownames(coldata)) 
ord.coldata_all <- coldata[idx, ]

all(rownames(ord.coldata_all) %in% colnames(counts_all))
all(colnames(counts_all) == rownames(ord.coldata_all))

write.table(ord.coldata_all, "Cleaned_tables/models/aza/EZECohort_coldata_aza_allPatients.txt", sep = "\t", row.names = TRUE)

# Inactive disease patients: Filtering for inactive disease patients (R) with CRP < 5
coldata_inactive <- coldata %>% 
  filter(remission == "R") %>% 
  filter(crp < 5)

counts_inactive <- counts[,colnames(counts) %in% coldata_inactive$sample_id]
idx2 <- match(colnames(counts_inactive), rownames(coldata_inactive)) 
ord.coldata_inactive <- coldata_inactive[idx2, ]

all(rownames(ord.coldata_inactive) %in% colnames(counts_inactive))
all(colnames(counts_inactive) == rownames(ord.coldata_inactive))

write.table(ord.coldata_inactive, "Cleaned_tables/models/aza/EZECohort_coldata_aza.txt", sep = "\t", row.names = TRUE)

# Normalizing counts ----
deseq_norm <- function(counts, coldata){
  dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                       colData = coldata,
                                       design = ~ diagnosis_class + sex + bmi_class + age_group + aza_vs_noSyst)
  
  dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
  dds_counts <- estimateSizeFactors(dds_counts)
  return(dds_counts)
}

# all patients
dds_counts <- deseq_norm(counts = counts_all, coldata = ord.coldata_all)
norm_counts <- counts(dds_counts, normalized = TRUE)
df_normalized_counts <- as.data.frame(norm_counts)
write.table(df_normalized_counts, "Cleaned_tables/models/aza/normalized_counts_aza_allPatients.txt", sep = "\t",  quote = FALSE)

vst_counts <- vst(dds_counts, blind=TRUE) 
vst_counts <- assay(vst_counts)
write.table(vst_counts, "Cleaned_tables/models/aza/vst_counts_aza_allPatients.txt", sep = "\t",  quote = FALSE)

# inactive patients
dds_counts_inactive <- deseq_norm(counts = counts_inactive, coldata = ord.coldata_inactive)
norm_counts_inactive <- counts(dds_counts_inactive, normalized = TRUE)
df_normalized_counts_inactive <- as.data.frame(norm_counts_inactive)
write.table(df_normalized_counts_inactive, "Cleaned_tables/models/aza/normalized_counts_aza.txt", sep = "\t",  quote = FALSE)

vst_counts_inactive <- vst(dds_counts_inactive, blind=TRUE) 
vst_counts_inactive <- assay(vst_counts_inactive)
write.table(vst_counts_inactive, "Cleaned_tables/models/aza/vst_counts_aza.txt", sep = "\t",  quote = FALSE)
