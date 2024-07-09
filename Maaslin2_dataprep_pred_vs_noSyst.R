# EZE cohort: Therapy Signatures
  # Linear Mixed Model: Data preparation
  # Prednisolone vs No Systemic Therapy
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(tidyverse)
library(DESeq2)

folder <- "Cleaned_tables/models/pred"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading data ----
coldata <- read.csv("Cleaned_tables/EZECohort_coldata_clean_ord.txt", sep="\t")
counts <-  read.csv("Cleaned_tables/EZECohort_counts_ord.txt", sep="\t")

# Exploring, cleaning and formatting coldata ----
# Filtering for pred and no systemic patients and excluding Pso patients 
coldata <- coldata %>% 
  filter(Prednisolon == "1" | No_syst == "1") %>%  #420 patients either use prednisolon or no other systemic therapy
  filter(!(diagnosis_class == "Pso"))

coldata <- coldata %>% 
  mutate(pred_vs_noSyst = case_when(
    (Prednisolon == "1")~ "pred",
    (No_syst == "1") ~ "no_syst"))

# Filtering for remitters with CRP < 5
coldata <- coldata %>% 
  filter(remission == "R") %>% 
  filter(crp < 5)

# Transforming categorical columns as factor
coldata <- coldata %>% 
  mutate(across(c(diagnosis_class, age_group, 
                  sex, bmi_class, pred_vs_noSyst), as.factor))

# Matching coldata and counts samples ----
counts <- counts[,colnames(counts) %in% coldata$sample_id]
idx <- match(colnames(counts), rownames(coldata)) 
ord.coldata <- coldata[idx, ]

all(rownames(ord.coldata) %in% colnames(counts))
all(colnames(counts) == rownames(ord.coldata))

write.table(ord.coldata, "Cleaned_tables/models/pred/EZECohort_coldata_pred.txt", sep = "\t", row.names = TRUE)

# Normalizing counts ----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = ord.coldata,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group +
                                       pred_vs_noSyst)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized = TRUE)
df_normalized_counts <- as.data.frame(normalized_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

write.table(df_normalized_counts, "Cleaned_tables/models/pred/normalized_counts_pred.txt", sep = "\t",  quote = FALSE)
write.table(vst_counts_norm, "Cleaned_tables/models/pred/vst_counts_pred.txt", sep = "\t",  quote = FALSE)
