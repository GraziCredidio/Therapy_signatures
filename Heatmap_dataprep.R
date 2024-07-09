# EZE cohort: Therapy Signatures
  # Heatmaps: data preprocessing
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(DESeq2)
library(tidyverse)

# Loading files ----
coldata_antiTNF <- read.table("Cleaned_tables/models/antiTNF/EZECohort_coldata_antiTNF.txt", sep = "\t")
coldata_pred <- read.table("Cleaned_tables/models/pred/EZECohort_coldata_pred.txt", sep = "\t")
coldata_aza <- read.table("Cleaned_tables/models/aza/EZECohort_coldata_aza.txt", sep = "\t")
counts <- read.table("Cleaned_tables/EZECohort_counts_ord.txt")

# healthy subjects data
healthy_coldata <- read.table("Raw_tables/DZHK_EZE_col_data_wo_outliers.txt", sep = "\t", header = TRUE)
healthy_counts <- read.table("Raw_tables/DZHK_EZE_counts_wo_outliers.txt", sep = "\t")

# Manipulating healthy coldata ----
# Changing col names and recoding
healthy_coldata_filtered <- healthy_coldata %>% 
  filter(Diagnosis == "Healthy") %>% 
  dplyr::rename(sample_id = Sample_ID,
         age = Age,
         sex = Gender,
         diagnosis_class = Diagnosis) %>% 
  dplyr::select(!(Age_group))

healthy_coldata_filtered <- healthy_coldata_filtered %>% 
  mutate(sex = recode( sex,
                       "M" = "Male",
                       "W" = "Female"))

rownames(healthy_coldata_filtered) <- healthy_coldata_filtered$sample_id

# Adding age group column
healthy_coldata_filtered <- healthy_coldata_filtered %>% 
  mutate(age_group = cut(age, breaks=c(0,20,35,50,65,Inf), include.lowest = FALSE, labels = c("0_20", "21_35", "36_50", "51_65", "65_inf")))

# Merging coldatas ----
coldata_antiTNF_full <- full_join(coldata_antiTNF, healthy_coldata_filtered) 
rownames(coldata_antiTNF_full) <- coldata_antiTNF_full$sample_id

coldata_pred_full <- full_join(coldata_pred, healthy_coldata_filtered) 
rownames(coldata_pred_full) <- coldata_pred_full$sample_id

coldata_aza_full <- full_join(coldata_aza, healthy_coldata_filtered) 
rownames(coldata_aza_full) <- coldata_aza_full$sample_id

# Merging counts ----
# matching order of genes
idx_genes <- match(rownames(counts), rownames(healthy_counts))
healthy_counts <- healthy_counts[idx_genes,]

# filtering healthy counts to have only healthy individuals' counts
healthy_counts <- healthy_counts[,colnames(healthy_counts) %in% rownames(healthy_coldata_filtered)]

# merging tables
counts_full <- merge(counts, healthy_counts, by="row.names", all=TRUE)
rownames(counts_full) <- counts_full$Row.names
counts_full <- counts_full %>% 
  dplyr::select(!(Row.names))


# Data preprocessing for normalization ----
# mathcing patients order in coldata and counts 
counts_full_antiTNF <- counts_full[,colnames(counts_full) %in% coldata_antiTNF_full$sample_id]
idx_antiTNF <- match(colnames(counts_full_antiTNF), rownames(coldata_antiTNF_full)) 
ord.coldata_antiTNF_full <- coldata_antiTNF_full[idx_antiTNF, ]
ord.coldata_antiTNF_full <- ord.coldata_antiTNF_full %>%
  mutate_at('antiTNF_vs_noBiologics', ~replace_na(.,"Healthy"))
write.table(ord.coldata_antiTNF_full, "Cleaned_tables/heatmaps/hm_coldata_antiTNF.txt", sep = "\t", row.names = TRUE)

counts_full_pred <- counts_full[,colnames(counts_full) %in% coldata_pred_full$sample_id]
idx_pred <- match(colnames(counts_full_pred), rownames(coldata_pred_full)) 
ord.coldata_pred_full <- coldata_pred_full[idx_pred, ]
ord.coldata_pred_full <- ord.coldata_pred_full %>%
  mutate_at('pred_vs_noSyst', ~replace_na(.,"Healthy"))
write.table(ord.coldata_pred_full, "Cleaned_tables/heatmaps/hm_coldata_pred.txt", sep = "\t", row.names = TRUE)

counts_full_aza <- counts_full[,colnames(counts_full) %in% coldata_aza_full$sample_id]
idx_aza <- match(colnames(counts_full_aza), rownames(coldata_aza_full)) 
ord.coldata_aza_full <- coldata_aza_full[idx_aza, ]
ord.coldata_aza_full <- ord.coldata_aza_full %>%
  mutate_at('aza_vs_noSyst', ~replace_na(.,"Healthy"))
write.table(ord.coldata_aza_full, "Cleaned_tables/heatmaps/hm_coldata_aza.txt", sep = "\t", row.names = TRUE)

# Normalizing counts ---- 
normalized_full <- function(counts, ordColdata, pathNorm){
  dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                       colData = ordColdata,
                                       design = ~ 1)
  
  dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
  dds_counts <- estimateSizeFactors(dds_counts)
  normalized_counts <- counts(dds_counts, normalized = TRUE)
  df_normalized_counts <- as.data.frame(normalized_counts)

  write.table(df_normalized_counts, pathNorm, sep = "\t",  quote = FALSE)
}

normalized_full(counts = counts_full_antiTNF,
                ordColdata = ord.coldata_antiTNF_full,
                "Cleaned_tables/heatmaps/hm_normalized_counts_antiTNF.txt")

normalized_full(counts = counts_full_pred,
                ordColdata = ord.coldata_pred_full,
                "Cleaned_tables/heatmaps/hm_normalized_counts_pred.txt")

normalized_full(counts = counts_full_aza,
                ordColdata = ord.coldata_aza_full,
                "Cleaned_tables/heatmaps/hm_normalized_counts_aza.txt")
