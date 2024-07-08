# Modeling EZE cohort
# MaAslin2 (gene expression model): Prednisolon x No systemic therapy
# Data preparation

graphics.off()
rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("D:\\Documentos/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

library(tidyverse)
library(DESeq2)

# Loading data ----
coldata <- read.csv("Cleaned_tables/EZECohort_ord.coldata_maaslin_16.03.txt", sep="\t")
counts <-  read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_16.03.txt", sep="\t")

# Metadata ----
coldata <- coldata %>% 
  filter(Prednisolon == "1" | No_syst == "1") %>%  #420 patients either use prednisolon or no other systemic therapy
  filter(!(diagnosis_class == "Pso"))

coldata <- coldata %>% 
  mutate(pred_vs_noSyst = case_when(
    (Prednisolon == "1")~ "pred",
    (No_syst == "1") ~ "no_syst"))

# Formating coldata ----
## As factor
coldata <- coldata %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, pred_vs_noSyst), as.factor))

# Matching coldata and counts samples ----
counts <- counts[,colnames(counts) %in% coldata$sample_id]
idx <- match(colnames(counts), rownames(coldata)) 
ord.coldata <- coldata[idx, ]

all(rownames(ord.coldata) %in% colnames(counts))
all(colnames(counts) == rownames(ord.coldata))

write.table(ord.coldata, "Output_files/Maaslin2/Tables/coldata_maaslin2_pred_vs_noSyst_16.03.txt", sep = "\t", row.names = TRUE)

# DDS ----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = ord.coldata,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group2 +
                                       pred_vs_noSyst)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized = TRUE)
df_normalized_counts <- as.data.frame(normalized_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)


write.table(df_normalized_counts, "Output_files/DESeq2/Maaslin2/DESeq2_df_normalized_counts_pred_vs_noSyst_16.03.txt", sep = "\t",  quote = FALSE)

write.table(vst_counts_norm, "Output_files/DESeq2/Maaslin2/DESeq2_vst_counts_pred_vs_noSyst_29.03.txt", sep = "\t",  quote = FALSE)


######## Data prep: only remitters

# Loading data ----
coldata_pred <- read.csv("Output_files/Maaslin2/Tables/coldata_maaslin2_pred_vs_noSyst_16.03.txt", sep = "\t")
counts <-  read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_16.03.txt", sep="\t")

# Formating coldata ----
coldata <- coldata_pred %>% 
  filter(remission == "R")

## As factor
coldata <- coldata %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, pred_vs_noSyst), as.factor))

# Matching coldata and counts samples ----
counts <- counts[,colnames(counts) %in% coldata$sample_id]
idx <- match(colnames(counts), rownames(coldata)) 
ord.coldata <- coldata[idx, ]

all(rownames(ord.coldata) %in% colnames(counts))
all(colnames(counts) == rownames(ord.coldata))

write.table(ord.coldata, "Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_pred_vs_noSyst_R_17.03.txt", sep = "\t", row.names = TRUE)

# DDS ----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = ord.coldata,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group2 +
                                       pred_vs_noSyst)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized = TRUE)
df_normalized_counts <- as.data.frame(normalized_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

write.table(df_normalized_counts, "Output_files/DESeq2/Maaslin2/Remission/DESeq2_df_normalized_counts_pred_vs_noSyst_R_17.03.txt", sep = "\t",  quote = FALSE)

write.table(vst_counts_norm, "Output_files/DESeq2/Maaslin2/Remission/DESeq2_vst_counts_pred_vs_noSyst_R_29.03.txt", sep = "\t",  quote = FALSE)



######## Data prep: remitters that use ONLY aza as systemic therapy

# Loading data
coldata_R_pred <- read.table("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_pred_vs_noSyst_R_17.03.txt", sep = "\t")
counts <- read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_16.03.txt", sep="\t")
names(coldata_R_pred)
# Formating coldata 
## excluding patients that use systemic therapies other than Pred
coldata_R_onlyPred <- coldata_R_pred %>% 
  filter(!(MTX == 1 | X6.MP == 1 | Sulfasalazin ==1 | Chloroquin ==1 | Hydroxychloroquin ==1 | Azathioprin ==1 |
             Leflunomid ==1 | Ciclosporin ==1 | Mycophenolat.Mofetil ==1 | Cyclophosphamid ==1 | Tacrolimus ==1 |
             Apremilast ==1 | Mepacrine ==1))

table(coldata_R_onlyPred$pred_vs_noSyst)

## As factor
coldata_R_onlyPred <- coldata_R_onlyPred %>% 
  mutate(across(c(diagnosis_class, age_group2, 
                  sex, bmi_class, pred_vs_noSyst), as.factor))

# Matching coldata and counts samples ----
counts <- counts[,colnames(counts) %in% coldata_R_onlyPred$sample_id]
idx <- match(colnames(counts), rownames(coldata_R_onlyPred)) 
ord.coldata <- coldata_R_onlyPred[idx, ]

all(rownames(ord.coldata) %in% colnames(counts))
all(colnames(counts) == rownames(ord.coldata))

write.table(ord.coldata, "Output_files/Maaslin2/Tables/Remission/only_aza_pred/coldata_maaslin2_only_pred_R_26.04.txt", sep = "\t", row.names = TRUE)

# DDS ----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = ord.coldata,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group2 +
                                       pred_vs_noSyst)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized = TRUE)
df_normalized_counts <- as.data.frame(normalized_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

write.table(df_normalized_counts, "Output_files/DESeq2/Maaslin2/Remission/only_aza_pred/DESeq2_normalized_only_pred_R_26.04.txt", sep = "\t",  quote = FALSE)

write.table(vst_counts_norm, "Output_files/DESeq2/Maaslin2/Remission/only_aza_pred/DESeq2_vst_only_pred_R_26.04.txt", sep = "\t",  quote = FALSE)



######## Dataprep: inactive disease + crp < 5

coldata_R_pred <- read.table("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_pred_vs_noSyst_R_17.03.txt", sep = "\t")
counts <- read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_16.03.txt", sep="\t")
names(coldata_R_pred)
# Formating coldata 
coldata_R_pred_crp <- coldata_R_pred %>% 
  filter(crp < 5)

## As factor
coldata_R_pred_crp <- coldata_R_pred_crp %>% 
  mutate(across(c(diagnosis_class, age_group2, 
                  sex, bmi_class, pred_vs_noSyst), as.factor))

# Matching coldata and counts samples ----
counts <- counts[,colnames(counts) %in% coldata_R_pred_crp$sample_id]
idx <- match(colnames(counts), rownames(coldata_R_pred_crp)) 
ord.coldata <- coldata_R_pred_crp[idx, ]

all(rownames(ord.coldata) %in% colnames(counts))
all(colnames(counts) == rownames(ord.coldata))

write.table(ord.coldata, "Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_pred_R_crp_04.05.txt", sep = "\t", row.names = TRUE)

# DDS ----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = ord.coldata,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group2 +
                                       pred_vs_noSyst)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized = TRUE)
df_normalized_counts <- as.data.frame(normalized_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

write.table(df_normalized_counts, "Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_normalized_pred_R_crp_04.05.txt", sep = "\t",  quote = FALSE)

write.table(vst_counts_norm, "Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_pred_R_crp_04.05.txt", sep = "\t",  quote = FALSE)






