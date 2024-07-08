# Modeling EZE cohort
# MaAslin2 (gene expression model): Anti-TNF x No Biologics
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
  filter(biologics_TNF == "anti_tnf" | biologics_TNF == "no_biologics") %>%  #549 patients
  dplyr::rename(antiTNF_vs_noBiologics = biologics_TNF) %>% 
  filter(!(diagnosis_class == "Arthrosis" | diagnosis_class == "SLE" | diagnosis_class == "Pso"))

# Formating coldata ----
## As factor
coldata <- coldata %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, antiTNF_vs_noBiologics), as.factor))

# Matching coldata and counts samples ----
counts <- counts[,colnames(counts) %in% coldata$sample_id]
idx <- match(colnames(counts), rownames(coldata)) 
ord.coldata <- coldata[idx, ]

all(rownames(ord.coldata) %in% colnames(counts))
all(colnames(counts) == rownames(ord.coldata))

write.table(ord.coldata, "Output_files/Maaslin2/Tables/coldata_maaslin2_antiTNF_vs_noBiologics_diag_20.03.txt", sep = "\t", row.names = TRUE)

# DDS ----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = ord.coldata,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group2 +
                                       antiTNF_vs_noBiologics)


dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized = TRUE)
df_normalized_counts <- as.data.frame(normalized_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)


write.table(df_normalized_counts, "Output_files/DESeq2/Maaslin2/DESeq2_df_normalized_counts_antiTNF_vs_noBiologics_diag_20.03.txt", sep = "\t",  quote = FALSE)

write.table(vst_counts_norm, "Output_files/DESeq2/Maaslin2/DESeq2_vst_counts_antiTNF_vs_noBiologics_29.03.txt", sep = "\t",  quote = FALSE)

###################### Data prep: without patients using systemic therapies
# Loading data ----
# coldata <- read.table("Output_files/Maaslin2/Tables/coldata_maaslin2_antiTNF_vs_noBiologics_16.03.txt", sep = "\t")
# counts <-  read.csv("Output_files/DESeq2/Maaslin2/DESeq2_df_normalized_counts_antiTNF_vs_noBiologics_16.03.txt", sep = "\t")

# Metadata with systemic therapy (mtx, aza, pred) column as dummy variable
# coldata2 <- coldata %>% 
#   mutate(syst_therapy = ifelse(Azathioprin == 1, 1,
#                                ifelse(MTX == 1, 1,
#                                       ifelse(Prednisolon == 1, 1, 0))))
# 
# table(coldata2$syst_therapy)

# Metadata using No_syst column
# coldata <- coldata %>% 
#   filter(No_syst == 1)
# 
# table(coldata$antiTNF_vs_noBiologics)
# 
# # Formating coldata ----
# ## As factor
# coldata <- coldata %>% 
#   mutate(across(c(diagnosis_class, age_group, age_group2, 
#                   sex, bmi_class, antiTNF_vs_noBiologics), as.factor))
# 
# # Matching coldata and counts samples ----
# counts <- counts[,colnames(counts) %in% coldata$sample_id]
# idx <- match(colnames(counts), rownames(coldata)) 
# ord.coldata <- coldata[idx, ]
# 
# all(rownames(ord.coldata) %in% colnames(counts))
# all(colnames(counts) == rownames(ord.coldata))
# 
# write.table(ord.coldata, "Output_files/Maaslin2/Tables/coldata_maaslin2_antiTNF_vs_noBiologics_NoSyst_10.03.txt", sep = "\t", row.names = TRUE)
# 
# # DDS ----
# dds_counts <- DESeqDataSetFromMatrix(countData = counts,
#                                      colData = ord.coldata,
#                                      design = ~ diagnosis_class + sex + bmi_class + age_group2 +
#                                        antiTNF_vs_noBiologics)
# 
# dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
# dds_counts <- estimateSizeFactors(dds_counts)
# normalized_counts <- counts(dds_counts, normalized = TRUE)
# df_normalized_counts <- as.data.frame(normalized_counts)
# write.table(df_normalized_counts, "Output_files/DESeq2/Maaslin2/DESeq2_df_normalized_counts_antiTNF_vs_noBiologics_NoSyst_10.03.txt", sep = "\t",  quote = FALSE)


################################################## Data prep: only remitters
# Loading data ----
coldata <- read.csv("Output_files/Maaslin2/Tables/coldata_maaslin2_antiTNF_vs_noBiologics_diag_20.03.txt", sep = "\t")
counts <-  read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_16.03.txt", sep="\t")

# coldata %>%
#   filter(remission == "R" & antiTNF_vs_noBiologics == "no_biologics") %>%
#   nrow()

# Formating coldata ----
coldata <- coldata %>% 
  filter(remission == "R")

## As factor
coldata <- coldata %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, antiTNF_vs_noBiologics), as.factor))

# Matching coldata and counts samples ----
counts <- counts[,colnames(counts) %in% coldata$sample_id]
idx <- match(colnames(counts), rownames(coldata)) 
ord.coldata <- coldata[idx, ]

all(rownames(ord.coldata) %in% colnames(counts))
all(colnames(counts) == rownames(ord.coldata))

write.table(ord.coldata, "Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_antiTNF_vs_noBio_diag_R_20.03.txt", sep = "\t", row.names = TRUE)

# DDS ----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = ord.coldata,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group2 +
                                       antiTNF_vs_noBiologics)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized = TRUE)
df_normalized_counts <- as.data.frame(normalized_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

write.table(df_normalized_counts, "Output_files/DESeq2/Maaslin2/Remission/DESeq2_df_normalized_counts_antiTNF_vs_noBio_diag_R_20.03.txt", sep = "\t",  quote = FALSE)

write.table(vst_counts_norm, "Output_files/DESeq2/Maaslin2/Remission/DESeq2_vst_counts_antiTNF_vs_noBiologics_R_29.03.txt", sep = "\t",  quote = FALSE)




###### Dataprep: inactive disease + crp < 5
coldata_noBio_R <- read.table("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_antiTNF_vs_noBio_diag_R_20.03.txt", sep = "\t")
counts <-  read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_16.03.txt", sep="\t")



# Formating coldata 
coldata_noBio_R_crp <- coldata_noBio_R %>% 
  filter(crp < 5)

coldata_noBio_R_crp <- coldata_noBio_R_crp %>% 
  mutate(across(c(diagnosis_class, age_group, age_group2, 
                  sex, bmi_class, antiTNF_vs_noBiologics), as.factor))

# Matching coldata and counts samples ----
counts <- counts[,colnames(counts) %in% coldata_noBio_R_crp$sample_id]
idx <- match(colnames(counts), rownames(coldata_noBio_R_crp)) 
ord.coldata <- coldata_noBio_R_crp[idx, ]

all(rownames(ord.coldata) %in% colnames(counts))
all(colnames(counts) == rownames(ord.coldata))

write.table(ord.coldata, "Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_antiTNF_R_crp_04.05.txt", sep = "\t", row.names = TRUE)

# DDS ----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = ord.coldata,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group2 +
                                       antiTNF_vs_noBiologics)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized = TRUE)
df_normalized_counts <- as.data.frame(normalized_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

write.table(df_normalized_counts, "Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_normalized_antiTNF_R_crp_04.05.txt", sep = "\t",  quote = FALSE)

write.table(vst_counts_norm, "Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_antiTNF_R_crp_04.05.txt", sep = "\t",  quote = FALSE)

