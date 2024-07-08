# EZE cohort: Therapy Signatures
  # Samples ID matching: col data and count data
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(tidyverse)

# Loading data ----
coldata_EZE <- read.csv("Cleaned_tables/EZECohort_coldata_clean.txt", sep="\t")
counts <- read.csv("Raw_tables/all_merged_gene_counts_EZE.txt",  sep = "\t")

seq_id <- read.csv("Raw_tables/DZHK_and_EZE_Sequenced_sample_and_pat_info_090223.txt", sep = "\t")
popgen_id <- read.csv("Raw_tables/Study_ID_Popgen_ID_match.txt", sep = "\t")

# Extracting patient ID included in coldata that were sequenced ----
popgen_id_included <- popgen_id[popgen_id$Study.ID %in% coldata_EZE$study_id,]
seq_id_included <- seq_id[seq_id$Study.ID %in% coldata_EZE$study_id,]

setdiff(coldata_EZE$study_id, seq_id_included$Study.ID) #no sequencing: "EZE097" "EZE240" "EZE273" "EZE498" "EZE508" "EZE514" "EZE636"

coldata_EZE <- coldata_EZE[coldata_EZE$study_id %in% seq_id_included$Study.ID, ] #excluding patients without sequencing data

patient_matching <- seq_id_included %>% 
  select("Sample_ID", "Patient", "Study.ID") %>% 
  dplyr::rename(study_id = "Study.ID",
                sample_id = "Sample_ID")

# Matching coldata study_id with sample ID ----
coldata_EZE <- inner_join(coldata_EZE, patient_matching, by= "study_id") 

coldata_EZE <- coldata_EZE %>% 
  relocate(c(sample_id, Patient), .after = study_id) 

# Matching coldata and counts data IDs ----
for (col in 1:ncol(counts)){
  colnames(counts)[col] <-  sub(".L.*", "", colnames(counts)[col]) }

coldata_id <- coldata_EZE$sample_id
counts_id <- colnames(counts)

diff <- setdiff(coldata_id, counts_id) #exclude "H18681" = EZE536

coldata_EZE <- coldata_EZE %>% 
  filter(!(sample_id == "H18681"))

rownames(coldata_EZE) <- coldata_EZE$sample_id
counts <- counts[,colnames(counts) %in% coldata_EZE$sample_id] #maching counts sample_ids and coldata sample_id
idx <- match(colnames(counts), rownames(coldata_EZE)) 
ord.coldata <- coldata_EZE[idx, ]

all(rownames(ord.coldata) %in% colnames(counts))
all(colnames(counts) == rownames(ord.coldata))

write.table(ord.coldata, "Cleaned_tables/EZECohort_coldata_clean_ord.txt", sep="\t", row.names = TRUE)
write.table(counts, "Cleaned_tables/EZECohort_counts_ord.txt", sep="\t", row.names = TRUE)