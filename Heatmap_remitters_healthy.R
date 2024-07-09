# Therapy signature analysis - EZE cohort
# Heatmap: significant genes from Maaslin2
    # Plus: preprocessing and inclusion of healthy subjects on heatmaps
    # Date: 05.04


graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install("DESeq2")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install("ComplexHeatmap")
install.packages("circlize")

library(ComplexHeatmap)
library(DESeq2)
library(tidyverse)
library(circlize)

# Loading files ----
# raw EZE coldata
healthy_coldata <- read.table("Raw_tables/DZHK_EZE_col_data_wo_outliers.txt", sep = "\t", header = TRUE)

# coldatas used in maaslin2: remitters
coldata_noBio <- read.table("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_antiTNF_vs_noBio_diag_R_20.03.txt", sep = "\t")
coldata_pred_noSyst <- read.table("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_pred_vs_noSyst_R_17.03.txt", sep = "\t")
coldata_aza_noSyst <- read.table("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_aza_vs_noSyst_R_17.03.txt", sep = "\t")

# counts
healthy_counts <- read.table("Raw_tables/DZHK_EZE_counts_wo_outliers.txt", sep = "\t")
counts <-  read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_16.03.txt", sep="\t")

# Merging coldatas ----
# Changing col names and recoding healthy coldata
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
                       "W" = "Female"
  ))

rownames(healthy_coldata_filtered) <- healthy_coldata_filtered$sample_id

# adding columns that exist in the maaslin coldatas
healthy_coldata_filtered <- healthy_coldata_filtered %>% 
  mutate(age_group = cut(age, breaks=c(0,25,65,Inf), include.lowest = FALSE, labels = c("Young", "Adult", "Senior"))) %>% 
  relocate(age_group, .after = age)


healthy_coldata_filtered <- healthy_coldata_filtered %>% 
  mutate(age_group2 = cut(age, breaks=c(0,20,35,50,65,Inf), include.lowest = FALSE, labels = c("0_20", "21_35", "36_50", "51_65", "65_inf"))) %>% 
  relocate(age_group2, .after = age_group)


# Merging coldatas with healthy coldata
coldata_noBio_full <- full_join(coldata_noBio, healthy_coldata_filtered) 
rownames(coldata_noBio_full) <- coldata_noBio_full$sample_id

coldata_pred_noSyst_full <- full_join(coldata_pred_noSyst, healthy_coldata_filtered) 
rownames(coldata_pred_noSyst_full) <- coldata_pred_noSyst_full$sample_id

coldata_aza_noSyst_full <- full_join(coldata_aza_noSyst, healthy_coldata_filtered) 
rownames(coldata_aza_noSyst_full) <- coldata_aza_noSyst_full$sample_id


# Merging counts tables----
# matching order of genes
idx_genes <- match(rownames(counts), rownames(healthy_counts))
healthy_counts <- healthy_counts[idx_genes,]

# filtering healthy counts to have only healthy individuals' counts
healthy_counts <- healthy_counts[,colnames(healthy_counts) %in% rownames(healthy_coldata_filtered)]

# merging tables
counts_full <- merge(counts, healthy_counts, by="row.names", all=TRUE)  # merge by row names (by=0 or by="row.names")
rownames(counts_full) <- counts_full$Row.names
counts_full <- counts_full %>% 
  dplyr::select(!(Row.names))


# Creating normalized and vst counts for heatmaps ----
# mathcing patients order 
counts_full_noBio <- counts_full[,colnames(counts_full) %in% coldata_noBio_full$sample_id]
idx1 <- match(colnames(counts_full_noBio), rownames(coldata_noBio_full)) 
ord.coldata_noBio_full <- coldata_noBio_full[idx1, ]
ord.coldata_noBio_full <- ord.coldata_noBio_full %>%
  mutate_at('antiTNF_vs_noBiologics', ~replace_na(.,"Healthy"))
write.table(ord.coldata_noBio_full, "Output_files/Maaslin2/Tables/Remission/coldata_heatmap_antiTNF_noBio_healthy_R_04.04.txt", sep = "\t", row.names = TRUE)


counts_full_pred_noSyst <- counts_full[,colnames(counts_full) %in% coldata_pred_noSyst_full$sample_id]
idx2 <- match(colnames(counts_full_pred_noSyst), rownames(coldata_pred_noSyst_full)) 
ord.coldata_pred_noSyst_full <- coldata_pred_noSyst_full[idx2, ]
ord.coldata_pred_noSyst_full <- ord.coldata_pred_noSyst_full %>%
  mutate_at('pred_vs_noSyst', ~replace_na(.,"Healthy"))
write.table(ord.coldata_pred_noSyst_full, "Output_files/Maaslin2/Tables/Remission/coldata_heatmap_pred_noSyst_healthy_R_04.04.txt", sep = "\t", row.names = TRUE)


counts_full_aza_noSyst <- counts_full[,colnames(counts_full) %in% coldata_aza_noSyst_full$sample_id]
idx3 <- match(colnames(counts_full_aza_noSyst), rownames(coldata_aza_noSyst_full)) 
ord.coldata_aza_noSyst_full <- coldata_aza_noSyst_full[idx3, ]
ord.coldata_aza_noSyst_full <- ord.coldata_aza_noSyst_full %>%
  mutate_at('aza_vs_noSyst', ~replace_na(.,"Healthy"))
write.table(ord.coldata_aza_noSyst_full, "Output_files/Maaslin2/Tables/Remission/coldata_heatmap_aza_noSyst_healthy_R_04.04.txt", sep = "\t", row.names = TRUE)


# DESeq ----
normalized_vst_full <- function(counts, ordColdata, pathNorm, pathVst){
  dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                       colData = ordColdata,
                                       design = ~ 1)
  
  dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
  dds_counts <- estimateSizeFactors(dds_counts)
  normalized_counts <- counts(dds_counts, normalized = TRUE)
  df_normalized_counts <- as.data.frame(normalized_counts)
  
  vst <- vst(dds_counts, blind=TRUE) 
  vst_counts_norm <- assay(vst)
  
  write.table(df_normalized_counts, pathNorm, sep = "\t",  quote = FALSE)
  write.table(vst_counts_norm, pathVst, sep = "\t",  quote = FALSE)
}

normalized_vst_full(counts = counts_full_noBio, 
                    ordColdata = ord.coldata_noBio_full, 
                    "Output_files/DESeq2/Heatmap/Remission/DESeq2_normalized_counts_antiTNF_vs_noBiologics_healthy_R_04.04.txt",
                    "Output_files/DESeq2/Heatmap/Remission/DESeq2_vst_counts_antiTNF_vs_noBiologics_healthy_R_04.04.txt")


normalized_vst_full(counts = counts_full_pred_noSyst, 
                    ordColdata = ord.coldata_pred_noSyst_full, 
                    "Output_files/DESeq2/Heatmap/Remission/DESeq2_normalized_counts_pred_noSyst_healthy_R_04.04.txt",
                    "Output_files/DESeq2/Heatmap/Remission/DESeq2_vst_counts_pred_noSyst_healthy_R_04.04.txt")


normalized_vst_full(counts = counts_full_aza_noSyst, 
                    ordColdata = ord.coldata_aza_noSyst_full, 
                    "Output_files/DESeq2/Heatmap/Remission/DESeq2_normalized_counts_aza_noSyst_healthy_R_04.04.txt",
                    "Output_files/DESeq2/Heatmap/Remission/DESeq2_vst_counts_aza_noSyst_healthy_R_04.04.txt")




# Data preprocessing heatmap with healthy individuals ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = T, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# top 100 genes
top100_genes <- c("antiTNF_vs_noBio_R_vst_q0.05",
                  "aza_vs_noSyst_R_vst_q0.05",
                  "pred_vs_noSyst_R_vst_q0.05")

# vst counts
# R_vst_counts <- c("vst_counts_antiTNF_vs_noBiologics_healthy_R_04.04",  
#                    "vst_counts_aza_noSyst_healthy_R_04.04",
#                    "vst_counts_pred_noSyst_healthy_R_04.04") 


# normalized counts
R_norm_counts <- c("normalized_counts_antiTNF_vs_noBiologics_healthy_R_04.04",
                    "normalized_counts_aza_noSyst_healthy_R_04.04",
                    "normalized_counts_pred_noSyst_healthy_R_04.04")


for (i in 1:3){
  R_sig_genes <- top100_genes[i]
  R_deseq_genes <- R_norm_counts[i]
  
  top100 <- read.table(file.path(paste("Output_files/Maaslin2/Results/Remission/significant_results/top100/maaslin2_top100_", 
                                       R_sig_genes, ".txt", sep = "")), sep = "\t")
  
  # Creating top100 expression data matrix for the top100 genes
  counts <- read.table(file.path(paste("Output_files/DESeq2/Heatmap/Remission/DESeq2_",
                                                  R_deseq_genes, ".txt", sep = "" )), sep = "\t")
  
  top100_ed <- counts[top100$feature, ]
  top100_ed_scaled <- as.data.frame(t(base::scale(t(top100_ed))))
  
  top100_ed_scaled$gene <- unique_ensg2gene[rownames(top100_ed_scaled), ]$hgnc_symbol
  #top100_ed_scaled <- top100_ed_scaled[!(is.na(top100_ed_scaled$genes) | top100_ed_scaled$genes == ""), ] #excluding genes without name on Biomart
  rownames(top100_ed_scaled) <- make.names(top100_ed_scaled$gene, unique=TRUE)
  
  top100_ed_scaled <- top100_ed_scaled %>% #excluding gene names column
    dplyr::select(!("gene"))
  
  m.top100_ed_scaled <- data.matrix(top100_ed_scaled, rownames.force = NA)
  m.top100_ed_scaled <- na.omit(m.top100_ed_scaled) # removed rows with NAs
  
  write.table(m.top100_ed_scaled, file = file.path(paste("Output_files/Maaslin2/Results/Remission/significant_results/top100/Expression_matrix/heatmap_healthy/maaslin2_top100_normalized_expressionDataMatrix_",
                                                         R_deseq_genes, ".txt", sep = "")), sep = "\t", quote = FALSE)
  
  
}


# Heatmaps ----

## By comparison
coldata_maaslin <- c("antiTNF_noBio_healthy_R_04.04", 
                     "aza_noSyst_healthy_R_04.04",
                     "pred_noSyst_healthy_R_04.04")

norm_ed_scaled <- c("normalized_counts_antiTNF_vs_noBiologics_healthy_R_04.04",
                    "normalized_counts_aza_noSyst_healthy_R_04.04",
                    "normalized_counts_pred_noSyst_healthy_R_04.04")

split <- c("antiTNF_vs_noBiologics", "aza_vs_noSyst", "pred_vs_noSyst")


for (i in 1:3) { 
  coldata_loading <- coldata_maaslin[i]
  #deseq_genes <- vst_ed_scaled[i]
  deseq_genes <- norm_ed_scaled[i]
  
  coldata <- data.frame()
  coldata <- read.table(file.path(paste("Output_files/Maaslin2/Tables/Remission/coldata_heatmap_",
                                        coldata_loading, ".txt", sep = "")), sep = "\t")
  
  top100_ed_scaled <- data.frame()
  top100_ed_scaled <- read.table(file.path(paste("Output_files/Maaslin2/Results/Remission/significant_results/top100/Expression_matrix/heatmap_healthy/maaslin2_top100_normalized_expressionDataMatrix_",
                                                 deseq_genes, ".txt", sep = "")), sep = "\t")
  
  # Coldata preprocessing
  coldata <- coldata %>%
    mutate(Syst_therapy = recode(No_syst,
                                 "0" = "Yes",
                                 "1" = "No"))

  coldata <- coldata %>%
    mutate(Prednisolon = recode(Prednisolon,
                                 "0" = "No",
                                 "1" = "Yes"))

  coldata <- coldata %>%
    mutate(across(c(diagnosis_class, sex, remission, Syst_therapy, split[i]), as.factor))

  
  # Expression data matrix preprocessing
  m.top100_ed_scaled <- data.matrix(top100_ed_scaled, rownames.force = NA)
  
  # Heatmaps
  col_ha <- HeatmapAnnotation(
    df = data.frame(
      Sex = coldata$sex,
      BMI = coldata$bmi,
      Age = coldata$age,
      Biologics = coldata$biologics,
      Systemic_therapy = coldata$Syst_therapy,
      Diagnosis = coldata$diagnosis_class,
      Pred_dose = coldata$prednisolone_dose,
      Pred = coldata$Prednisolon
    ),
    col = list(Sex = c("Female" = "#b492c4",
                       "Male" = "#fbaa51"),
              Systemic_therapy = c("No" = "#ca493d",
                                  "Yes" = "#d7b66a"),
               Diagnosis = c("Arthrosis" =  "#66c2a5" ,
                             "CD" = "#fc8d62",
                             "UC" = "#8da0cb",
                             "PsA" = "#e78ac3",
                             "Pso" = "#a6d854", 
                             "RA" = "#ffd92f",
                             "SLE" = "#e5c494",
                             "Healthy" = "green3"),
              Biologics = c("biologics" = "#df536b",
                           "no_biologics" = "#f5c710"),
              Pred = c("No" = "#984464" ,
                     "Yes" = "#c0affb")
    ),
    CRP = anno_barplot(coldata$crp),
    annotation_name_gp= gpar(fontsize = 20, fontface = "bold"),
    annotation_legend_param = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 15),
                                   title_position = "topcenter", width = unit(20, "mm")),
    simple_anno_size = unit(1, "cm")
  )
  
  colors_hm <- colorRamp2(c(-4, 0, 4), c("#333b93", "white", "#b61728"))
  
  lgd = Legend(col_fun = colors_hm, legend_height = unit(6, "cm"), name = " ")

  
  # normalized
  png(file = file.path(paste("Output_files/Heatmap/vst_sig_genes_maaslin2_remission/Healthy/Heatmap_byComparison_norm_R_", split[i], ".png", sep = "")),
      width = 1728, height = 1494, units = "px")
  
  
  draw(ComplexHeatmap::Heatmap(m.top100_ed_scaled,
                               top_annotation = col_ha, show_column_names = FALSE, name = " ",
                               show_row_names = TRUE, border = TRUE, 
                               column_split = coldata[ , split[i]],
                               column_gap = unit(3, "mm"),
                               cluster_row_slices = TRUE, cluster_column_slices = TRUE,
                               show_column_dend = TRUE, show_row_dend = TRUE, 
                               col = colors_hm, column_title_gp = gpar(fontsize = 25, fontface = "bold"), 
                               heatmap_legend_param = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 15),
                                                           title_position = "topcenter")
  ))
  dev.off()
}


## By diagnosis and no split
coldata_maaslin <- c("antiTNF_noBio_healthy_R_04.04", 
                     "aza_noSyst_healthy_R_04.04",
                     "pred_noSyst_healthy_R_04.04")

norm_ed_scaled <- c("normalized_counts_antiTNF_vs_noBiologics_healthy_R_04.04",
                    "normalized_counts_aza_noSyst_healthy_R_04.04",
                    "normalized_counts_pred_noSyst_healthy_R_04.04")

for (i in 1:3) {
  coldata_loading <- coldata_maaslin[i]
  deseq_genes <- norm_ed_scaled[i]
  
  
  coldata_loading <- coldata_maaslin[i]
  #deseq_genes <- vst_ed_scaled[i]
  deseq_genes <- norm_ed_scaled[i]
  
  coldata <- data.frame()
  coldata <- read.table(file.path(paste("Output_files/Maaslin2/Tables/Remission/coldata_heatmap_",
                                        coldata_loading, ".txt", sep = "")), sep = "\t")
  
  top100_ed_scaled <- data.frame()
  top100_ed_scaled <- read.table(file.path(paste("Output_files/Maaslin2/Results/Remission/significant_results/top100/Expression_matrix/heatmap_healthy/maaslin2_top100_normalized_expressionDataMatrix_",
                                                 deseq_genes, ".txt", sep = "")), sep = "\t")
  
  # Coldata preprocessing
  coldata <- coldata %>%
    mutate(Syst_therapy = recode(No_syst,
                                 "0" = "Yes",
                                 "1" = "No"))
  
  coldata <- coldata %>%
    mutate(Prednisolon = recode(Prednisolon,
                                "0" = "No",
                                "1" = "Yes"))
  
  coldata <- coldata %>% 
    mutate(across(c(diagnosis_class, sex, remission, Syst_therapy), as.factor))
  
  # Expression data matrix preprocessing
  m.top100_ed_scaled <- data.matrix(top100_ed_scaled, rownames.force = NA)
  
  # Heatmaps
    col_ha2 <- HeatmapAnnotation(
      df = data.frame(
        Sex = coldata$sex,
        BMI = coldata$bmi,
        Age = coldata$age,
        Biologics = coldata$biologics,
        Systemic_therapy = coldata$Syst_therapy,
        Diagnosis = coldata$diagnosis_class,
        Pred_dose = coldata$prednisolone_dose,
        Pred = coldata$Prednisolon
      ),
      col = list(Sex = c("Female" = "#b492c4",
                         "Male" = "#fbaa51"),
                 Systemic_therapy = c("No" = "#ca493d",
                                      "Yes" = "#d7b66a"),
                 Diagnosis = c("Arthrosis" =  "#66c2a5" ,
                               "CD" = "#fc8d62",
                               "UC" = "#8da0cb",
                               "PsA" = "#e78ac3",
                               "Pso" = "#a6d854", 
                               "RA" = "#ffd92f",
                               "SLE" = "#e5c494",
                               "Healthy" = "green3"),
                 Biologics = c("biologics" = "#df536b",
                               "no_biologics" = "#f5c710"),
                 Pred = c("No" = "#984464" ,
                          "Yes" = "#c0affb")
      ),
      CRP = anno_barplot(coldata$crp),
      annotation_name_gp= gpar(fontsize = 20, fontface = "bold"),
      annotation_legend_param = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 15),
                                     title_position = "topcenter", width = unit(20, "mm")),
      simple_anno_size = unit(1, "cm")
    )
  
  colors_hm <- colorRamp2(c(-4, 0, 4), c("#333b93", "white", "#b61728"))
  
  lgd = Legend(col_fun = colors_hm, legend_height = unit(6, "cm"), name = "Z-score")
  
  
  ## no column split
  png(file = file.path(paste("Output_files/Heatmap/vst_sig_genes_maaslin2_remission/Healthy/Heatmap_noSplit_norm_R_", deseq_genes, ".png", sep = "")),
      width = 1728, height = 1494, units = "px")
  
  draw(ComplexHeatmap::Heatmap(m.top100_ed_scaled,
                               top_annotation = col_ha2, show_column_names = FALSE, name = " ",
                               show_row_names = TRUE, border = TRUE, column_gap = unit(3, "mm"),
                               cluster_row_slices = TRUE, cluster_column_slices = TRUE, 
                               show_column_dend = TRUE, show_row_dend = TRUE, 
                               col = colors_hm, column_title_gp = gpar(fontsize = 25, fontface = "bold"), 
                               heatmap_legend_param = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 15),
                                                           title_position = "topcenter")
  ))
  dev.off()
  
  
  ## split by diagnosis
  png(file = file.path(paste("Output_files/Heatmap/vst_sig_genes_maaslin2_remission/Healthy/Heatmap_byDiagnosis_norm_R_", deseq_genes, ".png", sep = "")),
      width = 1728, height = 1494, units = "px")

  draw(ComplexHeatmap::Heatmap(m.top100_ed_scaled,
                               top_annotation = col_ha2, show_column_names = FALSE, name = " ",
                               show_row_names = TRUE, border = TRUE, column_split = coldata$diagnosis_class, 
                               column_gap = unit(3, "mm"),
                               cluster_row_slices = TRUE, cluster_column_slices = TRUE, 
                               show_column_dend = TRUE, show_row_dend = TRUE, 
                               col = colors_hm, column_title_gp = gpar(fontsize = 25, fontface = "bold"), 
                               heatmap_legend_param = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 15),
                                                           title_position = "topcenter")
  ))
  
  dev.off()
  
}


############################### Remitters + crp ##################################
# Loading files ----
# raw EZE coldata
healthy_coldata <- read.table("Raw_tables/DZHK_EZE_col_data_wo_outliers.txt", sep = "\t", header = TRUE)

# coldatas used in maaslin2: remitters
coldata_noBio <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_antiTNF_R_crp_04.05.txt", sep = "\t")
coldata_pred_noSyst <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_pred_R_crp_04.05.txt", sep = "\t")
coldata_aza_noSyst <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_aza_R_crp_04.05.txt", sep = "\t")

# counts
healthy_counts <- read.table("Raw_tables/DZHK_EZE_counts_wo_outliers.txt", sep = "\t")
counts <-  read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_16.03.txt", sep="\t")

# Merging coldatas ----
# Changing col names and recoding healthy coldata
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
                       "W" = "Female"
  ))

rownames(healthy_coldata_filtered) <- healthy_coldata_filtered$sample_id

# adding columns that exist in the maaslin coldatas
healthy_coldata_filtered <- healthy_coldata_filtered %>% 
  mutate(age_group = cut(age, breaks=c(0,25,65,Inf), include.lowest = FALSE, labels = c("Young", "Adult", "Senior"))) %>% 
  relocate(age_group, .after = age)


healthy_coldata_filtered <- healthy_coldata_filtered %>% 
  mutate(age_group2 = cut(age, breaks=c(0,20,35,50,65,Inf), include.lowest = FALSE, labels = c("0_20", "21_35", "36_50", "51_65", "65_inf"))) %>% 
  relocate(age_group2, .after = age_group)


# Merging coldatas with healthy coldata
coldata_noBio_full <- full_join(coldata_noBio, healthy_coldata_filtered) 
rownames(coldata_noBio_full) <- coldata_noBio_full$sample_id

coldata_pred_noSyst_full <- full_join(coldata_pred_noSyst, healthy_coldata_filtered) 
rownames(coldata_pred_noSyst_full) <- coldata_pred_noSyst_full$sample_id

coldata_aza_noSyst_full <- full_join(coldata_aza_noSyst, healthy_coldata_filtered) 
rownames(coldata_aza_noSyst_full) <- coldata_aza_noSyst_full$sample_id


# Merging counts tables----
# matching order of genes
idx_genes <- match(rownames(counts), rownames(healthy_counts))
healthy_counts <- healthy_counts[idx_genes,]

# filtering healthy counts to have only healthy individuals' counts
healthy_counts <- healthy_counts[,colnames(healthy_counts) %in% rownames(healthy_coldata_filtered)]

# merging tables
counts_full <- merge(counts, healthy_counts, by="row.names", all=TRUE)  # merge by row names (by=0 or by="row.names")
rownames(counts_full) <- counts_full$Row.names
counts_full <- counts_full %>% 
  dplyr::select(!(Row.names))


# Creating normalized and vst counts for heatmaps ----
# mathcing patients order 
counts_full_noBio <- counts_full[,colnames(counts_full) %in% coldata_noBio_full$sample_id]
idx1 <- match(colnames(counts_full_noBio), rownames(coldata_noBio_full)) 
ord.coldata_noBio_full <- coldata_noBio_full[idx1, ]
ord.coldata_noBio_full <- ord.coldata_noBio_full %>%
  mutate_at('antiTNF_vs_noBiologics', ~replace_na(.,"Healthy"))
write.table(ord.coldata_noBio_full, "Output_files/Maaslin2/Tables/Remission_crp/coldata_heatmap_antiTNF_healthy_R_crp.txt", sep = "\t", row.names = TRUE)


counts_full_pred_noSyst <- counts_full[,colnames(counts_full) %in% coldata_pred_noSyst_full$sample_id]
idx2 <- match(colnames(counts_full_pred_noSyst), rownames(coldata_pred_noSyst_full)) 
ord.coldata_pred_noSyst_full <- coldata_pred_noSyst_full[idx2, ]
ord.coldata_pred_noSyst_full <- ord.coldata_pred_noSyst_full %>%
  mutate_at('pred_vs_noSyst', ~replace_na(.,"Healthy"))
write.table(ord.coldata_pred_noSyst_full, "Output_files/Maaslin2/Tables/Remission_crp/coldata_heatmap_pred_healthy_R_crp.txt", sep = "\t", row.names = TRUE)


counts_full_aza_noSyst <- counts_full[,colnames(counts_full) %in% coldata_aza_noSyst_full$sample_id]
idx3 <- match(colnames(counts_full_aza_noSyst), rownames(coldata_aza_noSyst_full)) 
ord.coldata_aza_noSyst_full <- coldata_aza_noSyst_full[idx3, ]
ord.coldata_aza_noSyst_full <- ord.coldata_aza_noSyst_full %>%
  mutate_at('aza_vs_noSyst', ~replace_na(.,"Healthy"))
write.table(ord.coldata_aza_noSyst_full, "Output_files/Maaslin2/Tables/Remission_crp/coldata_heatmap_aza_healthy_R_crp.txt", sep = "\t", row.names = TRUE)


# DESeq ----
normalized_vst_full <- function(counts, ordColdata, pathNorm, pathVst){
  dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                       colData = ordColdata,
                                       design = ~ 1)
  
  dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
  dds_counts <- estimateSizeFactors(dds_counts)
  normalized_counts <- counts(dds_counts, normalized = TRUE)
  df_normalized_counts <- as.data.frame(normalized_counts)
  
  vst <- vst(dds_counts, blind=TRUE) 
  vst_counts_norm <- assay(vst)
  
  write.table(df_normalized_counts, pathNorm, sep = "\t",  quote = FALSE)
  write.table(vst_counts_norm, pathVst, sep = "\t",  quote = FALSE)
}

normalized_vst_full(counts = counts_full_noBio, 
                    ordColdata = ord.coldata_noBio_full, 
                    "Output_files/DESeq2/Heatmap/Remission_crp/DESeq2_normalized_antiTNF_healthy_R_crp.txt",
                    "Output_files/DESeq2/Heatmap/Remission_crp/DESeq2_vst_antiTNF_healthy_R_crp.txt")


normalized_vst_full(counts = counts_full_pred_noSyst, 
                    ordColdata = ord.coldata_pred_noSyst_full, 
                    "Output_files/DESeq2/Heatmap/Remission_crp/DESeq2_normalized_pred_healthy_R_crp.txt",
                    "Output_files/DESeq2/Heatmap/Remission_crp/DESeq2_vst_pred_healthy_R_crp.txt")


normalized_vst_full(counts = counts_full_aza_noSyst, 
                    ordColdata = ord.coldata_aza_noSyst_full, 
                    "Output_files/DESeq2/Heatmap/Remission_crp/DESeq2_normalized_aza_healthy_R_crp.txt",
                    "Output_files/DESeq2/Heatmap/Remission_crp/DESeq2_vst_counts_aza_healthy_R_crp.txt")




# Data preprocessing heatmap with healthy individuals ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = T, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# top 100 genes
top100_genes <- c("antiTNF_vs_noBio_R_crp",
                  "aza_vs_noSyst_R_crp",
                  "pred_vs_noSyst_R_crp")

# normalized counts
R_norm_counts <- c("antiTNF_healthy_R_crp",
                   "aza_healthy_R_crp",
                   "pred_healthy_R_crp")


for (i in 1:3){
  R_sig_genes <- top100_genes[i]
  R_deseq_genes <- R_norm_counts[i]
  
  top100 <- read.table(file.path(paste("Output_files/Maaslin2/Results/Remission_crp/significant_results/top100/maaslin2_top100_", 
                                       R_sig_genes, ".txt", sep = "")), sep = "\t")
  
  # Creating top100 expression data matrix for the top100 genes
  counts <- read.table(file.path(paste("Output_files/DESeq2/Heatmap/Remission_crp/DESeq2_",
                                       R_deseq_genes, ".txt", sep = "" )), sep = "\t")
  
  top100_ed <- counts[top100$feature, ]
  top100_ed_scaled <- as.data.frame(t(base::scale(t(top100_ed))))
  
  top100_ed_scaled$gene <- unique_ensg2gene[rownames(top100_ed_scaled), ]$hgnc_symbol
  #top100_ed_scaled <- top100_ed_scaled[!(is.na(top100_ed_scaled$genes) | top100_ed_scaled$genes == ""), ] #excluding genes without name on Biomart
  rownames(top100_ed_scaled) <- make.names(top100_ed_scaled$gene, unique=TRUE)
  
  top100_ed_scaled <- top100_ed_scaled %>% #excluding gene names column
    dplyr::select(!("gene"))
  
  m.top100_ed_scaled <- data.matrix(top100_ed_scaled, rownames.force = NA)
  m.top100_ed_scaled <- na.omit(m.top100_ed_scaled) # removed rows with NAs
  
  write.table(m.top100_ed_scaled, file = file.path(paste("Output_files/Maaslin2/Results/Remission_crp/significant_results/top100/Expression_matrix/heatmap_healthy/maaslin2_top100_normalized_expressionDataMatrix_",
                                                         R_deseq_genes, ".txt", sep = "")), sep = "\t", quote = FALSE)
  
  
}


# Heatmaps ----
## By comparison
coldata_maaslin <- c("antiTNF_healthy_R_crp", 
                     "aza_healthy_R_crp",
                     "pred_healthy_R_crp")

norm_ed_scaled <- c("antiTNF_healthy_R_crp",
                    "aza_healthy_R_crp",
                    "pred_healthy_R_crp")

split <- c("antiTNF_vs_noBiologics", "aza_vs_noSyst", "pred_vs_noSyst")


for (i in 1:3) { 
  coldata_loading <- coldata_maaslin[i]
  #deseq_genes <- vst_ed_scaled[i]
  deseq_genes <- norm_ed_scaled[i]
  
  coldata <- data.frame()
  coldata <- read.table(file.path(paste("Output_files/Maaslin2/Tables/Remission_crp/coldata_heatmap_",
                                        coldata_loading, ".txt", sep = "")), sep = "\t")
  
  top100_ed_scaled <- data.frame()
  top100_ed_scaled <- read.table(file.path(paste("Output_files/Maaslin2/Results/Remission_crp/significant_results/top100/Expression_matrix/heatmap_healthy/maaslin2_top100_normalized_expressionDataMatrix_",
                                                 deseq_genes, ".txt", sep = "")), sep = "\t")
  
  # Coldata preprocessing
  coldata <- coldata %>%
    mutate(Syst_therapy = recode(No_syst,
                                 "0" = "Yes",
                                 "1" = "No"))
  
  coldata <- coldata %>%
    mutate(Prednisolon = recode(Prednisolon,
                                "0" = "No",
                                "1" = "Yes"))
  
  coldata <- coldata %>%
    mutate(across(c(diagnosis_class, sex, remission, Syst_therapy, split[i]), as.factor))
  
  
  # Expression data matrix preprocessing
  m.top100_ed_scaled <- data.matrix(top100_ed_scaled, rownames.force = NA)
  
  # Heatmaps
  col_ha <- HeatmapAnnotation(
    df = data.frame(
      Sex = coldata$sex,
      BMI = coldata$bmi,
      Age = coldata$age,
      Biologics = coldata$biologics,
      Systemic_therapy = coldata$Syst_therapy,
      Diagnosis = coldata$diagnosis_class,
      Pred_dose = coldata$prednisolone_dose,
      Pred = coldata$Prednisolon
    ),
    col = list(Sex = c("Female" = "#b492c4",
                       "Male" = "#fbaa51"),
               Systemic_therapy = c("No" = "#ca493d",
                                    "Yes" = "#d7b66a"),
               Diagnosis = c("Arthrosis" =  "#66c2a5" ,
                             "CD" = "#fc8d62",
                             "UC" = "#8da0cb",
                             "PsA" = "#e78ac3",
                             "Pso" = "#a6d854", 
                             "RA" = "#ffd92f",
                             "SLE" = "#e5c494",
                             "Healthy" = "green3"),
               Biologics = c("biologics" = "#df536b",
                             "no_biologics" = "#f5c710"),
               Pred = c("No" = "#984464" ,
                        "Yes" = "#c0affb")
    ),
    CRP = anno_barplot(coldata$crp),
    annotation_name_gp= gpar(fontsize = 20, fontface = "bold"),
    annotation_legend_param = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 15),
                                   title_position = "topcenter", width = unit(20, "mm")),
    simple_anno_size = unit(1, "cm")
  )
  
  colors_hm <- colorRamp2(c(-4, 0, 4), c("#333b93", "white", "#b61728"))
  
  lgd = Legend(col_fun = colors_hm, legend_height = unit(6, "cm"), name = " ")
  
  
  # normalized
  png(file = file.path(paste("Output_files/Heatmap/vst_sig_genes_maaslin2_remission_crp/Heatmap_byComparison_norm_R_crp", split[i], ".png", sep = "")),
      width = 1728, height = 1494, units = "px")
  
  
  draw(ComplexHeatmap::Heatmap(m.top100_ed_scaled,
                               top_annotation = col_ha, show_column_names = FALSE, name = " ",
                               show_row_names = TRUE, border = TRUE, 
                               column_split = coldata[ , split[i]],
                               column_gap = unit(3, "mm"),
                               cluster_row_slices = TRUE, cluster_column_slices = TRUE,
                               show_column_dend = TRUE, show_row_dend = TRUE, 
                               col = colors_hm, column_title_gp = gpar(fontsize = 25, fontface = "bold"), 
                               heatmap_legend_param = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 15),
                                                           title_position = "topcenter")
  ))
  dev.off()
}
