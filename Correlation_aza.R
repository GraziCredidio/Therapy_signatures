# EZE cohort 
# Correlation clin parameters and aza gene expression (inactive)
# Date: 08.05


graphics.off()
rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

install.packages("ggpubr")
library("ggpubr")
library(tidyverse)


# Remitters
coldata_R_aza <- read.table("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_aza_vs_noSyst_R_17.03.txt", sep = "\t")
vst_counts_R_aza <- read.table("Output_files/DESeq2/Maaslin2/Remission/DESeq2_vst_counts_aza_vs_noSyst_R_29.03.txt", sep = "\t")
sig_genes_R_aza <- read.csv("Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_vst_q0.05.txt",
                            sep = "\t")
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting aza_dose, leukocytes, thrombocytes and erythrocytes from (updated) redcap coldata 
redcap_patients <- redcap[redcap$study_id %in% coldata_R_aza$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, aza_dose, leucocytes, erythrocytes, thrombocytes, endoscopy_performed, 
                partial_mayo, complete_mayo, hbi, cdai, sescd, endo_mayo)

idx <- match(coldata_R_aza$study_id, redcap_patients$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_aza <- left_join(coldata_R_aza, ord_redcap_patients, by = "study_id")
coldata_R_aza <- coldata_R_aza %>% 
  filter(Azathioprin == 1)

# Checking number of patients with info of clinical paramenters
CD <- coldata_R_aza %>% #16 CD
  filter(diagnosis_class == "CD") # 11 patients without sescd (no endoscopy); 1 without aza dose info (exclude)
table(is.na(CD$aza_dose))

# coldata_R_aza <- coldata_R_aza  %>% # excluding 1 CD without aza dose info (exclude)
#   filter(!(is.na(aza_dose)))

UC <- coldata_R_aza %>% #8 UC
  filter(diagnosis_class == "UC") # all have endomayo
table(UC$endoscopy_performed)
table(is.na(UC$endo_mayo))


# subsetting vst counts to patients with pred dose info
vst_counts_R_aza <- vst_counts_R_aza[,colnames(vst_counts_R_aza) %in% coldata_R_aza$sample_id] 

# excluding genes without names from sig genes
sig_genes_R_aza <- sig_genes_R_aza %>% 
  filter(!(genes == ""))

# Correlation
genes <- sig_genes_R_aza[1:30, 1]
all_correlations <- data.frame()

for (i in 1:30) {
  sig_gene <- genes[i]
  
  gene_counts <- vst_counts_R_aza[rownames(vst_counts_R_aza) == sig_gene,]
  coldata_gene <- as.data.frame(t(gene_counts))
  coldata_gene$sample_id <- rownames(coldata_gene)
  coldata_gene <- left_join(coldata_gene, coldata_R_aza, by = "sample_id")
  
  correlation_aza_dose <- cor.test(coldata_gene[,1], coldata_gene$aza_dose,  method = "spearman", exact=FALSE)
  correlation_hbi <- cor.test(coldata_gene[,1], coldata_gene$hbi, method = "spearman", exact=FALSE) # CD patients only
  correlation_pMayo <- cor.test(coldata_gene[,1], coldata_gene$partial_mayo, method = "spearman", exact=FALSE)# UC patients only
  correlation_age <- cor.test(coldata_gene[,1], coldata_gene$age,  method = "spearman", exact=FALSE)
  correlation_bmi <- cor.test(coldata_gene[,1], coldata_gene$bmi,  method = "spearman", exact=FALSE)
  correlation_crp <- cor.test(coldata_gene[,1], coldata_gene$crp,  method = "spearman", exact=FALSE)
  correlation_leuco <- cor.test(coldata_gene[,1], coldata_gene$leucocytes,  method = "spearman", exact=FALSE)
  correlation_erythro <- cor.test(coldata_gene[,1], coldata_gene$erythrocytes,  method = "spearman", exact=FALSE)
  correlation_thrombo <- cor.test(coldata_gene[,1], coldata_gene$thrombocytes,  method = "spearman", exact=FALSE)
  
  all_correlations <- rbind(all_correlations, data.frame(row.names = sig_gene, 
                                                         Aza_dose = correlation_aza_dose$estimate, P_value_aza_dose = correlation_aza_dose$p.value,
                                                         HBI = correlation_hbi$estimate, P_value_hbi = correlation_hbi$p.value,
                                                         pMayo = correlation_pMayo$estimate, P_value_pMayo = correlation_pMayo$p.value,
                                                         Age = correlation_age$estimate, P_value_age = correlation_age$p.value,
                                                         BMI = correlation_bmi$estimate, P_value_bmi = correlation_bmi$p.value,
                                                         CRP = correlation_crp$estimate, P_value_crp = correlation_crp$p.value,
                                                         Leukocytes = correlation_leuco$estimate, P_value_leuco = correlation_leuco$p.value,
                                                         Erythrocytes = correlation_erythro$estimate, P_value_erythro = correlation_erythro$p.value,
                                                         Thrombocytes = correlation_thrombo$estimate, P_value_thrombo = correlation_thrombo$p.value))
  
}

# Separating data frames with Rho and p-values
all_correlations_rho <- all_correlations %>% 
  dplyr::select(Aza_dose, HBI, pMayo, Age, BMI, CRP, Leukocytes, Erythrocytes, Thrombocytes)

all_correlations_p <- all_correlations %>% 
  dplyr::select(P_value_aza_dose, P_value_hbi,P_value_pMayo, P_value_age, P_value_bmi, P_value_crp, P_value_leuco, P_value_erythro, P_value_thrombo)

# Calculating adjusted p-values
all_correlations_p$Padj_aza_dose <- p.adjust(all_correlations_p$P_value_aza_dose, method = "BH")
all_correlations_p$Padj_hbi <- p.adjust(all_correlations_p$P_value_hbi, method = "BH")
all_correlations_p$Padj_pMayo <- p.adjust(all_correlations_p$P_value_pMayo, method = "BH")
all_correlations_p$Padj_age <- p.adjust(all_correlations_p$P_value_age, method = "BH")
all_correlations_p$Padj_bmi <- p.adjust(all_correlations_p$P_value_bmi, method = "BH")
all_correlations_p$Padj_crp <- p.adjust(all_correlations_p$P_value_crp, method = "BH")
all_correlations_p$Padj_leuco <- p.adjust(all_correlations_p$P_value_leuco, method = "BH")
all_correlations_p$Padj_erythro <- p.adjust(all_correlations_p$P_value_erythro, method = "BH")
all_correlations_p$Padj_thrombo <- p.adjust(all_correlations_p$P_value_thrombo, method = "BH")

any(all_correlations_p[,10:18] < 0.05) # No significant p-value after adjusting for multiple comparisons

# replacing gene symbols by gene names
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

rownames(all_correlations_rho) <- unique_ensg2gene[rownames(all_correlations_rho), ]$hgnc_symbol
rownames(all_correlations_p) <- unique_ensg2gene[rownames(all_correlations_p), ]$hgnc_symbol

# Heatmap
library(ComplexHeatmap)
library(RColorBrewer)

palette <- colorRampPalette(c("#440d57", "#20928c", "#efe51c"))
palette_hm <- palette(10)

heatmap_a <- Heatmap(t(all_correlations_rho), show_row_dend = FALSE, show_column_dend = FALSE, name = "Spearman's Rho", height = unit(5, "cm"), width = unit(16, "cm"),
                     col = palette_hm, row_names_side = "left",  column_names_side = "top", column_names_rot = 45, column_names_gp = gpar(fontsize = 10),
                     heatmap_legend_param = list(
                       title = "Spearman's Rho", at = c(-1, -0.8, 0, 0.8, 1),
                       title_position = "topcenter", legend_width = unit(7, "cm"),
                       legend_direction = "horizontal"
                     ))

draw(heatmap_a, heatmap_legend_side = "bottom")







########################################## Remitters + crp ###################################################

graphics.off()
rm(list = ls())

coldata_R_aza <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_aza_R_crp_04.05.txt", sep = "\t")
vst_counts_R_aza <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_aza_R_crp_04.05.txt", sep = "\t")
sig_genes_R_aza <-  read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_crp.txt", sep = "\t")
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting aza_dose, leukocytes, thrombocytes and erythrocytes from (updated) redcap coldata 
redcap_patients <- redcap[redcap$study_id %in% coldata_R_aza$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, aza_dose, leucocytes, erythrocytes, thrombocytes, endoscopy_performed, 
                partial_mayo, complete_mayo, hbi, cdai, sescd, endo_mayo)

idx <- match(coldata_R_aza$study_id, redcap_patients$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_aza <- left_join(coldata_R_aza, ord_redcap_patients, by = "study_id")

coldata_R_aza <- coldata_R_aza %>% 
  filter(Azathioprin == 1)

# Checking number of patients with info of clinical paramenters
CD <- coldata_R_aza %>% #14 CD
  filter(diagnosis_class == "CD") # 9 patients without sescd (no endoscopy); 5 had endoscopy
#table(CD$endoscopy_performed)

UC <- coldata_R_aza %>% #6 UC
  filter(diagnosis_class == "UC") # all have endomayo
# table(UC$endoscopy_performed)
# table(is.na(UC$endo_mayo))

# coldata_R_aza <- coldata_R_aza  %>% # excluding 1 CD without aza dose info (exclude)
#   filter(!(is.na(aza_dose)))


# subsetting vst counts to aza patients
vst_counts_R_aza <- vst_counts_R_aza[,colnames(vst_counts_R_aza) %in% coldata_R_aza$sample_id] 

# excluding genes without names from sig genes
sig_genes_R_aza <- sig_genes_R_aza %>% 
  filter(!(genes == ""))

# filtering to the genes in topvar+varpart
top25percent_varPart <- read.table("Output_files/VariancePartition/Remission_crp/aza/top25percent_varPart_aza_R_crp.txt")
top25percent_varPart_genes <- rownames(top25percent_varPart)

topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVar_genes <- rownames(topVar) 

intersect_topVar_varPart25 <- base::intersect(top25percent_varPart_genes, topVar_genes)

# Correlation
#genes <- sig_genes_R_aza[1:30, 1] #first 30 sig genes

genes <- intersect_topVar_varPart25
all_correlations <- data.frame()

for (i in 1:length(genes)) {
  sig_gene <- genes[i]
  
  gene_counts <- vst_counts_R_aza[rownames(vst_counts_R_aza) == sig_gene,]
  coldata_gene <- as.data.frame(t(gene_counts))
  coldata_gene$sample_id <- rownames(coldata_gene)
  coldata_gene <- left_join(coldata_gene, coldata_R_aza, by = "sample_id")
  
  correlation_aza_dose <- cor.test(coldata_gene[,1], coldata_gene$aza_dose,  method = "spearman", exact=FALSE)
  #correlation_hbi <- cor.test(coldata_gene[,1], coldata_gene$hbi, method = "spearman", exact=FALSE) # CD patients only
  #correlation_pMayo <- cor.test(coldata_gene[,1], coldata_gene$partial_mayo, method = "spearman", exact=FALSE)# UC patients only
  correlation_age <- cor.test(coldata_gene[,1], coldata_gene$age,  method = "spearman", exact=FALSE)
  correlation_bmi <- cor.test(coldata_gene[,1], coldata_gene$bmi,  method = "spearman", exact=FALSE)
  correlation_crp <- cor.test(coldata_gene[,1], coldata_gene$crp,  method = "spearman", exact=FALSE)
  correlation_leuco <- cor.test(coldata_gene[,1], coldata_gene$leucocytes,  method = "spearman", exact=FALSE)
  correlation_erythro <- cor.test(coldata_gene[,1], coldata_gene$erythrocytes,  method = "spearman", exact=FALSE)
  correlation_thrombo <- cor.test(coldata_gene[,1], coldata_gene$thrombocytes,  method = "spearman", exact=FALSE)
  
  all_correlations <- rbind(all_correlations, data.frame(row.names = sig_gene, 
                                                         Aza_dose = correlation_aza_dose$estimate, P_value_aza_dose = correlation_aza_dose$p.value,
                                                         #HBI = correlation_hbi$estimate, P_value_hbi = correlation_hbi$p.value,
                                                         #pMayo = correlation_pMayo$estimate, P_value_pMayo = correlation_pMayo$p.value,
                                                         Age = correlation_age$estimate, P_value_age = correlation_age$p.value,
                                                         BMI = correlation_bmi$estimate, P_value_bmi = correlation_bmi$p.value,
                                                         CRP = correlation_crp$estimate, P_value_crp = correlation_crp$p.value,
                                                         Leukocytes = correlation_leuco$estimate, P_value_leuco = correlation_leuco$p.value,
                                                         Erythrocytes = correlation_erythro$estimate, P_value_erythro = correlation_erythro$p.value,
                                                         Thrombocytes = correlation_thrombo$estimate, P_value_thrombo = correlation_thrombo$p.value))
  
}

# Separating data frames with Rho and p-values
all_correlations_rho <- all_correlations %>% 
  dplyr::select(Aza_dose, Age, BMI, CRP, Leukocytes, Erythrocytes, Thrombocytes) # HBI, pMayo,

all_correlations_p <- all_correlations %>% 
  dplyr::select(P_value_aza_dose,  P_value_age, P_value_bmi, P_value_crp, P_value_leuco, P_value_erythro, P_value_thrombo)#P_value_hbi,P_value_pMayo,

# Calculating adjusted p-values
all_correlations_padj <- all_correlations_p %>% 
  as.matrix %>% 
  as.vector %>% 
  p.adjust(method='BH') %>% 
  matrix(ncol=7)
all_correlations_padj <- as.data.frame(all_correlations_padj)
rownames(all_correlations_padj) <- rownames(all_correlations_p)
colnames(all_correlations_padj) <- colnames(all_correlations_p)

which(any(all_correlations_padj < 0.05)) #no significant


# replacing gene symbols by gene names
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

all_correlations_rho$gene <- unique_ensg2gene[rownames(all_correlations_rho), ]$hgnc_symbol
all_correlations_rho <- all_correlations_rho[!(is.na(all_correlations_rho$gene) | all_correlations_rho$gene == ""), ]
rownames(all_correlations_rho) <- all_correlations_rho$gene
all_correlations_rho <- all_correlations_rho %>% #excluding gene names column
  dplyr::select(!("gene"))

all_correlations_padj$gene <- unique_ensg2gene[rownames(all_correlations_padj), ]$hgnc_symbol
all_correlations_padj <- all_correlations_padj[!(is.na(all_correlations_padj$gene) | all_correlations_padj$gene == ""), ]
rownames(all_correlations_padj) <- all_correlations_padj$gene
all_correlations_padj <- all_correlations_padj %>% #excluding gene names column
  dplyr::select(!("gene"))

# Heatmap
library(ComplexHeatmap)
library(RColorBrewer)

palette <- colorRampPalette(c("#440d57", "#20928c", "#efe51c"))
palette_hm <- palette(10)

heatmap_a <- Heatmap(t(all_correlations_rho), show_row_dend = FALSE, show_column_dend = FALSE, name = "Spearman's Rho", height = unit(5, "cm"), width = unit(25, "cm"),
                     col = palette_hm, row_names_side = "left",  column_names_side = "top", column_names_rot = 45, column_names_gp = gpar(fontsize = 8),
                     heatmap_legend_param = list(
                       title = "Spearman's Rho", at = c(-1, -0.8, 0, 0.8, 1),
                       title_position = "topcenter", legend_width = unit(7, "cm"),
                       legend_direction = "horizontal"
                     ))

draw(heatmap_a, heatmap_legend_side = "bottom")



# reordering to match heatmap


order <- c("BOK", "CLIC3", "MLC1", "SPON2", "PRF1", "NMUR1","KIR2DL3","KIR2DL1","KIR3DL1","SLC1A7","GNLY","CLDND2","NKG7","FGFBP2","GZMB",   
           "PRSS23", "ADGRG1","TBX21","S1PR5","C1orf21","FCRL6" ,   "CCL4"  ,   "XCL2"    , "FASLG"  ,  "CD160"   , "GZMA" ,    "KLRF1"  ,  
           "KLRD1"  ,  "HOPX"  ,   "NCAM1", "SH2D1B",   "RNF165" ,  "DGKK",     "AKR1C3",   "GRIK4",    "TSPEAR",   "CTSE",     "GADD45A",
           "HEPACAM2", "TMCC2",    "RHAG",     "CR1L",     "C2orf88",  "RPL3L" ,   "SHISA4",  "GDF15" ,   "ITLN1" ,   "KEL",      "AQP1",
           "KCNH2"  ,  "PGF"   ,   "YPEL4"  ,  "DYRK3"  ,  "KRT79"  ,  "ETV7"  ,   "KIR2DL4" , "TCL1A"  ,  "CD72"   ,  "CORO2B")

all_correlations_rho <- all_correlations_rho[order(match(rownames(all_correlations_rho), order)), , drop = FALSE]
all_correlations_padj <- all_correlations_padj[order(match(rownames(all_correlations_padj), order)), , drop = FALSE]


# Inserting * on significant correlations
padj_all_correlations <- t(as.matrix(all_correlations_padj))
all_correlations_rho <- t(as.matrix(all_correlations_rho))

heatmap_significant <- Heatmap(all_correlations_rho, 
                               cell_fun = function(j, i, x, y, w, h, fill) {
                                 if(padj_all_correlations[i, j] < 0.001) {
                                   grid.text("***", x, y)
                                 } else if(padj_all_correlations[i, j] < 0.01) {
                                   grid.text("**", x, y)
                                 } else if(padj_all_correlations[i, j] < 0.05) {
                                   grid.text("*", x, y)
                                 }
                               },
                               show_row_dend = FALSE, show_column_dend = FALSE, name = "Spearman's Rho", height = unit(5, "cm"), width = unit(30, "cm"),
                               col = palette_hm, row_names_side = "left",  column_names_side = "top", column_names_rot = 45, column_names_gp = gpar(fontsize = 10),
                               column_order= order,
                               heatmap_legend_param = list(
                                 title = "Spearman's Rho", at = c(-1, -0.8, 0, 0.8, 1),
                                 title_position = "topcenter", legend_width = unit(7, "cm"),
                                 legend_direction = "horizontal"
                               ))

draw(heatmap_significant, heatmap_legend_side = "bottom")


# vertical
padj_all_correlations <- as.matrix(all_correlations_padj)
all_correlations_rho <- as.matrix(all_correlations_rho)


heatmap_significant <- Heatmap(t(all_correlations_rho), 
                               cell_fun = function(j, i, x, y, w, h, fill) {
                                 if(padj_all_correlations[i, j] < 0.001) {
                                   grid.text("***", x, y)
                                 } else if(padj_all_correlations[i, j] < 0.01) {
                                   grid.text("**", x, y)
                                 } else if(padj_all_correlations[i, j] < 0.05) {
                                   grid.text("*", x, y)
                                 }
                               },
                               show_row_dend = FALSE, show_column_dend = FALSE, name = "Spearman's Rho", height = unit(30, "cm"), width = unit(7, "cm"),
                               col = palette_hm, row_names_side = "left",  column_names_side = "top", column_names_rot = 45, #column_names_gp = gpar(fontsize = 10),
                               row_order = order,
                               heatmap_legend_param = list(
                                 title = "Spearman's Rho", at = c(-1, -0.8, 0, 0.8, 1),
                                 title_position = "topcenter", legend_width = unit(7, "cm"),
                                 legend_direction = "horizontal"
                               ))

draw(heatmap_significant, heatmap_legend_side = "bottom")






#padj ----
# # wrong
# p_aza_dose <- all_correlations_p$P_value_aza_dose
# p_hbi <- all_correlations_p$P_value_hbi
# p_Mayo <- all_correlations_p$P_value_pMayo
# p_age <- all_correlations_p$P_value_age
# p_bmi <- all_correlations_p$P_value_bmi
# p_crp <- all_correlations_p$P_value_crp
# p_leuco <- all_correlations_p$P_value_leuco
# p_erythro <- all_correlations_p$P_value_erythro
# p_thrombo <- all_correlations_p$P_value_thrombo
# 
# all_p <- c(p_aza_dose, p_hbi, p_Mayo, p_age, p_bmi, p_crp, p_leuco, p_erythro, p_thrombo)
# all_correlations_padj <- p.adjust(all_p, method = "BH")
# 
# 
# # before:
# all_correlations_p$Padj_aza_dose <- p.adjust(all_correlations_p$P_value_aza_dose, method = "BH")
# #all_correlations_p$Padj_hbi <- p.adjust(all_correlations_p$P_value_hbi, method = "BH")
# #all_correlations_p$Padj_pMayo <- p.adjust(all_correlations_p$P_value_pMayo, method = "BH")
# all_correlations_p$Padj_age <- p.adjust(all_correlations_p$P_value_age, method = "BH")
# all_correlations_p$Padj_bmi <- p.adjust(all_correlations_p$P_value_bmi, method = "BH")
# all_correlations_p$Padj_crp <- p.adjust(all_correlations_p$P_value_crp, method = "BH")
# all_correlations_p$Padj_leuco <- p.adjust(all_correlations_p$P_value_leuco, method = "BH")
# all_correlations_p$Padj_erythro <- p.adjust(all_correlations_p$P_value_erythro, method = "BH")
# all_correlations_p$Padj_thrombo <- p.adjust(all_correlations_p$P_value_thrombo, method = "BH")
# 
# which(any(all_correlations_p[,10:18] < 0.05))
# all_correlations_p[all_correlations_p$Padj_aza_dose < 0.05,] #PGF and GDF15




male <- coldata_R_aza %>% 
  filter(sex == "Male")

mean(male$thrombocytes)
sd(male$thrombocytes)

mean(male$erythrocytes)
sd(male$erythrocytes)

mean(male$leucocytes)
sd(male$leucocytes)

female <- coldata_R_aza %>% 
  filter(sex == "Female")

mean(female$thrombocytes)
sd(female$thrombocytes)

mean(female$erythrocytes)
sd(female$erythrocytes)

mean(female$leucocytes)
sd(female$leucocytes)

write.table(male, "Output_files/Lab_parameters/Males_aza.txt", sep = "\t")
write.table(female, "Output_files/Lab_parameters/Females_aza.txt", sep = "\t")
