# EZE cohort 
# Correlation clin parameters (patients with pred dose info) and gene expression
  # Date: 08.05


graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

install.packages("ggpubr")
library("ggpubr")
library(tidyverse)


# Remitters
coldata_R_pred <- read.table("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_pred_vs_noSyst_R_17.03.txt", sep = "\t")
vst_counts_R_pred <- read.table("Output_files/DESeq2/Maaslin2/Remission/DESeq2_vst_counts_pred_vs_noSyst_R_29.03.txt", sep = "\t")
sig_genes_R_pred <- read.csv("Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_vst_q0.05.txt",
                      sep = "\t")
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# #Filtering coldata to have only patients with pred dose information
# coldata_R_pred <- coldata_R_pred %>% 
#   filter(!is.na(prednisolone_dose))

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_pred$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)

idx <- match(coldata_R_pred$study_id, redcap_patients$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_pred <- left_join(coldata_R_pred, ord_redcap_patients, by = "study_id")

# subsetting vst counts to patients with pred dose info
vst_counts_R_pred <- vst_counts_R_pred[,colnames(vst_counts_R_pred) %in% coldata_R_pred$sample_id] #do I need to vst re-transform?

# excluding genes without names from sig genes
sig_genes_R_pred <- sig_genes_R_pred %>% 
  filter(!(genes == ""))

# Correlation
genes <- sig_genes_R_pred[1:30, 1]
all_correlations <- data.frame()

for (i in 1:30) {
 sig_gene <- genes[i]
  
  gene_counts <- vst_counts_R_pred[rownames(vst_counts_R_pred) == sig_gene,]
  coldata_gene <- as.data.frame(t(gene_counts))
  coldata_gene$sample_id <- rownames(coldata_gene)
  coldata_gene <- left_join(coldata_gene, coldata_R_pred, by = "sample_id")
  
  correlation_pred <- cor.test(coldata_gene[,1], coldata_gene$prednisolone_dose,  method = "spearman", exact=FALSE)
  correlation_age <- cor.test(coldata_gene[,1], coldata_gene$age,  method = "spearman", exact=FALSE)
  correlation_bmi <- cor.test(coldata_gene[,1], coldata_gene$bmi,  method = "spearman", exact=FALSE)
  correlation_crp <- cor.test(coldata_gene[,1], coldata_gene$crp,  method = "spearman", exact=FALSE)
  correlation_leuco <- cor.test(coldata_gene[,1], coldata_gene$leucocytes,  method = "spearman", exact=FALSE)
  correlation_erythro <- cor.test(coldata_gene[,1], coldata_gene$erythrocytes,  method = "spearman", exact=FALSE)
  correlation_thrombo <- cor.test(coldata_gene[,1], coldata_gene$thrombocytes,  method = "spearman", exact=FALSE)
  
  all_correlations <- rbind(all_correlations, data.frame(row.names = sig_gene, 
                                                         Pred_dose = correlation_pred$estimate, P_value_pred = correlation_pred$p.value,
                                                         Age = correlation_age$estimate, P_value_age = correlation_age$p.value,
                                                         BMI = correlation_bmi$estimate, P_value_bmi = correlation_bmi$p.value,
                                                         CRP = correlation_crp$estimate, P_value_crp = correlation_crp$p.value,
                                                         Leukocytes = correlation_leuco$estimate, P_value_leuco = correlation_leuco$p.value,
                                                         Erythrocytes = correlation_erythro$estimate, P_value_erythro = correlation_erythro$p.value,
                                                         Thrombocytes = correlation_thrombo$estimate, P_value_thrombo = correlation_thrombo$p.value))
  
}


# Separating data frames with Rho and p-values
all_correlations_rho <- all_correlations %>% 
  dplyr::select(Pred_dose, Age, BMI, CRP, Leukocytes, Erythrocytes, Thrombocytes)

all_correlations_p <- all_correlations %>% 
  dplyr::select(P_value_pred, P_value_age, P_value_bmi, P_value_crp, P_value_leuco, P_value_erythro, P_value_thrombo)

# Calculating adjusted p-values
all_correlations_p$Padj_value_pred <- p.adjust(all_correlations_p$P_value_pred, method = "BH")
all_correlations_p$Padj_value_age <- p.adjust(all_correlations_p$P_value_age, method = "BH")
all_correlations_p$Padj_value_bmi <- p.adjust(all_correlations_p$P_value_bmi, method = "BH")
all_correlations_p$Padj_value_crp <- p.adjust(all_correlations_p$P_value_crp, method = "BH")
all_correlations_p$Padj_value_leuco <- p.adjust(all_correlations_p$P_value_leuco, method = "BH")
all_correlations_p$Padj_value_erythro <- p.adjust(all_correlations_p$P_value_erythro, method = "BH")
all_correlations_p$Padj_value_thrombo <- p.adjust(all_correlations_p$P_value_thrombo, method = "BH")

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


any(all_correlations_p[,7:11] < 0.05) # No significant p-value after adjusting for multiple comparisons









############ Remitters + crp ##############

graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

install.packages("ggpubr")
library("ggpubr")
library(tidyverse)


# Remitters
coldata_R_pred <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_pred_R_crp_04.05.txt", sep = "\t")
vst_counts_R_pred <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_vst_pred_R_crp_04.05.txt", sep = "\t")
sig_genes_R_pred <- read.table( "Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt", sep = "\t")
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

sig_genes_topVar_lfc0.5 <- read.table("Output_files/Heatmap/vst_sig_genes_maaslin2_remission_crp/pred/Pred_top50_siggenes_lfc0.5_topVar_HM.txt", sep = "\t")

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_pred$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes, neutrophils)

idx <- match(coldata_R_pred$study_id, redcap_patients$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_pred <- left_join(coldata_R_pred, ord_redcap_patients, by = "study_id")

# subsetting vst counts to patients with pred dose info
vst_counts_R_pred <- vst_counts_R_pred[,colnames(vst_counts_R_pred) %in% coldata_R_pred$sample_id] 

# excluding genes without names from sig genes
sig_genes_R_pred <- sig_genes_R_pred %>% 
  filter(!(genes == ""))

# Correlation
genes <- sig_genes_R_pred[sig_genes_R_pred$feature %in% sig_genes_topVar_lfc0.5$feature,1]
#genes <- sig_genes_R_pred[1:30, 1]
all_correlations <- data.frame()

for (i in 1:length(genes)) {
  sig_gene <- genes[i]
  
  gene_counts <- vst_counts_R_pred[rownames(vst_counts_R_pred) == sig_gene,]
  coldata_gene <- as.data.frame(t(gene_counts))
  coldata_gene$sample_id <- rownames(coldata_gene)
  coldata_gene <- left_join(coldata_gene, coldata_R_pred, by = "sample_id")
  
  correlation_pred <- cor.test(coldata_gene[,1], coldata_gene$prednisolone_dose,  method = "spearman", exact=FALSE)
  correlation_age <- cor.test(coldata_gene[,1], coldata_gene$age,  method = "spearman", exact=FALSE)
  correlation_bmi <- cor.test(coldata_gene[,1], coldata_gene$bmi,  method = "spearman", exact=FALSE)
  correlation_crp <- cor.test(coldata_gene[,1], coldata_gene$crp,  method = "spearman", exact=FALSE)
  correlation_leuco <- cor.test(coldata_gene[,1], coldata_gene$leucocytes,  method = "spearman", exact=FALSE)
  correlation_erythro <- cor.test(coldata_gene[,1], coldata_gene$erythrocytes,  method = "spearman", exact=FALSE)
  correlation_thrombo <- cor.test(coldata_gene[,1], coldata_gene$thrombocytes,  method = "spearman", exact=FALSE)
  correlation_neutrophils <- cor.test(coldata_gene[,1], coldata_gene$neutrophils,  method = "spearman", exact=FALSE)
  
  
  all_correlations <- rbind(all_correlations, data.frame(row.names = sig_gene, 
                                                         Pred_dose = correlation_pred$estimate, P_value_pred = correlation_pred$p.value,
                                                         Age = correlation_age$estimate, P_value_age = correlation_age$p.value,
                                                         BMI = correlation_bmi$estimate, P_value_bmi = correlation_bmi$p.value,
                                                         CRP = correlation_crp$estimate, P_value_crp = correlation_crp$p.value,
                                                         Leukocytes = correlation_leuco$estimate, P_value_leuco = correlation_leuco$p.value,
                                                         Erythrocytes = correlation_erythro$estimate, P_value_erythro = correlation_erythro$p.value,
                                                         Thrombocytes = correlation_thrombo$estimate, P_value_thrombo = correlation_thrombo$p.value,
                                                         Neutrophils = correlation_neutrophils$estimate, P_value_neutrophils = correlation_neutrophils$p.value))
  
  
}


# Separating data frames with Rho and p-values
all_correlations_rho <- all_correlations %>% 
  dplyr::select(Pred_dose, Age, BMI, CRP, Leukocytes, Erythrocytes, Thrombocytes, Neutrophils)

all_correlations_p <- all_correlations %>% 
  dplyr::select(P_value_pred, P_value_age, P_value_bmi, P_value_crp, P_value_leuco, P_value_erythro, P_value_thrombo, P_value_neutrophils)

# Calculating adjusted p-values
all_correlations_padj <- all_correlations_p %>% 
  as.matrix %>% 
  as.vector %>% 
  p.adjust(method='BH') %>% 
  matrix(ncol=8)

all_correlations_padj <- as.data.frame(all_correlations_padj)
rownames(all_correlations_padj) <- rownames(all_correlations_p)
colnames(all_correlations_padj) <- colnames(all_correlations_p)

# which(all_correlations_p < 0.05)
# which(all_correlations_padj < 0.05)

# replacing gene symbols by gene names
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

rownames(all_correlations_rho) <- unique_ensg2gene[rownames(all_correlations_rho), ]$hgnc_symbol
rownames(all_correlations_padj) <- unique_ensg2gene[rownames(all_correlations_padj), ]$hgnc_symbol

# Heatmap
library(ComplexHeatmap)
library(RColorBrewer)

palette <- colorRampPalette(c("#440d57", "#20928c", "#efe51c"))
palette_hm <- palette(10)

heatmap_a <- Heatmap(t(as.matrix(all_correlations_rho)), show_row_dend = FALSE, show_column_dend = FALSE, name = "Spearman's Rho", height = unit(5, "cm"), width = unit(16, "cm"),
                     col = palette_hm, row_names_side = "left",  column_names_side = "top", column_names_rot = 45, column_names_gp = gpar(fontsize = 10),
                     heatmap_legend_param = list(
                       title = "Spearman's Rho", at = c(-1, -0.8, 0, 0.8, 1),
                       title_position = "topcenter", legend_width = unit(7, "cm"),
                       legend_direction = "horizontal"
                     ))

draw(heatmap_a, heatmap_legend_side = "bottom")



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
show_row_dend = FALSE, show_column_dend = FALSE, name = "Spearman's Rho", height = unit(5, "cm"), width = unit(16, "cm"),
col = palette_hm, row_names_side = "left",  column_names_side = "top", column_names_rot = 45, column_names_gp = gpar(fontsize = 10),
heatmap_legend_param = list(
  title = "Spearman's Rho", at = c(-1, -0.8, 0, 0.8, 1),
  title_position = "topcenter", legend_width = unit(7, "cm"),
  legend_direction = "horizontal"
))

draw(heatmap_significant, heatmap_legend_side = "bottom")


### Reordering to match heatmap order

order <- c("DBH-AS1",   "LIMS2",     "ID3",       "SPIB" ,     "IRS2"  ,    "MIR223HG" , "CLEC4E"  ,  "TLR2" ,     "FKBP5",
           "KCNE1"  ,   "ECHDC3"  ,  "TPST1"  ,   "INHBB"   ,  "NSUN7"   , "DUSP1"  ,   "TSC22D3",   "GADD45A" ,  "ASPH" ,     
           "IRAK3",     "ZNF608" ,   "VNN1" ,     "GPR141" ,   "SIPA1L2" ,  "IL18R1" ,   "PFKFB2" ,   "SLC8A1" ,   "COL9A2"  ,  "VSIG4" ,   
           "DAAM2"  ,   "MAOA"  ,    "FLT3"  ,    "AMPH"  ,    "IL1R2"  ,   "ADAMTS2" ,  "OLAH"  ,    "GRB10"  ,   "PCSK9" ,    "ARG1",
           "LINC02207", "ABCG1"  ,   "PPIAP29" ,  "SLC5A9"  ,  "DUSP13"  ,  "IFITM3P2" , "CD163"  ,   "GPER1"  ,   "C5orf67"  , "SLC22A4"  , "TLR5"   ,   "NAIP")

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







male <- coldata_R_pred %>% 
  filter(sex == "Male")

mean(male$thrombocytes)
sd(male$thrombocytes)

mean(male$erythrocytes)
sd(male$erythrocytes)

mean(male$leucocytes)
sd(male$leucocytes)



female <- coldata_R_pred %>% 
  filter(sex == "Female") %>% 
  filter(!(is.na(thrombocytes)))

mean(female$thrombocytes)
sd(female$thrombocytes)

mean(female$erythrocytes)
sd(female$erythrocytes)

mean(female$leucocytes)
sd(female$leucocytes)



write.table(male, "Output_files/Lab_parameters/Males_pred.txt", sep = "\t")
write.table(female, "Output_files/Lab_parameters/Females_pred.txt", sep = "\t")
