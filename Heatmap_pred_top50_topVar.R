# DE Analysis - EZE cohort
# Heatmaps - inactive disease- Aza: topvar + varPart 25%

graphics.off()
rm(list = ls())

setwd("C:\\Documents\\Masters thesis\\EZE_cohort") #laptop
setwd("D:\\Documentos\\Workspace\\Masters-Thesis\\EZE\\EZE_cohort") #PC

# Loading packages ----
library(ComplexHeatmap)
library(tidyverse)
library(RColorBrewer)
library(circlize)

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = T, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Loading data ----
topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVar_genes <- rownames(topVar)

coldata_pred <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_heatmap_pred_healthy_R_crp.txt", sep = "\t")
normalized_counts <- read.table("Output_files/DESeq2/Heatmap/Remission_crp/DESeq2_normalized_pred_healthy_R_crp.txt", sep = "\t")

sig_genes_R_pred <- read.table( "Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt", sep = "\t")

sig_genes_lfc0.5 <- sig_genes_R_pred %>% 
  filter(abs(coef) > 0.5)

sig_genes_TVG <- intersect(topVar_genes,sig_genes_lfc0.5$feature )

sig_genes_lfc0.5_TVG <- sig_genes_lfc0.5[sig_genes_lfc0.5$feature %in% sig_genes_TVG,]
sig_genes_lfc0.5_TVG <- sig_genes_lfc0.5_TVG %>% 
  arrange(qval)
sig_genes_lfc0.5_TVG <- sig_genes_lfc0.5_TVG[!(is.na(sig_genes_lfc0.5_TVG$genes) | sig_genes_lfc0.5_TVG$genes == ""), ]

top50_sig_genes_lfc0.5 <- sig_genes_lfc0.5_TVG[1:50,]

write.table(top50_sig_genes_lfc0.5, "Output_files/Heatmap/vst_sig_genes_maaslin2_remission_crp/pred/Pred_top50_siggenes_topVar_HM.txt", sep = "\t")
# Counts preprocessing ----
heatmap_ed <- normalized_counts[top50_sig_genes_lfc0.5$feature, ]
heatmap_ed_scaled <- as.data.frame(t(base::scale(t(heatmap_ed))))

heatmap_ed_scaled$gene <- unique_ensg2gene[rownames(heatmap_ed_scaled), ]$hgnc_symbol
rownames(heatmap_ed_scaled) <- heatmap_ed_scaled$gene
heatmap_ed_scaled <- heatmap_ed_scaled %>% #excluding gene names column
  dplyr::select(!("gene"))

# Coldata preprocessing to be plotted ----
coldata_pred <- coldata_pred %>%
  mutate(Syst_therapy = recode(No_syst,
                               "0" = "Yes",
                               "1" = "No")) %>% 
  mutate(across(c(diagnosis_class, sex, remission, Syst_therapy, pred_vs_noSyst), as.factor))



# Expression data matrix preprocessing ----
m.top50_ed_R_pred <- data.matrix(heatmap_ed_scaled, rownames.force = NA)



# Heatmap annotations ----
colors_hm <- colorRamp2(c(-4, 0, 4), c("#333b93", "white", "#b61728"))

col_ha_pred = HeatmapAnnotation(df = data.frame(
  Sex = coldata_pred$sex,
  BMI = coldata_pred$bmi,
  Age = coldata_pred$age,
  #systemic_therapy = coldata_R_pred$Syst_therapy,
  Diagnosis = coldata_pred$diagnosis_class,
  Pred_dose = coldata_pred$prednisolone_dose,
  Biologics = coldata_pred$biologics),
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
                           "no_biologics" = "#f5c710")
  ),
  CRP = anno_barplot(coldata_pred$crp)#,
  #annotation_name_gp= gpar(fontsize = 12, fontface = "bold")
)

# Heatmap
coldata_pred$pred_vs_noSyst <- ordered(coldata_pred$pred_vs_noSyst, levels = c("Healthy", "no_syst","pred"))

hm = ComplexHeatmap::Heatmap(m.top50_ed_R_pred, top_annotation = col_ha_pred, show_column_names = FALSE, name = "Z score", #height = unit(30, "cm"), width = unit(15, "cm"),
                             show_row_names = TRUE, border = TRUE, column_split = coldata_pred$pred_vs_noSyst, column_gap = unit(2, "mm"),
                             cluster_row_slices = T,cluster_column_slices = FALSE, row_title_rot = 0,
                             show_column_dend = FALSE, show_row_dend = TRUE, #row_names_gp = gpar(fontsize=7), 
                             col = colors_hm, heatmap_legend_param = list(
                               title = "Z score", 
                               legend_direction = "horizontal"
                             ))

hm = draw(hm, heatmap_legend_side = "bottom")

index <- row_order(hm)
gene_order <- rownames(m.top50_ed_R_pred)[index]  
