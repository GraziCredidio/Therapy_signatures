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
coldata_R_aza <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_heatmap_aza_healthy_R_crp.txt", sep = "\t")
normalized_counts <- read.table("Output_files/DESeq2/Heatmap/Remission_crp/DESeq2_normalized_aza_healthy_R_crp.txt", sep = "\t")

#sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_crp.txt", sep = "\t")

top25percent_varPart <- read.table("Output_files/VariancePartition/Remission_crp/aza/top25percent_varPart_aza_R_crp.txt")
top25percent_varPart_genes <- rownames(top25percent_varPart)

topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVar_genes <- rownames(topVar) 

intersect_topVar_varPart25 <- base::intersect(top25percent_varPart_genes, topVar_genes)

# Counts preprocessing ----
heatmap_ed <- normalized_counts[intersect_topVar_varPart25, ]
heatmap_ed_scaled <- as.data.frame(t(base::scale(t(heatmap_ed))))

heatmap_ed_scaled$gene <- unique_ensg2gene[rownames(heatmap_ed_scaled), ]$hgnc_symbol
heatmap_ed_scaled <- heatmap_ed_scaled[!(is.na(heatmap_ed_scaled$gene) | heatmap_ed_scaled$gene == ""), ]
rownames(heatmap_ed_scaled) <- heatmap_ed_scaled$gene
heatmap_ed_scaled <- heatmap_ed_scaled %>% #excluding gene names column
  dplyr::select(!("gene"))

# Coldata preprocessing to be plotted ----
coldata_R_aza <- coldata_R_aza %>%
  mutate(Syst_therapy = recode(No_syst,
                               "0" = "Yes",
                               "1" = "No")) %>% 
  mutate(Prednisolon = recode(Prednisolon,
                              "0" = "No",
                              "1" = "Yes")) %>% 
  mutate(across(c(diagnosis_class, sex, remission, Syst_therapy, aza_vs_noSyst), as.factor))



# Expression data matrix preprocessing ----
m.top100_ed_R_aza <- data.matrix(heatmap_ed_scaled, rownames.force = NA)

# Heatmap annotations ----
colors_hm <- colorRamp2(c(-4, 0, 4), c("#333b93", "white", "#b61728"))

col_ha_aza = HeatmapAnnotation(df = data.frame(
  Sex = coldata_R_aza$sex,
  BMI = coldata_R_aza$bmi,
  Age = coldata_R_aza$age,
  #Systemic_therapy = coldata_R_aza$Syst_therapy,
  Diagnosis = coldata_R_aza$diagnosis_class,
  #Pred_dose = coldata_R_aza$prednisolone_dose,
  #Pred = coldata_R_aza$Prednisolon,
  Biologics = coldata_R_aza$biologics),
  col = list(Sex = c("Female" = "#b492c4",
                     "Male" = "#fbaa51"),
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
  CRP = anno_barplot(coldata_R_aza$crp)#,
  #annotation_name_gp= gpar(fontsize = 12, fontface = "bold")
)


# Heatmap
coldata_R_aza$aza_vs_noSyst <- ordered(coldata_R_aza$aza_vs_noSyst, levels = c("Healthy","no_syst","aza"))
hm = ComplexHeatmap::Heatmap(m.top100_ed_R_aza, top_annotation = col_ha_aza, show_column_names = FALSE, name = "Z score", height = unit(30, "cm"), width = unit(15, "cm"),
                        show_row_names = TRUE, border = TRUE, column_split = coldata_R_aza$aza_vs_noSyst, column_gap = unit(2, "mm"),
                        cluster_row_slices = T,cluster_column_slices = FALSE, row_title_rot = 0,
                        show_column_dend = FALSE, show_row_dend = TRUE, #row_names_gp = gpar(fontsize=7), 
                        col = colors_hm, heatmap_legend_param = list(
                          title = "Z score", 
                          legend_direction = "horizontal"
                        ))

hm = draw(hm, heatmap_legend_side = "bottom")
index <- row_order(hm)
gene_order <- rownames(m.top100_ed_R_aza)[index]  
