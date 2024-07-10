# EZE cohort: Therapy Signatures
  # Heatmap: azathioprine x no systemic therapies
  # DEGs, TVGs, variance partition > 25% 
  # DEGs shared with TVGs and not shared with TVGs (unique)
  # Figure 15A
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(ComplexHeatmap)
library(tidyverse)
library(colorRamp2)
library(circlize)

# Loading files ----
# gene names
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = T, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# TVG
topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")

# heatmap coldata and counts
coldata_aza <- read.table("Cleaned_tables/heatmaps/hm_coldata_aza.txt", sep = "\t")
normalized_counts <- read.table("Cleaned_tables/heatmaps/hm_normalized_counts_aza.txt", sep = "\t")

# variance > 25% (varPart)
top25percent_varPart <- read.table("Output_files/Variance_partition/aza/top25percent_varPart_aza.txt")

# Shared genes (TVGs, DEGs, varpart>25%)----
topVar_genes <- rownames(topVar)
top25percent_varPart_genes <- rownames(top25percent_varPart)
intersect_topVar_varPart25 <- base::intersect(top25percent_varPart_genes, topVar_genes)

# Coldata preprocessing ----
coldata_aza <- coldata_aza %>%
  mutate(Syst_therapy = recode(No_syst,
                               "0" = "Yes",
                               "1" = "No")) %>% 
  mutate(Prednisolon = recode(Prednisolon,
                              "0" = "No",
                              "1" = "Yes")) %>% 
  mutate(across(c(diagnosis_class, sex, remission, Syst_therapy, aza_vs_noSyst), as.factor))

coldata_aza$aza_vs_noSyst <- ordered(coldata_aza$aza_vs_noSyst, levels = c("Healthy","no_syst","aza"))

# Counts preprocessing ----
heatmap_ed <- normalized_counts[intersect_topVar_varPart25, ]
heatmap_ed_scaled <- as.data.frame(t(base::scale(t(heatmap_ed))))

heatmap_ed_scaled$gene <- unique_ensg2gene[rownames(heatmap_ed_scaled), ]$hgnc_symbol
heatmap_ed_scaled <- heatmap_ed_scaled[!(is.na(heatmap_ed_scaled$gene) | heatmap_ed_scaled$gene == ""), ]
rownames(heatmap_ed_scaled) <- heatmap_ed_scaled$gene
heatmap_ed_scaled <- heatmap_ed_scaled %>%
  dplyr::select(!("gene"))

m.heatmap_ed_scaled <- data.matrix(heatmap_ed_scaled, rownames.force = NA)

# Heatmap plotting ----
colors_hm <- colorRamp2(c(-4, 0, 4), c("#333b93", "white", "#b61728"))

col_ha_aza = HeatmapAnnotation(df = data.frame(
  Sex = coldata_aza$sex,
  BMI = coldata_aza$bmi,
  Age = coldata_aza$age,
  Diagnosis = coldata_aza$diagnosis_class,
  Biologics = coldata_aza$biologics),
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
  CRP = anno_barplot(coldata_aza$crp)
  )

pdf(file = "Output_files/Heatmaps/aza/heatmap_aza_healthy.pdf", width = 10, height= 15)
hm = ComplexHeatmap::Heatmap(m.heatmap_ed_scaled, top_annotation = col_ha_aza, show_column_names = FALSE, name = "Z score", 
                             height = unit(30, "cm"), width = unit(15, "cm"),show_row_names = TRUE, border = TRUE, 
                             column_split = coldata_aza$aza_vs_noSyst, column_gap = unit(2, "mm"),
                             cluster_row_slices = TRUE, cluster_column_slices = FALSE, row_title_rot = 0,
                             show_column_dend = FALSE, show_row_dend = TRUE,
                             col = colors_hm, heatmap_legend_param = list(title = "Z score", legend_direction = "horizontal"))
hm = draw(hm, heatmap_legend_side = "bottom")
dev.off()
