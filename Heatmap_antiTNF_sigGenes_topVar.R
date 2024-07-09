# DE Analysis - EZE cohort
# Heatmap - inactive disease - Anti-TNF: sig + topVar

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

coldata_noBio <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_heatmap_antiTNF_healthy_R_crp.txt", sep = "\t")
ed_noBio <- read.table("Output_files/DESeq2/Heatmap/Remission_crp/DESeq2_normalized_antiTNF_healthy_R_crp.txt", sep = "\t")
sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_R_crp.txt", sep = "\t")
topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")

# Counts matrix
hm_ed <- ed_noBio[rownames(ed_noBio) %in% sig_genes$feature, ] #ENSG00000225982 not in ed_noBio
hm_ed_scaled <- as.data.frame(t(base::scale(t(hm_ed))))
hm_ed_scaled$gene <- unique_ensg2gene[rownames(hm_ed_scaled), ]$hgnc_symbol 
hm_ed_scaled <- hm_ed_scaled[!(is.na(hm_ed_scaled$gene) | hm_ed_scaled$gene == ""), ] #exclude genes without name
rownames(hm_ed_scaled) <- hm_ed_scaled$gene
hm_ed_scaled <- hm_ed_scaled %>% #excluding gene names column
  dplyr::select(!("gene"))

m.hm_ed_scaled <- data.matrix(hm_ed_scaled)


# Setting up left annotation
intersect_sig_topvar <- base::intersect(sig_genes$feature, rownames(topVar))
sig_topVar <- sig_genes[sig_genes$feature %in% intersect_sig_topvar,]
unique_sig <- setdiff(sig_genes$feature, sig_topVar$feature)
only_sig <- sig_genes[sig_genes$feature %in% unique_sig,]


sig_topVar <- sig_topVar[!(is.na(sig_topVar$gene) | sig_topVar$gene == ""), ]
intersect_sig_topvar <- data.frame(sig_topVar$genes)
intersect_sig_topvar$Comparison <- "Sig genes + TVGs"
colnames(intersect_sig_topvar) <- c("Gene", "Comparison") 

only_sig <- only_sig[!(is.na(only_sig$gene) | only_sig$gene == ""), ]
unique_sig <- data.frame(only_sig$genes)
unique_sig$Comparison <- "Unique sig gene"
colnames(unique_sig) <- c("Gene", "Comparison") 

hm_row_annotation <- rbind(unique_sig, intersect_sig_topvar)


# Coldata
coldata_noBio <- coldata_noBio %>%
  mutate(Syst_therapy = recode(No_syst,
                               "0" = "Yes",
                               "1" = "No")) %>% 
  mutate(Prednisolon = recode(Prednisolon,
                              "0" = "No",
                              "1" = "Yes"))%>% 
  mutate(across(c(diagnosis_class, sex, remission, Syst_therapy, antiTNF_vs_noBiologics), as.factor))

# Heatmap
colors_hm <- colorRamp2(c(-4, 0, 4), c("#333b93", "white", "#b61728"))

col_ha_noBio = HeatmapAnnotation(df = data.frame(
  Sex = coldata_noBio$sex,
  BMI = coldata_noBio$bmi,
  Age = coldata_noBio$age,
  Systemic_therapy = coldata_noBio$Syst_therapy,
  Diagnosis = coldata_noBio$diagnosis_class),
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
                           "Healthy" = "green3")),
  CRP = anno_barplot(coldata_noBio$crp))



row_ha = rowAnnotation(df = data.frame(Comparison = hm_row_annotation$Comparison),
                            col = list(Comparison=c("Sig genes + TVGs"="#E69F00", 
                                                    "Unique sig gene" = "#56B4E9")),
                            border = FALSE, gap = unit(2, "mm"))

coldata_noBio$antiTNF_vs_noBiologics <- ordered(coldata_noBio$antiTNF_vs_noBiologics, levels = c("Healthy", "no_biologics","anti_tnf"))
hm = ComplexHeatmap::Heatmap(m.hm_ed_scaled, top_annotation = col_ha_noBio, show_column_names = FALSE, name = "Z score",
                             show_row_names = TRUE, border = TRUE, column_split = coldata_noBio$antiTNF_vs_noBiologics, column_gap = unit(2, "mm"),
                             cluster_column_slices = F, row_title_rot = 90, 
                             left_annotation = row_ha, row_split = hm_row_annotation$Comparison, cluster_row_slices = T,
                             show_column_dend = FALSE, show_row_dend = TRUE, #row_names_gp = gpar(fontsize=9), 
                             col = colors_hm, heatmap_legend_param = list(legend_direction = "horizontal"))

hm = draw(hm, heatmap_legend_side = "bottom")