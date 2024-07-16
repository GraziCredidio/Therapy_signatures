# EZE cohort: Therapy Signatures
  # Heatmap: anti-TNF x no biologics
  # DEGs shared with TVGs and not shared with TVGs (unique)
  # Figure 5A
  # Author: Graziella Credidio

rm(list = ls())

folder <- "Output_files/Heatmaps/antiTNF"
if (!dir.exists(folder)) {
  dir.create(folder)
}

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
topVar_genes <- rownames(topVar)

# heatmap coldata and counts
coldata_antiTNF <- read.table("Cleaned_tables/heatmaps/hm_coldata_antiTNF.txt", sep = "\t")
normalized_counts <- read.table("Cleaned_tables/heatmaps/hm_normalized_counts_antiTNF.txt", sep = "\t") 

# DEGs
sig_genes <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_antiTNF_vs_noBio.txt", sep = "\t")

# Coldata preprocessing ----
coldata_antiTNF <- coldata_antiTNF %>%
  mutate(Syst_therapy = recode(No_syst,
                               "0" = "Yes",
                               "1" = "No")) %>% 
  mutate(Prednisolon = recode(Prednisolon,
                              "0" = "No",
                              "1" = "Yes"))%>% 
  mutate(across(c(diagnosis_class, sex, remission, Syst_therapy, antiTNF_vs_noBiologics), as.factor))

coldata_antiTNF$antiTNF_vs_noBiologics <- ordered(coldata_antiTNF$antiTNF_vs_noBiologics, levels = c("Healthy", "no_biologics","anti_tnf"))

# Counts preprocessing ----
hm_ed <- normalized_counts[rownames(normalized_counts) %in% sig_genes$feature, ] #ENSG00000225982 not in ed_noBio
hm_ed_scaled <- as.data.frame(t(base::scale(t(hm_ed))))
hm_ed_scaled$gene <- unique_ensg2gene[rownames(hm_ed_scaled), ]$hgnc_symbol 
hm_ed_scaled <- hm_ed_scaled[!(is.na(hm_ed_scaled$gene) | hm_ed_scaled$gene == ""), ] #exclude genes without name
rownames(hm_ed_scaled) <- hm_ed_scaled$gene
hm_ed_scaled <- hm_ed_scaled %>% 
  dplyr::select(!("gene"))

m.hm_ed_scaled <- data.matrix(hm_ed_scaled)

# Heatmap row annotation ----
# DEGs shared with TVGs
intersect_sig_topvar <- base::intersect(sig_genes$feature, rownames(topVar))
sig_topVar <- sig_genes[sig_genes$feature %in% intersect_sig_topvar,]
sig_topVar <- sig_topVar[!(is.na(sig_topVar$gene) | sig_topVar$gene == ""), ]
intersect_sig_topvar <- data.frame(sig_topVar$genes)
intersect_sig_topvar$Comparison <- "Shared (DEGs + TVGs)"
colnames(intersect_sig_topvar) <- c("Gene", "Comparison") 

# DEGs not shared with TVGs
unique_sig <- setdiff(sig_genes$feature, sig_topVar$feature)
only_sig <- sig_genes[sig_genes$feature %in% unique_sig,]
only_sig <- only_sig[!(is.na(only_sig$gene) | only_sig$gene == ""), ]
unique_sig <- data.frame(only_sig$genes)
unique_sig$Comparison <- "Unique DEGs"
colnames(unique_sig) <- c("Gene", "Comparison") 

hm_row_annotation <- rbind(unique_sig, intersect_sig_topvar)

# Heatmap plotting ----
colors_hm <- colorRamp2(c(-4, 0, 4), c("#333b93", "white", "#b61728"))

col_ha_antiTNF = HeatmapAnnotation(df = data.frame(
  Sex = coldata_antiTNF$sex,
  BMI = coldata_antiTNF$bmi,
  Age = coldata_antiTNF$age,
  Systemic_therapy = coldata_antiTNF$Syst_therapy,
  Diagnosis = coldata_antiTNF$diagnosis_class),
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
  CRP = anno_barplot(coldata_antiTNF$crp))


pdf(file = "Output_files/Heatmaps/antiTNF/heatmap_antiTNF_healthy.pdf", width = 10, height= 12)
row_ha = rowAnnotation(df = data.frame(Comparison = hm_row_annotation$Comparison),
                            col = list(Comparison=c("Shared (DEGs + TVGs)"="#E69F00", 
                                                    "Unique DEGs" = "#56B4E9")),
                            border = FALSE, gap = unit(2, "mm"))

hm = ComplexHeatmap::Heatmap(m.hm_ed_scaled, top_annotation = col_ha_antiTNF, show_column_names = FALSE, name = "Z score",
                             show_row_names = TRUE, border = TRUE, column_split = coldata_antiTNF$antiTNF_vs_noBiologics, column_gap = unit(2, "mm"),
                             cluster_column_slices = F, row_title_rot = 90, 
                             left_annotation = row_ha, row_split = hm_row_annotation$Comparison, cluster_row_slices = T,
                             show_column_dend = FALSE, show_row_dend = TRUE,  
                             col = colors_hm, heatmap_legend_param = list(legend_direction = "horizontal"))

hm = draw(hm, heatmap_legend_side = "bottom")
hm
dev.off()
