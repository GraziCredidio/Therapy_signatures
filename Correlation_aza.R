# EZE cohort: Therapy Signatures
  # Correlation of clinical parameters and DEGs: azathioprine x no systemic therapies
  # DEGs, TVGs, variance partition > 25% 
  # Figure 15B (Master's thesis)
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)

# Loading data ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_aza <- read.table("Cleaned_tables/models/aza/EZECohort_coldata_aza.txt", sep = "\t")
vst_counts <- read.table("Cleaned_tables/models/aza/vst_counts_aza.txt", sep = "\t")
sig_genes <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_aza_vs_noSyst.txt", sep = "\t")

top25percent_varPart <- read.table("Output_files/Variance_partition/aza/top25percent_varPart_aza.txt", sep = "\t")

topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")

# Data cleaning ----
# Filtering for patients using Azathioprin
coldata_aza <- coldata_aza %>% 
  filter(Azathioprin == 1)

# Exclude genes without gene names
sig_genes <- sig_genes %>% 
  filter(!(genes == ""))

# DEGs in topvar + varpart>25%
topVar_genes <- rownames(topVar) 
top25percent_varPart_genes <- rownames(top25percent_varPart)

intersect_topVar_varPart25 <- base::intersect(top25percent_varPart_genes, topVar_genes)

# Correlation dataframe----
genes <- intersect_topVar_varPart25
all_correlations <- data.frame()

for (i in 1:length(genes)) {
  deg <- genes[i]
  
  gene_counts <- vst_counts[rownames(vst_counts) == deg,]
  coldata_gene <- as.data.frame(t(gene_counts))
  coldata_gene$sample_id <- rownames(coldata_gene)
  coldata_gene <- left_join(coldata_gene, coldata_aza, by = "sample_id")
  
  correlation_aza_dose <- cor.test(coldata_gene[,1], coldata_gene$aza_dose,  method = "spearman", exact=FALSE)
  correlation_age <- cor.test(coldata_gene[,1], coldata_gene$age,  method = "spearman", exact=FALSE)
  correlation_bmi <- cor.test(coldata_gene[,1], coldata_gene$bmi,  method = "spearman", exact=FALSE)
  correlation_crp <- cor.test(coldata_gene[,1], coldata_gene$crp,  method = "spearman", exact=FALSE)
  correlation_leuco <- cor.test(coldata_gene[,1], coldata_gene$leucocytes,  method = "spearman", exact=FALSE)
  correlation_erythro <- cor.test(coldata_gene[,1], coldata_gene$erythrocytes,  method = "spearman", exact=FALSE)
  correlation_thrombo <- cor.test(coldata_gene[,1], coldata_gene$thrombocytes,  method = "spearman", exact=FALSE)
  
  all_correlations <- rbind(all_correlations, data.frame(row.names = deg, 
                                                         Aza_dose = correlation_aza_dose$estimate, P_value_aza_dose = correlation_aza_dose$p.value,
                                                         Age = correlation_age$estimate, P_value_age = correlation_age$p.value,
                                                         BMI = correlation_bmi$estimate, P_value_bmi = correlation_bmi$p.value,
                                                         CRP = correlation_crp$estimate, P_value_crp = correlation_crp$p.value,
                                                         Leukocytes = correlation_leuco$estimate, P_value_leuco = correlation_leuco$p.value,
                                                         Erythrocytes = correlation_erythro$estimate, P_value_erythro = correlation_erythro$p.value,
                                                         Thrombocytes = correlation_thrombo$estimate, P_value_thrombo = correlation_thrombo$p.value))
  
}

# Rho and p-value data frames
all_correlations_rho <- all_correlations %>% 
  dplyr::select(Aza_dose, Age, BMI, CRP, Leukocytes, Erythrocytes, Thrombocytes) 

all_correlations_p <- all_correlations %>% 
  dplyr::select(P_value_aza_dose,  P_value_age, P_value_bmi, P_value_crp, P_value_leuco, P_value_erythro, P_value_thrombo)

# Calculating adjusted p-values
all_correlations_padj <- all_correlations_p %>% 
  as.matrix %>% 
  as.vector %>% 
  p.adjust(method='BH') %>% 
  matrix(ncol=7)
all_correlations_padj <- as.data.frame(all_correlations_padj)
rownames(all_correlations_padj) <- rownames(all_correlations_p)
colnames(all_correlations_padj) <- colnames(all_correlations_p)

# Replacing gene symbols by gene names----
# rho
all_correlations_rho$gene <- unique_ensg2gene[rownames(all_correlations_rho), ]$hgnc_symbol
all_correlations_rho <- all_correlations_rho[!(is.na(all_correlations_rho$gene) | all_correlations_rho$gene == ""), ]
rownames(all_correlations_rho) <- all_correlations_rho$gene

all_correlations_rho <- all_correlations_rho %>% 
  dplyr::select(!("gene"))

#padj
all_correlations_padj$gene <- unique_ensg2gene[rownames(all_correlations_padj), ]$hgnc_symbol
all_correlations_padj <- all_correlations_padj[!(is.na(all_correlations_padj$gene) | all_correlations_padj$gene == ""), ]
rownames(all_correlations_padj) <- all_correlations_padj$gene

all_correlations_padj <- all_correlations_padj %>% #excluding gene names column
  dplyr::select(!("gene"))

# Heatmap preprocessing ----
# Reordering genes to match heatmap (figure 15A) and transforming as matrix
order <- c("BOK", "CLIC3", "MLC1", "SPON2", "PRF1", "NMUR1","KIR2DL3","KIR2DL1","KIR3DL1","SLC1A7","GNLY","CLDND2","NKG7","FGFBP2","GZMB",   
           "PRSS23", "ADGRG1","TBX21","S1PR5","C1orf21","FCRL6" ,   "CCL4"  ,   "XCL2"    , "FASLG"  ,  "CD160"   , "GZMA" ,    "KLRF1"  ,  
           "KLRD1"  ,  "HOPX"  ,   "NCAM1", "SH2D1B",   "RNF165" ,  "DGKK",     "AKR1C3",   "GRIK4",    "TSPEAR",   "CTSE",     "GADD45A",
           "HEPACAM2", "TMCC2",    "RHAG",     "CR1L",     "C2orf88",  "RPL3L" ,   "SHISA4",  "GDF15" ,   "ITLN1" ,   "KEL",      "AQP1",
           "KCNH2"  ,  "PGF"   ,   "YPEL4"  ,  "DYRK3"  ,  "KRT79"  ,  "ETV7"  ,   "KIR2DL4" , "TCL1A"  ,  "CD72"   ,  "CORO2B")

all_correlations_rho <- all_correlations_rho[order(match(rownames(all_correlations_rho), order)), , drop = FALSE]
all_correlations_rho <- as.matrix(all_correlations_rho)

all_correlations_padj <- all_correlations_padj[order(match(rownames(all_correlations_padj), order)), , drop = FALSE]
all_correlations_padj <- as.matrix(all_correlations_padj)

# Heatmap ----
palette <- colorRampPalette(c("#440d57", "#20928c", "#efe51c"))
palette_hm <- palette(10)

pdf(file = "Output_files/Heatmaps/aza/correlation_heatmap_aza.pdf", width = 10, height= 15)
heatmap_corr <- (Heatmap(all_correlations_rho,
                         cell_fun = function(j, i, x, y, w, h, fill) {
                           if(all_correlations_padj[i, j] < 0.001) {
                             grid.text("***", x, y)
                             } else if(all_correlations_padj[i, j] < 0.01) {
                               grid.text("**", x, y)
                               } else if(all_correlations_padj[i, j] < 0.05) {
                                 grid.text("*", x, y)
                                 }
                           },
                         show_row_dend = FALSE, show_column_dend = FALSE, name = "Spearman's Rho", height = unit(30, "cm"), 
                         width = unit(7, "cm"),col = palette_hm, row_names_side = "left",  column_names_side = "top", 
                         column_names_rot = 45,row_order = order, heatmap_legend_param = list(title = "Spearman's Rho", 
                                                                                              at = c(-1, -0.8, 0, 0.8, 1),
                                                                                              title_position = "topcenter", 
                                                                                              legend_width = unit(7, "cm"),
                                                                                              legend_direction = "horizontal"
                                                                                              )
                         )
                 )

draw(heatmap_corr, heatmap_legend_side = "bottom")
dev.off()
