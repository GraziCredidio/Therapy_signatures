# EZE cohort: Therapy Signatures
  # Correlation of clinical parameters and DEGs: prednisolone x no systemic therapies
  # Top 50 DEGs with absLFC > 0.5 shared with TVGs 
  # Figure 10B (Master's thesis)
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

coldata_pred <- read.table("Cleaned_tables/models/pred/EZECohort_coldata_pred.txt", sep = "\t")
vst_counts <- read.table("Cleaned_tables/models/pred/vst_counts_pred.txt", sep = "\t")
sig_genes <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_pred_vs_noSyst.txt", sep = "\t")

# TVG
topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVar_genes <- rownames(topVar)

# Extracting top50 DEGs (absLFC > 0.5) shared with TVG----
sig_genes_lfc0.5 <- sig_genes %>% 
  filter(abs(coef) > 0.5)

sig_genes_lfc0.5_TVG <- sig_genes_lfc0.5[sig_genes_lfc0.5$feature %in% topVar_genes,]
sig_genes_lfc0.5_TVG <- sig_genes_lfc0.5_TVG %>% 
  arrange(qval)
sig_genes_lfc0.5_TVG <- sig_genes_lfc0.5_TVG[!(is.na(sig_genes_lfc0.5_TVG$genes) | sig_genes_lfc0.5_TVG$genes == ""), ]

top50_sig_genes_lfc0.5 <- sig_genes_lfc0.5_TVG[1:50,]

# Filtering vst_counts to top50 DEGs
top50_counts <- vst_counts[top50_sig_genes_lfc0.5$feature, ]

# Matching gene counts samples to prednisolone patients ----
coldata_pred <- coldata_pred %>% 
  filter(Prednisolon == 1)

top50_counts <- top50_counts[,colnames(top50_counts) %in% coldata_pred$sample_id] 

# Correlation dataframe ----
#genes <- sig_genes[sig_genes$feature %in% sig_genes_topVar_lfc0.5$feature,1]
genes <- row.names(top50_counts)
all_correlations <- data.frame()

for (i in 1:length(genes)) {
  sig_gene <- genes[i]
  
  gene_counts <- top50_counts[rownames(top50_counts) == sig_gene,]
    #vst_counts_R_pred[rownames(vst_counts_R_pred) == sig_gene,]
  coldata_gene <- as.data.frame(t(gene_counts))
  coldata_gene$sample_id <- rownames(coldata_gene)
  coldata_gene <- left_join(coldata_gene, coldata_pred, by = "sample_id")
  
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

# Rho and p-value data frames
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

# Replacing gene symbols by gene names----
rownames(all_correlations_rho) <- unique_ensg2gene[rownames(all_correlations_rho), ]$hgnc_symbol
rownames(all_correlations_padj) <- unique_ensg2gene[rownames(all_correlations_padj), ]$hgnc_symbol

# Heatmap preprocessing----
# Reordering genes to match heatmap (figure 10A) and transforming as matrix
order <- c("DBH-AS1", "LIMS2", "ID3", "SPIB", "IRS2", "MIR223HG", "CLEC4E", "TLR2", "FKBP5",
           "KCNE1", "ECHDC3", "TPST1", "INHBB", "NSUN7", "DUSP1", "TSC22D3", "GADD45A", "ASPH",     
           "IRAK3", "ZNF608", "VNN1", "GPR141", "SIPA1L2", "IL18R1", "PFKFB2", "SLC8A1", "COL9A2",  
           "VSIG4", "DAAM2", "MAOA", "FLT3", "AMPH", "IL1R2", "ADAMTS2", "OLAH", "GRB10", "PCSK9", 
           "ARG1", "LINC02207", "ABCG1", "PPIAP29", "SLC5A9", "DUSP13", "IFITM3P2", "CD163", "GPER1", 
           "C5orf67", "SLC22A4", "TLR5", "NAIP")

all_correlations_rho <- all_correlations_rho[order(match(rownames(all_correlations_rho), order)), , drop = FALSE]
all_correlations_rho <- as.matrix(all_correlations_rho)

all_correlations_padj <- all_correlations_padj[order(match(rownames(all_correlations_padj), order)), , drop = FALSE]
all_correlations_padj <- as.matrix(all_correlations_padj)

# Heatmap ----
palette <- colorRampPalette(c("#440d57", "#20928c", "#efe51c"))
palette_hm <- palette(10)

pdf(file = "Output_files/Heatmaps/pred/correlation_heatmap_pred.pdf", width = 10, height= 15)
heatmap_corr <- Heatmap(all_correlations_rho,
                        cell_fun = function(j, i, x, y, w, h, fill) {
                          if(all_correlations_padj[i, j] < 0.001) {
                            grid.text("***", x, y)
                            } else if(all_correlations_padj[i, j] < 0.01) {
                              grid.text("**", x, y)
                              } else if(all_correlations_padj[i, j] < 0.05) {
                                grid.text("*", x, y)
                                }
                          },
                        show_row_dend = FALSE, show_column_dend = FALSE, name = "Spearman's Rho", 
                        height = unit(30, "cm"), width = unit(7, "cm"), col = palette_hm, 
                        row_names_side = "left",  column_names_side = "top", column_names_rot = 45,
                        row_order = order,heatmap_legend_param = list(title = "Spearman's Rho", 
                                                                      at = c(-1, -0.8, 0, 0.8, 1),
                                                                      title_position = "topcenter", 
                                                                      legend_width = unit(7, "cm"),
                                                                      legend_direction = "horizontal")
                        )

draw(heatmap_corr, heatmap_legend_side = "bottom")
dev.off()
