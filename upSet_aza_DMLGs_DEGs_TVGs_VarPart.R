# EZE cohort - Aza model (R + crp)
# Upset plot (1 up + 1 down): all DEGs + DMLGs + Varpart 25% + overall TVGs + inactive TVGs

library(UpSetR)
library(tidyverse)


# overall TVGs:
topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVarGenes_names <- rownames(topVarGenes) 

# DMLGs:
aza_up_correlation <- read.table("Output_files/Methylation/aza/aza_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')  #498
aza_up_correlation <- aza_up_correlation$Gene
aza_up_correlation <- unique(aza_up_correlation)

aza_down_correlation <- read.table("Output_files/Methylation/aza/aza_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t') #1055
aza_down_correlation <- aza_down_correlation$Gene
aza_down_correlation <- unique(aza_down_correlation)

# DEGs:
R_crp_aza <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_crp.txt", sep = "\t")
pos_sig_aza <- R_crp_aza %>%
  filter(coef > 0)%>% 
  pull(feature)

neg_sig_aza <- R_crp_aza %>%
  filter(coef < 0)%>% 
  pull(feature)

# VarPart > 25%:
aza_varPart_25 <- read.table("Output_files/VariancePartition/Remission_crp/aza/top25percent_varPart_aza_R_crp.txt", sep = "\t")
aza_varPart_25_up <- intersect(pos_sig_aza, rownames(aza_varPart_25))

aza_varPart_25_down <- intersect(neg_sig_aza, rownames(aza_varPart_25))

# Selecting just the "features" column and plotting ----

upregulated_genes <- list(
  `DMLGs` = aza_up_correlation,
  `DEGs` = pos_sig_aza,
  `VarPart>25%` = aza_varPart_25_up,
  `TVGs` = topVarGenes_names
)

pdf("Output_files/Methylation/aza/Upset/Upset_aza_up_DMLGs_DEGs_VarPart_TVGs.pdf", width = 8, height = 7)
up <- upset(fromList(upregulated_genes), order.by = "freq", nintersects = NA, 
              point.size = 3,  text.scale = 2, sets.x.label = "# of genes",
              shade.color = "#FF9999", matrix.color = "#990000",
              main.bar.color = "#990000", mainbar.y.label = " ",
              sets.bar.color = "#990000", att.color = "#990000",
              set_size.show =TRUE)
up
dev.off()



downregulated_genes <- list(
  `DMLGs` = aza_down_correlation,
  `DEGs` = neg_sig_aza,
  `VarPart>25%` = aza_varPart_25_down,
  `TVGs` = topVarGenes_names
)

pdf("Output_files/Methylation/aza/Upset/Upset_aza_down_DMLGs_DEGs_VarPart_TVGs.pdf", width = 8, height = 7)
down <- upset(fromList(downregulated_genes), order.by = "freq", nintersects = NA, 
      point.size = 3, text.scale = 2, sets.x.label = "# of genes",
      shade.color = "lightblue", matrix.color = "skyblue4",
      main.bar.color = "skyblue4", mainbar.y.label = " ",
      sets.bar.color = "skyblue4", att.color = "skyblue4",
      set_size.show =TRUE)

down
dev.off()

# Extracting genes from each intersection
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

down_all_intersections <- as.data.frame(intersect(neg_sig_aza, intersect(aza_down_correlation, intersect(topVarGenes_names, aza_varPart_25_down))))
down_all_intersections$genes <- unique_ensg2gene[down_all_intersections$`intersect(neg_sig_aza, intersect(aza_down_correlation, intersect(topVarGenes_names, aza_varPart_25_down)))`,]$hgnc_symbol

up_all_intersections <- as.data.frame(intersect(pos_sig_aza, intersect(aza_up_correlation, intersect(topVarGenes_names, aza_varPart_25_up))))
up_all_intersections$genes <- unique_ensg2gene[up_all_intersections$`intersect(pos_sig_aza, intersect(aza_up_correlation, intersect(topVarGenes_names, aza_varPart_25_up)))`,]$hgnc_symbol



# TVGs + Varpart but not DMLGs
TVG_varpart_up <- setdiff(intersect(topVarGenes_names, aza_varPart_25_up), aza_up_correlation)
TVG_varpart_down <- setdiff(intersect(topVarGenes_names, aza_varPart_25_down), aza_down_correlation)
