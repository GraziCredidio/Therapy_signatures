# EZE cohort: Therapy Signatures
  # UpSet plot: azathioprine DEGs, DNAm-DEGs, varPart > 25% and TVG gene sets
  # Figure 18
  # Author: Graziella Credidio


rm(list = ls())

# Loading packages ----
library(UpSetR)
library(tidyverse)

folder <- "Output_files/UpSet"
if (!dir.exists(folder)) {
    dir.create(folder)}


# Loading data ----
topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVar_genes <- rownames(topVar)

# DNAm DEGs
aza_dnam_up <- read.table("Output_files/Methylation/aza/aza_upregulated_DNAm_DEG.txt", sep = '\t')
aza_dnam_up <- unique(aza_dnam_up$Gene)

aza_dnam_down <- read.table("Output_files/Methylation/aza/aza_downregulated_DNAm_DEG.txt", sep = '\t')
aza_dnam_down <- unique(aza_dnam_down$Gene)

# DEGs
aza_degs <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_aza_vs_noSyst.txt", sep = '\t')
aza_up_degs <- aza_degs %>% 
  filter(coef > 0) %>% 
  pull(feature)

aza_down_degs <- aza_degs %>% 
  filter(coef < 0) %>% 
  pull(feature)

# VarPart > 25%
aza_top25percent_varPart <- read.table("Output_files/Variance_partition/aza/top25percent_varPart_aza.txt")
aza_up_top25percent_varPart <- intersect(aza_up_degs, rownames(aza_top25percent_varPart))
aza_down_top25percent_varPart <- intersect(aza_down_degs, rownames(aza_top25percent_varPart))


# Plotting ----
# Upregualted
upregulated_genes <- list(
  `DNAm-DEGs` = aza_dnam_up,
  `DEGs` = aza_up_degs,
  `VarPart>25%` = aza_up_top25percent_varPart,
  `TVGs` = topVar_genes
)

pdf("Output_files/UpSet/upset_aza_up.pdf", width = 10, height = 7)
up <- upset(fromList(upregulated_genes), order.by = "freq", nintersects = NA, 
              point.size = 3,  text.scale = 2, sets.x.label = "# of genes",
              shade.color = "#FF9999", matrix.color = "#990000",
              main.bar.color = "#990000", mainbar.y.label = " ",
              sets.bar.color = "#990000", att.color = "#990000",
              set_size.show =TRUE)
up
dev.off()

# Downregulated
downregulated_genes <- list(
  `DNAm-DEGs` = aza_dnam_down,
  `DEGs` = aza_down_degs,
  `VarPart>25%` = aza_down_top25percent_varPart,
  `TVGs` = topVar_genes
)

pdf("Output_files/UpSet/upset_aza_down.pdf", width = 10, height = 7)
down <- upset(fromList(downregulated_genes), order.by = "freq", nintersects = NA, 
      point.size = 3, text.scale = 2, sets.x.label = "# of genes",
      shade.color = "lightblue", matrix.color = "skyblue4",
      main.bar.color = "skyblue4", mainbar.y.label = " ",
      sets.bar.color = "skyblue4", att.color = "skyblue4",
      set_size.show =TRUE)

down
dev.off()
