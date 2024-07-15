# EZE Cohort - therapy signatures
# Methylation data integration with transcriptome
# Venn diagrams

rm(list = ls())
setwd("C:/Documents/Masters thesis/EZE_cohort") # laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC


aza_correlation <- read.table("Output_files/Methylation/aza/aza_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = "\t") #1553
pred_correlation <- read.table("Output_files/Methylation/pred/pred_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = "\t") #27045
antiTNF_correlation <- read.table("Output_files/Methylation/antiTNF/antiTNF_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = "\t") #118

# Canonical methylation (more meth, less gene expression and vice-versa)
aza_canonical <- aza_correlation %>% 
  filter(Rho < 0) #945 sites
aza_canonical <- unique(aza_canonical$Gene_name) #302 genes (277 with non empty gene names)


pred_canonical <- pred_correlation %>% 
  filter(Rho < 0) #15656 sites
pred_canonical <- unique(pred_canonical$Gene_name) #4005 genes (3729 with non empty gene names)


antiTNF_canonical <- antiTNF_correlation %>% 
  filter(Rho < 0) #72
antiTNF_canonical <- unique(antiTNF_canonical$Gene_name) #34 genes (33 with non empty gene names)

# Venn diagram of up and downregulated genes ----
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
install.packages("rlist")

library(tidyverse)
library(ggvenn)
library(circlize)
library(rlist)

topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVarGenes_names <- rownames(topVarGenes) 

# Methylation
aza_up_correlation <- read.table("Output_files/Methylation/aza/aza_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')  #498
aza_up_correlation <- aza_up_correlation$Gene
aza_up_correlation <- unique(aza_up_correlation)
aza_down_correlation <- read.table("Output_files/Methylation/aza/aza_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t') #1055
aza_down_correlation <- aza_down_correlation$Gene
aza_down_correlation <- unique(aza_down_correlation)

pred_up_correlation <- read.table("Output_files/Methylation/pred/pred_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t') #14675
pred_up_correlation <- pred_up_correlation$Gene
pred_up_correlation <- unique(pred_up_correlation)
pred_down_correlation <- read.table("Output_files/Methylation/pred/pred_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t') #12370
pred_down_correlation <- pred_down_correlation$Gene
pred_down_correlation <- unique(pred_down_correlation)

noBio_up_correlation <- read.table("Output_files/Methylation/antiTNF/antiTNF_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')
noBio_up_correlation <- noBio_up_correlation$Gene
noBio_up_correlation <- unique(noBio_up_correlation)
noBio_down_correlation <- read.table("Output_files/Methylation/antiTNF/antiTNF_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')
noBio_down_correlation <- noBio_down_correlation$Gene
noBio_down_correlation <- unique(noBio_down_correlation)

# Sig genes (inactive + CRP)
setwd("C:/Documents/Masters thesis/EZE_cohort/Output_files/Maaslin2/Results/Remission/significant_results") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort/Output_files/Maaslin2/Results/Remission_crp/significant_results") #PC
files_R_crp = list.files(pattern="*.txt")
myfiles_R_crp = lapply(files_R_crp, read.delim)
R_crp_noBio <- myfiles_R_crp[[1]]
R_crp_aza <- myfiles_R_crp[[2]]
R_crp_pred <- myfiles_R_crp[[3]]

pos_sig_aza <- R_crp_aza %>%
  filter(coef > 0)%>% 
  pull(feature)
neg_sig_aza <- R_crp_aza %>%
  filter(coef < 0)%>% 
  pull(feature)

pos_sig_pred <-  R_crp_pred %>% 
  filter(coef > 0)%>% 
  pull(feature)
neg_sig_pred <-   R_crp_pred %>% 
  filter(coef < 0)%>% 
  pull(feature)

pos_sig_noBio <- R_crp_noBio %>% 
  filter(coef > 0) %>% 
  pull(feature)
neg_sig_noBio <- R_crp_noBio %>% 
  filter(coef < 0)%>% 
  pull(feature)

# aza: varpart > 25% + DMLGs + TVGs
aza_varPart_25 <- read.table("Output_files/VariancePartition/Remission_crp/aza/top25percent_varPart_aza_R_crp.txt", sep = "\t")
aza_varPart_25_up <- intersect(pos_sig_aza, rownames(aza_varPart_25))
aza_varPart_25_down <- intersect(neg_sig_aza, rownames(aza_varPart_25))

# Lists
aza_down_corr_sig <- list(`Aza x no Syst (sig)` = neg_sig_aza,
                          `Aza x no Syst (corr meth)` = aza_down_correlation,
                          `Top var genes` = topVarGenes_names)
aza_pos_corr_sig <- list(`Aza x no Syst (sig)` = pos_sig_aza,
                         `Aza x no Syst (corr meth)` = aza_up_correlation,
                         `Top var genes` = topVarGenes_names)

aza_pos_corr_varPart <- list(`VarPart > 25%` = aza_varPart_25_up,
                             `Aza x no Syst (corr meth)` = aza_up_correlation,
                             `TVGs` = topVarGenes_names)
aza_neg_corr_varPart <- list(`VarPart > 25%` = aza_varPart_25_down,
                             `Aza x no Syst (corr meth)` = aza_down_correlation,
                             `TVGs` = topVarGenes_names)



pred_down_corr_sig <- list(`Pred x no Syst (sig)` = neg_sig_pred,
                           `Pred x no Syst (corr meth)` = pred_down_correlation,
                           `Top var genes` = topVarGenes_names)
pred_pos_corr_sig <- list(`Pred x no Syst (sig)` = pos_sig_pred,
                          `Pred x no Syst (corr meth)` = pred_up_correlation,
                          `Top var genes` = topVarGenes_names)


noBio_down_corr_sig <- list(`Anti TNF x no Bio (sig)` = neg_sig_noBio,
                            `Anti TNF x no Bio (corr meth)` = noBio_down_correlation,
                            `Top var genes` = topVarGenes_names)
noBio_pos_corr_sig <- list(`Anti TNF x no Bio (sig)` = pos_sig_noBio,
                           `Anti TNF x no Bio (corr meth)` = noBio_up_correlation,
                           `Top var genes` = topVarGenes_names)


# Plots
setwd("C:/Documents/Masters thesis/EZE_cohort") # laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

pos_fillColors <- c("#dfa693", "#dc6e55", "#bacb87")
neg_fillColors <- c("#8bbdd9", "#0079bf", "#bacb87")

Venn_plots <- function(data, fillColors, fileName) {
  Venn <- ggvenn(data = data, stroke_linetype = 2, stroke_size = 0.5, set_name_size = 8, text_size = 9,  show_percentage = FALSE, 
                 fill_color = fillColors)
  
  ggsave(filename = fileName, plot = Venn, 
         height = 10, width = 10, units = "in", dpi = 300, bg = "white")
  
  Venn
  
}


# Methyl.correlation + sig genes + top 2000 var
Venn_plots(aza_pos_corr_sig, pos_fillColors,"Output_files/Methylation/aza/Venn/venn_aza_corr_topVar_pos.pdf")
Venn_plots(aza_down_corr_sig, neg_fillColors,"Output_files/Methylation/aza/Venn/venn_aza_corr_topVar_neg.pdf")

Venn_plots(aza_pos_corr_varPart, pos_fillColors,"Output_files/Methylation/aza/Venn/varPart_TVGs/venn_aza_DMLGs_varPart_TVGs_pos.pdf")
Venn_plots(aza_neg_corr_varPart, neg_fillColors,"Output_files/Methylation/aza/Venn/varPart_TVGs/venn_aza_DMLGs_varPart_TVGs_neg.pdf")


Venn_plots(pred_pos_corr_sig, pos_fillColors,"Output_files/Methylation/pred/Venn/venn_pred_corr_topVar_pos.png")
Venn_plots(pred_down_corr_sig, neg_fillColors,"Output_files/Methylation/pred/Venn/venn_pred_corr_topVar_neg.png")

Venn_plots(noBio_pos_corr_sig, pos_fillColors,"Output_files/Methylation/antiTNF/Venn/venn_antiTNF_corr_topVar_pos.png")
Venn_plots(noBio_down_corr_sig, neg_fillColors,"Output_files/Methylation/antiTNF/Venn/venn_antiTNF_corr_topVar_neg.png")


# intersection genes: methylation + sig + topvar ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install("biomaRt")
library(biomaRt)
install.packages("kableExtra")
library(kableExtra)

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")

# Aza
aza_up_intersection <- intersect(aza_up_correlation, intersect(topVarGenes_names, pos_sig_aza))
aza_up_intersection_description <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name"), filters = "ensembl_gene_id", 
                                         values = aza_up_intersection, mart = mart)
write.table(aza_up_intersection_description, "Output_files/Methylation/aza/Intersection/intersect_genes_methylation_topvar_up_aza_R_crp_description.txt", sep = "\t")

aza_up_intersection_description %>%
  kbl() %>%
  kable_styling()

aza_down_intersection <- intersect(aza_down_correlation, intersect(topVarGenes_names, neg_sig_aza))
aza_down_intersection_description <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name"), filters = "ensembl_gene_id", 
                                           values = aza_down_intersection, mart = mart)
write.table(aza_down_intersection_description, "Output_files/Methylation/aza/Intersection/intersect_genes_methylation_topvar_down_aza_R_crp_description.txt", sep = "\t")

aza_down_intersection_description %>%
  kbl() %>%
  kable_styling()

# Pred
pred_up_intersection <- intersect(pred_up_correlation, intersect(topVarGenes_names, pos_sig_pred))
pred_up_intersection_description <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name"), filters = "ensembl_gene_id", 
                                          values = pred_up_intersection, mart = mart)
write.table(pred_up_intersection_description, "Output_files/Methylation/pred/Intersection/intersect_genes_methylation_topvar_up_pred_R_crp_description.txt", sep = "\t")
pred_up_intersection_description %>%
  kbl() %>%
  kable_styling()

pred_down_intersection <- intersect(pred_down_correlation, intersect(topVarGenes_names, neg_sig_pred))
pred_down_intersection_description <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name"), filters = "ensembl_gene_id", 
                                            values = pred_down_intersection, mart = mart)
write.table(pred_down_intersection_description, "Output_files/Methylation/pred/Intersection/intersect_genes_methylation_topvar_down_pred_R_crp_description.txt", sep = "\t")
pred_down_intersection_description %>%
  kbl() %>%
  kable_styling()


# anti-TNF
antiTNF_up_intersection <- intersect(noBio_up_correlation, intersect(topVarGenes_names, pos_sig_noBio))
antiTNF_up_intersection_description <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name"), filters = "ensembl_gene_id", 
                                             values = antiTNF_up_intersection, mart = mart)
write.table(antiTNF_up_intersection_description, "Output_files/Methylation/antiTNF/Intersection/intersect_genes_methylation_topvar_up_antiTNF_R_crp_description.txt", sep = "\t")
antiTNF_up_intersection_description %>%
  kbl() %>%
  kable_styling()

antiTNF_down_intersection <- intersect(noBio_down_correlation, intersect(topVarGenes_names, neg_sig_noBio))
antiTNF_down_intersection_description <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name"), filters = "ensembl_gene_id", 
                                               values = antiTNF_down_intersection, mart = mart)
write.table(antiTNF_down_intersection_description, "Output_files/Methylation/antiTNF/Intersection/intersect_genes_methylation_topvar_down_antiTNF_R_crp_description.txt", sep = "\t")
antiTNF_down_intersection_description %>%
  kbl() %>%
  kable_styling()







# intersection genes: correlated genes >25% var part + top var ----
rm(list = ls())

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")

# Aza
aza_top25percent <- read.table("Output_files/Methylation/aza/varPart/varPart_aza_corGenes_25percent.txt", sep = "\t")
aza_top25percent <- rownames(aza_top25percent)
aza_up_intersection <- intersect(aza_top25percent, intersect(topVarGenes_names, aza_up_correlation))

aza_up_intersection_description <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name"), filters = "ensembl_gene_id",
                                         values = aza_up_intersection, mart = mart)
write.table(aza_up_intersection_description, "Output_files/Methylation/aza/Intersection/top25percent/aza_intersect_up_corGenes_topVar_25percent_description.txt", sep = "\t")
aza_up_intersection_description %>%
  kbl() %>%
  kable_styling()

aza_down_intersection <- intersect(aza_top25percent, intersect(topVarGenes_names, aza_down_correlation))
aza_down_intersection_description <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name"), filters = "ensembl_gene_id", 
                                           values = aza_down_intersection, mart = mart)
write.table(aza_down_intersection_description, "Output_files/Methylation/aza/Intersection/top25percent/aza_intersect_down_corGenes_topVar_25percent_description.txt", sep = "\t")
aza_down_intersection_description %>%
  kbl() %>%
  kable_styling()

# pred
pred_top25percent <- read.table("Output_files/Methylation/pred/varPart/varPart_pred_corGenes_25percent.txt", sep = "\t")
pred_top25percent <- rownames(pred_top25percent)
pred_up_intersection <- intersect(pred_top25percent, intersect(topVarGenes_names, pred_up_correlation))

pred_up_intersection_description <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name"), filters = "ensembl_gene_id",
                                          values = pred_up_intersection, mart = mart)
write.table(pred_up_intersection_description, "Output_files/Methylation/pred/Intersection/top25percent/pred_intersect_up_corGenes_topVar_25percent_description.txt", sep = "\t")
pred_up_intersection_description %>%
  kbl() %>%
  kable_styling()

pred_down_intersection <- intersect(pred_top25percent, intersect(topVarGenes_names, pred_down_correlation)) #0




###
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

a <- data.frame(intersect(noBio_up_correlation, topVarGenes_names))
a$gene <- unique_ensg2gene[a$intersect.noBio_up_correlation..topVarGenes_names., ]$hgnc_symbol

b <- data.frame(intersect(noBio_down_correlation, topVarGenes_names))
b$gene <- unique_ensg2gene[b$intersect.noBio_down_correlation..topVarGenes_names., ]$hgnc_symbol
