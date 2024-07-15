## Therapy signature analysis - EZE cohort
# Venn diagram of comparisons

graphics.off()
rm(list = ls())

setwd("C:/Documents/Masters thesis/EZE_cohort") # laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

# Loading packages ----
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
install.packages("rlist")

library(tidyverse)
library(ggvenn)
library(circlize)
library(rlist)
library(DESeq2)

# Top var genes ----
# Loading files
coldata <- read.csv("Cleaned_tables/EZECohort_ord.coldata_maaslin_28.02.txt", sep = "\t")
counts <- read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_28.02.txt", sep = "\t")

# Tranforming as factor
coldata <- coldata%>%
  mutate(across(c(diagnosis_class, age_group2, sex, remission, bmi_class, No_syst,
                  smoking, biologics, biologics_TNF), as.factor))

# DESeq2 object
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = coldata,
                                     design = ~ 1)
dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- as.data.frame(assay(vst))

# Top variating genes
topVarGenes_idx <- head(order(rowVars(assay(vst)),decreasing=TRUE),5000) #5000 or 2000
topVarGenes <- vst_counts_norm[topVarGenes_idx,]
topVarGenes_names <- rownames(topVarGenes) 

write.table(topVarGenes, "Cleaned_tables/topVar_5000_genes.txt", sep = "\t", quote = FALSE)

# Loading files ----
topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVarGenes_names <- rownames(topVarGenes) 

# All patients ----
# Venn: anti tnf: q 0.1; mtx: q 0.1
setwd("C:/Documents/Masters thesis/EZE_cohort/Output_files/Maaslin2/significant_results") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort/Output_files/Maaslin2/significant_results") #pc

files = list.files(pattern="*.txt")
myfiles = lapply(files, read.delim)

noBio_q0.05 <- myfiles[[1]]$feature
nonTNF_q0.05 <- myfiles[[3]]$feature
  
noBio  <- myfiles[[2]]$feature
nonTNF  <- myfiles[[4]]$feature

aza_noAza <- myfiles[[5]]$feature
aza_noSyst <- myfiles[[6]]$feature

mtx_noMtx <- myfiles[[8]]$feature
mtx_noSyst <- myfiles[[10]]$feature

mtx_noMtx_q0.05  <- myfiles[[7]]$feature
mtx_noSyst_q0.05 <- myfiles[[9]]$feature

pred_noPred <- myfiles[[11]]$feature
pred_noSyst <- myfiles[[12]]$feature


# Lists with significant genes
setwd("C:/Documents/Masters thesis/EZE_cohort") # laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

sig_TNF <- list(`Anti-TNF vs Non-anti-TNF` = nonTNF,
             `Anti-TNF vs No biologics` = noBio) #q 0.1

sig_aza <- list(`Aza vs No aza` = aza_noAza,
                `Aza vs No systemic` = aza_noSyst)

sig_mtx <- list(`MTX vs No MTX` = mtx_noMtx,
                `MTX vs No systemic` = mtx_noSyst) #q 0.1

sig_pred <- list(`Pred vs No pred` = pred_noPred,
                 `Pred vs No systemic` = pred_noSyst)

# Lists with significant genes and topvar genes 
sig_TNF_topVar <- list(`Anti-TNF vs Non-anti-TNF` = nonTNF,
                       `Anti-TNF vs No biologics` = noBio,
                       `Top var genes` = topVarGenes_names)


sig_noBio_q0.05 <- list( `Anti-TNF vs Non-anti-TNF` = nonTNF_q0.05,
                         `Anti-TNF vs No biologics` = noBio_q0.05,
                        `Top var genes` = topVarGenes_names)

sig_aza_topVar <- list(`Aza vs No aza` = aza_noAza,
                       `Aza vs No systemic` = aza_noSyst,
                       `Top var genes` = topVarGenes_names)

sig_mtx_topVar <- list(`MTX vs No MTX` = mtx_noMtx,
                       `MTX vs No systemic` = mtx_noSyst,
                       `Top var genes` = topVarGenes_names)

sig_mtx_topVar_q0.05 <- list(`MTX vs No MTX` = mtx_noMtx_q0.05 ,
                       `MTX vs No systemic` = mtx_noSyst_q0.05 ,
                       `Top var genes` = topVarGenes_names)

sig_pred_topVar <- list(`Pred vs No pred` = pred_noPred,
                        `Pred vs No systemic` = pred_noSyst,
                        `Top var genes` = topVarGenes_names)

###################################### Remitters ##############################################
setwd("C:/Documents/Masters thesis/EZE_cohort/Output_files/Maaslin2/Results/Remission/significant_results")

files_R = list.files(pattern="*.txt")
myfiles_R = lapply(files_R, read.delim)

R_noBio <- myfiles_R[[1]]$feature

R_aza_noAza <- myfiles_R[[2]]$feature
R_aza_noSyst <- myfiles_R[[3]]$feature

R_mtx_noMtx <- myfiles_R[[4]]$feature
R_mtx_noSyst <- myfiles_R[[5]]$feature

R_pred_noPred <- myfiles_R[[6]]$feature
R_pred_noSyst <- myfiles_R[[7]]$feature

# Lists with significant genes
setwd("C:/Documents/Masters thesis/EZE_cohort") # laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

R_sig_aza <- list(`Aza vs No aza` = R_aza_noAza,
                `Aza vs No systemic` = R_aza_noSyst)

R_sig_mtx <- list(`MTX vs No MTX` = R_mtx_noMtx,
                `MTX vs No systemic` = R_mtx_noSyst)

R_sig_pred <- list(`Pred vs No pred` = R_pred_noPred,
                 `Pred vs No systemic` = R_pred_noSyst)

# Lists with significant genes and topvar genes 
R_sig_TNF_topVar <- list(`Anti-TNF vs No biologics` = R_noBio,
                       `Top var genes` = topVarGenes_names)

R_sig_aza_topVar <- list(`Aza vs No aza` = R_aza_noAza,
                       `Aza vs No systemic` = R_aza_noSyst,
                       `Top var genes` = topVarGenes_names)

R_sig_mtx_topVar <- list(`MTX vs No MTX` = R_mtx_noMtx,
                       `MTX vs No systemic` = R_mtx_noSyst,
                       `Top var genes` = topVarGenes_names)

R_sig_pred_topVar <- list(`Pred vs No pred` = R_pred_noPred,
                        `Pred vs No systemic` = R_pred_noSyst,
                        `Top var genes` = topVarGenes_names)           

##################################### All patients x  remitters ##################################
# significant genes
sig_R_all_noBio <- list(`Anti TNF x no Biologics (all)` = noBio_q0.05,
                      `Anti TNF x no Biologics (inactive)` = R_noBio)

sig_R_all_aza <- list(`Aza x No systemic (all)` = aza_noSyst,
                      `Aza x No systemic (inactive)` = R_aza_noSyst)

sig_R_all_mtx <- list(`MTX x No systemic (all)` = mtx_noSyst_q0.05,
                      `MTX x No systemic (inactive)` = R_mtx_noSyst)

sig_R_all_pred <- list(`Pred vs No systemic (all)` = pred_noSyst,
                      `Pred vs No systemic (inactive)` = R_pred_noSyst)

# significant genes and topVar
sig_R_all_noBio_topVar <- list(`Anti TNF x no Biologics (all)` = noBio_q0.05,
                        `Anti TNF x no Biologics (inactive)` = R_noBio,
                        `Top var genes` = topVarGenes_names)

sig_R_all_aza_topVar <- list(`Aza x No systemic (all)` = aza_noSyst,
                      `Aza x No systemic (inactive)` = R_aza_noSyst,
                      `Top var genes` = topVarGenes_names)

sig_R_all_mtx_topVar <- list(`MTX x No systemic (all)` = mtx_noSyst_q0.05,
                      `MTX x No systemic (inactive)` = R_mtx_noSyst,
                      `Top var genes` = topVarGenes_names)

sig_R_all_pred_topVar <- list(`Pred vs No systemic (all)` = pred_noSyst,
                       `Pred vs No systemic (inactive)` = R_pred_noSyst,
                       `Top var genes` = topVarGenes_names)
                    
# Venn diagrams ----
fillColors2 <- c("pink", "gold")
fillColors3 <- c("pink", "gold", "lightgreen")

Venn_sig <- function(data, fillColors, fileName) {
  Venn <- ggvenn(data = data, stroke_linetype = 2, stroke_size = 0.5, set_name_size = 8, text_size = 9,  show_percentage = FALSE, 
                 fill_color = fillColors)
  
  ggsave(filename = fileName, plot = Venn, 
         height = 10, width = 10, units = "in", dpi = 300, bg = "white")
  
  Venn
  
}

# All patients
Venn_sig(sig_TNF, fillColors2,"Output_files/VennDiagram/venn_TNF_vst_q0.1_05.04.png")
Venn_sig(sig_aza, fillColors2,"Output_files/VennDiagram/venn_aza_vst_05.04.png")
Venn_sig(sig_mtx, fillColors2,"Output_files/VennDiagram/venn_mtx_vst_q0.1_05.04.png")
Venn_sig(sig_pred, fillColors2,"Output_files/VennDiagram/venn_pred_vst_05.04.png")

Venn_sig(sig_TNF_topVar, fillColors3,"Output_files/VennDiagram/venn_TNF_topVar_vst_q0.1_05.04.png")
Venn_sig(sig_aza_topVar, fillColors3,"Output_files/VennDiagram/venn_aza_topVar_vst_05.04.png")
Venn_sig(sig_mtx_topVar, fillColors3,"Output_files/VennDiagram/venn_mtx_topVar_vst_q0.1_05.04.png")
Venn_sig(sig_pred_topVar, fillColors3,"Output_files/VennDiagram/venn_pred_topVar_vst_05.04.png")

Venn_sig(sig_noBio_q0.05, fillColors3, "Output_files/VennDiagram/venn_TNF_topVar_vst_q0.05_05.04.png")
Venn_sig(sig_mtx_topVar_q0.05, fillColors3, "Output_files/VennDiagram/venn_mtx_topVar_vst_q0.05_05.04.png")



# Remitters only 
Venn_sig(R_sig_aza, fillColors2,"Output_files/VennDiagram/Remission/venn_aza_R_vst_q0.05_05.04.png")
Venn_sig(R_sig_mtx, fillColors2,"Output_files/VennDiagram/Remission/venn_mtx_R_vst_q0.01_05.04.png")
Venn_sig(R_sig_pred, fillColors2,"Output_files/VennDiagram/Remission/venn_pred_R_vst_q0.05_05.04.png")

Venn_sig(R_sig_TNF_topVar, fillColors3,"Output_files/VennDiagram/Remission/venn_TNF_topVar_R_vst_q0.05_05.04.png")
Venn_sig(R_sig_aza_topVar, fillColors3,"Output_files/VennDiagram/Remission/venn_aza_topVar_R_vst_q0.05_05.04.png")
Venn_sig(R_sig_mtx_topVar, fillColors3,"Output_files/VennDiagram/Remission/venn_mtx_topVar_R_vst_q0.01_05.04.png")
Venn_sig(R_sig_pred_topVar, fillColors3,"Output_files/VennDiagram/Remission/venn_pred_topVar_R_vst_q0.05_05.04.png")





# All patients x remitters only (antiTNF x no Bio, pred/aza/mtx x no Syst)
Venn_sig(sig_R_all_noBio, fillColors2,"Output_files/VennDiagram/All_vs_Remitters/venn_noBio_all_vs_R_vst_05.04.png")
Venn_sig(sig_R_all_aza, fillColors2,"Output_files/VennDiagram/All_vs_Remitters/venn_aza_all_vs_R_vst_05.04.png")
Venn_sig(sig_R_all_mtx, fillColors2,"Output_files/VennDiagram/All_vs_Remitters/venn_mtx_all_vs_R_vst_05.04.png")
Venn_sig(sig_R_all_pred, fillColors2,"Output_files/VennDiagram/All_vs_Remitters/venn_pred_all_vs_R_vst_05.04.png")

Venn_sig(sig_R_all_noBio_topVar, fillColors3,"Output_files/VennDiagram/All_vs_Remitters/venn_TNF_topVar_all_vs_R_vst_05.04.png")
Venn_sig(sig_R_all_aza_topVar, fillColors3,"Output_files/VennDiagram/All_vs_Remitters/venn_aza_topVar_all_vs_R_vst_05.04.png")
Venn_sig(sig_R_all_mtx_topVar, fillColors3,"Output_files/VennDiagram/All_vs_Remitters/venn_mtx_topVar_all_vs_R_vst_05.04.png")
Venn_sig(sig_R_all_pred_topVar, fillColors3,"Output_files/VennDiagram/All_vs_Remitters/venn_pred_topVar_all_vs_R_vst_05.04.png")


# Saving files ----
# setwd("C:/Documents/Masters thesis/EZE_cohort") # laptop
# setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC
# 
# # List with all significant genes
# sig_all <- list(`Anti-TNF vs Non-anti-TNF` = nonTNF,
#   `Anti-TNF vs No biologics` = noBio,
# 
#   `Aza vs No aza` = aza_noAza,
#   `Aza vs No systemic` = aza_noSyst,
# 
#   `MTX vs No MTX` = mtx_noMtx,
#   `MTX vs No systemic` = mtx_noSyst,
# 
#   `Pred vs No pred` = pred_noPred,
#   `Pred vs No systemic` = pred_noSyst)
# 
# list.save(sig_all, "Output_files/Maaslin2/significant_results/list_sig_genes_all_21.03.rds")
# 
# sig_all_R <- list(`Anti-TNF vs No biologics remitters` = R_noBio,
#                   
#                   `Aza vs No aza remitters` = R_aza_noAza,
#                   `Aza vs No systemic remitters` = R_aza_noSyst,
#                   
#                   # `MTX vs No MTX remitters` = R_mtx_noMtx,
#                   # `MTX vs No systemic remitters` = R_mtx_noSyst,
#                   
#                   `Pred vs No pred remitters` = R_pred_noPred,
#                   `Pred vs No systemic remitters` = R_pred_noSyst)
# 
# list.save(sig_all_R, "Output_files/Maaslin2/significant_results/list_sig_genes_remitters_21.03.rds")
# 
# # Intersection and diff of significant genes and top var genes 
# sig_all_intersect_diff <- list(
#   `intersect_aza` = intersect(aza_noAza, aza_noSyst),
#   `exclusive_noAza` = setdiff(aza_noAza, aza_noSyst),
#   `exclusive_aza_noSyst` = setdiff(aza_noSyst, aza_noAza),
#   
#   `intersect_pred` = intersect(pred_noPred, pred_noSyst),
#   `exclusive_noPred` = setdiff(pred_noPred, pred_noSyst),
#   `exclusive_pred_noSyst` = setdiff(pred_noSyst, pred_noPred),
#   
#   `intersect_mtx` = intersect(mtx_noMtx, mtx_noSyst),
#   `exclusive_noMtx` = setdiff(mtx_noMtx, mtx_noSyst),
#   `exclusive_pred_noSyst` = setdiff(mtx_noSyst, mtx_noMtx),
#   
#   `intersect_TNF` = intersect(nonTNF, noBio),
#   `exclusive_nonTNF` = setdiff(nonTNF, noBio),
#   `exclusive_noBio` = setdiff(noBio, nonTNF),
# 
#   `top_var_genes` = topVarGenes_names
# )
#   
# list.save(sig_all_intersect_diff, "Output_files/Maaslin2/significant_results/list_sig_genes_all_intersect_diff_21.03.rds")
# 
# intersect(intersect(a,b),c)

# # Extracting genes that are common to all patients, remitters and topvar (top2000)
# sig_common_genes <- list(
#   `intersect_aza_topVar` = intersect(intersect(aza_noSyst, R_aza_noSyst),topVarGenes_names),
#   `intersect_pred_topVar` = intersect(intersect(pred_noSyst, R_pred_noSyst), topVarGenes_names),
#   `intersect_noBio_topVar` = intersect(intersect(noBio_q0.05, R_noBio), topVarGenes_names)
# )

# All intersections and exclusive genes
sig_all_remitters_intersect_diff <- list(
  `intersect_aza` = intersect(aza_noSyst, R_aza_noSyst),
  `exclusive_all_aza` = setdiff(aza_noSyst, R_aza_noSyst),
  `exclusive_R_aza` = setdiff(R_aza_noSyst, aza_noSyst),
  `intersect_aza_topVar` = intersect(intersect(aza_noSyst, R_aza_noSyst),topVarGenes_names),
  `intersect_R_aza_topVar` = intersect(R_aza_noSyst, topVarGenes_names),

  `intersect_pred` = intersect(pred_noSyst, R_pred_noSyst),
  `exclusive_all_pred` = setdiff(pred_noSyst, R_pred_noSyst),
  `exclusive_R_pred` = setdiff(R_pred_noSyst, pred_noSyst),
  `intersect_pred_topVar` = intersect(intersect(pred_noSyst, R_pred_noSyst), topVarGenes_names),
  `intersect_R_pred_topVar` = intersect(R_pred_noSyst, topVarGenes_names),

  `intersect_noBio` = intersect(noBio_q0.05, R_noBio),
  `exclusive_all_noBio` = setdiff(noBio_q0.05, R_noBio),
  `exclusive_R_noBio` = setdiff(R_noBio, noBio_q0.05),
  `intersect_noBio_topVar` = intersect(intersect(noBio_q0.05, R_noBio), topVarGenes_names),
  `intersect_R_noBio_topVar` = intersect(R_noBio, topVarGenes_names),
  

  `top_var_genes` = topVarGenes_names
)

list.save(sig_all_remitters_intersect_diff, "Output_files/Maaslin2/significant_results/list_sig_genes_all_remitters_intersect_diff_06.04.rds")

