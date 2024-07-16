# EZE cohort: Therapy Signatures
  # Venn diagrams DNAm-DEGs (inactive DEGs x DNAm-DEG x TVGs)
  # Figures 6A, 12A, 17A 
  # Author: Graziella Credidio

rm(list = ls())

folders <- c("Output_files/Venn_diagrams/DNAm_degs/aza", "Output_files/Venn_diagrams/DNAm_degs/pred", "Output_files/Venn_diagrams/DNAm_degs/antiTNF")
for (i in folders){
  if (!dir.exists(i)) {
    dir.create(i)}
}

# Loading packages ----
library(tidyverse)
library(ggvenn)

# Loading data ----
topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t") #TVG all patients
topVarGenes_names <- rownames(topVarGenes)

# Anti-TNF
degs_antiTNF <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_antiTNF_vs_noBio.txt", sep = '\t')

up_dnam_degs_antiTNF <- read.table("Output_files/Methylation/antiTNF/antiTNF_upregulated_DNAm_DEG.txt", sep = '\t')
up_dnam_degs_antiTNF <- unique(up_dnam_degs_antiTNF$Gene)

down_dnam_degs_antiTNF <- read.table("Output_files/Methylation/antiTNF/antiTNF_downregulated_DNAm_DEG.txt", sep = '\t')
down_dnam_degs_antiTNF <- unique(down_dnam_degs_antiTNF$Gene)

# Pred
degs_pred <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_pred_vs_noSyst.txt", sep = '\t')

up_dnam_degs_pred <- read.table("Output_files/Methylation/pred/pred_upregulated_DNAm_DEG.txt", sep = '\t')
up_dnam_degs_pred <- unique(up_dnam_degs_pred$Gene)

down_dnam_degs_pred <- read.table("Output_files/Methylation/pred/pred_downregulated_DNAm_DEG.txt", sep = '\t')
down_dnam_degs_pred <- unique(down_dnam_degs_pred$Gene)

# Aza
degs_aza <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_aza_vs_noSyst.txt", sep = '\t')

up_dnam_degs_aza <- read.table("Output_files/Methylation/aza/aza_upregulated_DNAm_DEG.txt", sep = '\t')
up_dnam_degs_aza <- unique(up_dnam_degs_aza$Gene)

down_dnam_degs_aza <- read.table("Output_files/Methylation/aza/aza_downregulated_DNAm_DEG.txt", sep = '\t')
down_dnam_degs_aza <- unique(down_dnam_degs_aza$Gene)

# Data preprocessing ----
filter_degs <- function(degs_df, greater_lower_than, filter_coef = 0){
  if (greater_lower_than == '>'){
    degs_filtered <- degs_df %>% 
      filter(coef > filter_coef) %>% 
      pull(feature)
  } 
  else{
    degs_filtered <- degs_df %>% 
      filter(coef < filter_coef) %>% 
      pull(feature)
  }
  
  return(degs_filtered)
}

# anti-TNF
degs_antiTNF_up <- filter_degs(degs_antiTNF, '>')
degs_antiTNF_down <- filter_degs(degs_antiTNF, '<')

# Azathioprine
degs_aza_up <- filter_degs(degs_aza, '>')
degs_aza_down <- filter_degs(degs_aza, '<')

# Prednisolone: absLFC > 0.5
degs_pred_up <- filter_degs(degs_pred, '>', 0.5)
degs_pred_down <- filter_degs(degs_pred, '<', -0.5)


# Creating lists with up and downregulated gene sets ----
# Anti-TNF
list_antiTNF_up <- list(`DEGs` = degs_antiTNF_up,
                        `DNAm-DEGs` = up_dnam_degs_antiTNF,
                        `TVGs` = topVarGenes_names)
list_antiTNF_down <- list(`DEGs` = degs_antiTNF_down,
                          `DNAm-DEGs` = down_dnam_degs_antiTNF,
                          `TVGs` = topVarGenes_names)

# Azathioprine
list_aza_up <- list(`DEGs` = degs_aza_up,
                    `DNAm-DEGs` = up_dnam_degs_aza,
                    `TVGs` = topVarGenes_names)
list_aza_down <- list(`DEGs` = degs_aza_down,
                      `DNAm-DEGs` = down_dnam_degs_aza,
                      `TVGs` = topVarGenes_names)

# Prednisolone
list_pred_up <- list(`DEGs` = degs_pred_up,
                     `DNAm-DEGs` = up_dnam_degs_pred,
                     `TVGs` = topVarGenes_names)
list_pred_down <- list(`DEGs` = degs_pred_down,
                       `DNAm-DEGs` = down_dnam_degs_pred,
                       `TVGs` = topVarGenes_names)

# Venn Diagrams ----
up_fillColors <- c("#dfa693", "#dc6e55", "#bacb87")
down_fillColors <- c("#8bbdd9", "#0079bf", "#bacb87")

venn_degs <- function(data, fillColors, filePath) {
  Venn <- ggvenn(data = data, stroke_linetype = 2, stroke_size = 0.5, set_name_size = 8, text_size = 9,  show_percentage = FALSE, 
                 fill_color = fillColors)
  ggsave(filename = filePath, plot = Venn, 
         height = 10, width = 10, units = "in", dpi = 300, bg = "white")
  Venn
}

# Creating venn diagrams
# Anti-TNF
dir.create("Output_files/Venn_diagrams/DNAm_degs/antiTNF")
venn_degs(list_antiTNF_up, up_fillColors, "Output_files/Venn_diagrams/DNAm_degs/antiTNF/venn_DNAm_up_antiTNF.pdf")
venn_degs(list_antiTNF_down, down_fillColors, "Output_files/Venn_diagrams/DNAm_degs/antiTNF/venn_DNAm_down_antiTNF.pdf")

# Aza
dir.create("Output_files/Venn_diagrams/DNAm_degs/aza")
venn_degs(list_aza_up, up_fillColors, "Output_files/Venn_diagrams/DNAm_degs/aza/venn_DNAm_up_aza.pdf")
venn_degs(list_aza_down, down_fillColors, "Output_files/Venn_diagrams/DNAm_degs/aza/venn_DNAm_down_aza.pdf")

# Pred
dir.create("Output_files/Venn_diagrams/DNAm_degs/pred")
venn_degs(list_pred_up, up_fillColors, "Output_files/Venn_diagrams/DNAm_degs/pred/venn_DNAm_up_pred.pdf")
venn_degs(list_pred_down, down_fillColors, "Output_files/Venn_diagrams/DNAm_degs/pred/venn_DNAm_down_pred.pdf")
