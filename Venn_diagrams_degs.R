# EZE cohort: Therapy Signatures
  # Venn diagrams DEGs (all patients x inactive patients x TVGs)
  # Figures 3A, 8A, 13A 
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(tidyverse)
library(ggvenn)

# Loading data ----
topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t") #TVG all patients
topVarGenes_names <- rownames(topVarGenes)

# Anti-TNF
degs_antiTNF <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_antiTNF_vs_noBio.txt", sep = '\t')
degs_antiTNF_all <- read.table("Output_files/Maaslin2/significant_results/allPatients/maaslin2_significant_results_antiTNF_vs_noBio_allPatients.txt", sep = '\t')

# Pred
degs_pred <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_pred_vs_noSyst.txt", sep = '\t')
degs_pred_all <- read.table("Output_files/Maaslin2/significant_results/allPatients/maaslin2_significant_results_pred_vs_noSyst_allPatients.txt", sep = '\t')
  
# Aza
degs_aza <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_aza_vs_noSyst.txt", sep = '\t')
degs_aza_all <- read.table("Output_files/Maaslin2/significant_results/allPatients/maaslin2_significant_results_aza_vs_noSyst_allPatients.txt", sep = '\t')

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
degs_antiTNF_up_all <- filter_degs(degs_antiTNF_all, '>')

degs_antiTNF_down <- filter_degs(degs_antiTNF, '<')
degs_antiTNF_down_all <- filter_degs(degs_antiTNF_all, '<')

# Azathioprine
degs_aza_up <- filter_degs(degs_aza, '>')
degs_aza_up_all <- filter_degs(degs_aza_all, '>')

degs_aza_down <- filter_degs(degs_aza, '<')
degs_aza_down_all <- filter_degs(degs_aza_all, '<')

# Prednisolone: absLFC > 0.5
degs_pred_up <- filter_degs(degs_pred, '>', 0.5)
degs_pred_up_all <- filter_degs(degs_pred_all, '>', 0.5)

degs_pred_down <- filter_degs(degs_pred, '<', -0.5)
degs_pred_down_all <- filter_degs(degs_pred_all, '<', -0.5)


# Creating lists with up and downregulated degs----
# Anti-TNF
list_antiTNF_up <- list(`All` = degs_antiTNF_up_all,
                     `Inactive` = degs_antiTNF_up,
                     `TVGs` = topVarGenes_names)
list_antiTNF_down <- list(`All` = degs_antiTNF_down_all,
                       `Inactive` = degs_antiTNF_down,
                       `TVGs` = topVarGenes_names)

# Azathioprine
list_aza_up <- list(`All` = degs_aza_up_all,
                     `Inactive` = degs_aza_up,
                     `TVGs` = topVarGenes_names)
list_aza_down <- list(`All` = degs_aza_down_all,
                       `Inactive` = degs_aza_down,
                       `TVGs` = topVarGenes_names)

# Prednisolone
list_pred_up <- list(`All` = degs_pred_up_all,
                     `Inactive` = degs_pred_up,
                     `TVGs` = topVarGenes_names)
list_pred_down <- list(`All` = degs_pred_down_all,
                       `Inactive` = degs_pred_down,
                       `TVGs` = topVarGenes_names)

# Venn diagrams ----
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
dir.create("Output_files/Venn_diagrams/degs/antiTNF")
venn_degs(list_antiTNF_up, up_fillColors, "Output_files/Venn_diagrams/degs/antiTNF/venn_up_antiTNF.pdf")
venn_degs(list_antiTNF_down, down_fillColors, "Output_files/Venn_diagrams/degs/antiTNF/venn_down_antiTNF.pdf")

# Aza
dir.create("Output_files/Venn_diagrams/degs/aza")
venn_degs(list_aza_up, up_fillColors, "Output_files/Venn_diagrams/degs/aza/venn_up_aza.pdf")
venn_degs(list_aza_down, down_fillColors, "Output_files/Venn_diagrams/degs/aza/venn_down_aza.pdf")

# Pred
dir.create("Output_files/Venn_diagrams/degs/pred")
venn_degs(list_pred_up, up_fillColors, "Output_files/Venn_diagrams/degs/pred/venn_up_pred.pdf")
venn_degs(list_pred_down, down_fillColors, "Output_files/Venn_diagrams/degs/pred/venn_down_pred.pdf")
