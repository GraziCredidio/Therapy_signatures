# EZE Cohort
# Aza(R+crp) with 0.5 L2FC threshold + TVGs


topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t") #TVG all patients
topVarGenes_names <- rownames(topVarGenes)

# All patients sig genes:
sig_genes_aza_all <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_aza_vs_noSyst_vst_q0.05.txt", sep = "\t")

aza_all_lfc0.5_up <- sig_genes_aza_all %>% 
  filter(coef > 0.5)
aza_all_lfc0.5_up <- aza_all_lfc0.5_up$feature

aza_all_lfc0.5_down <- sig_genes_aza_all %>% 
  filter(coef < -0.5)
aza_all_lfc0.5_down <- aza_all_lfc0.5_down$feature

# R_crp sig genes: 
sig_genes_aza <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_crp.txt", sep = "\t") 

aza_lfc05.5_up <- sig_genes_aza %>% 
  filter(coef > 0.5)
aza_lfc05.5_up <- aza_lfc05.5_up$feature

aza_lfc05.5_down <- sig_genes_aza %>% 
  filter(coef < -0.5)
aza_lfc05.5_down <- aza_lfc05.5_down$feature



# All x Remitters + crp  ----
aza_pos_R_crp_all_topVar <- list(`Aza vs No systemic (all)` = aza_all_lfc0.5_up,
                                  `Aza vs No systemic (inactive)` = aza_lfc05.5_up,
                                  `Top var genes` = topVarGenes_names)

aza_neg_R_crp_all_topVar <- list(`Aza vs No systemic (all)` = aza_all_lfc0.5_down,
                                  `Aza vs No systemic (inactive)` = aza_lfc05.5_down,
                                  `Top var genes` = topVarGenes_names)



# Venn diagrams ----
pos_fillColors <- c("#dfa693", "#dc6e55", "#bacb87")
neg_fillColors <- c("#8bbdd9", "#0079bf", "#bacb87")

Venn_sig <- function(data, fillColors, fileName) {
  Venn <- ggvenn(data = data, stroke_linetype = 2, stroke_size = 0.5, set_name_size = 8, text_size = 9,  show_percentage = FALSE, 
                 fill_color = fillColors)
  
  ggsave(filename = fileName, plot = Venn, 
         height = 10, width = 10, units = "in", dpi = 300, bg = "white")
  
  Venn
  
}


# All + remitters + TVGs 
Venn_sig(aza_pos_R_crp_all_topVar, pos_fillColors,"Output_files/VennDiagram/Remitters_crp/LFC0.5_allTVG/venn_aza_topVar_pos_lfc0.5.pdf")
Venn_sig(aza_neg_R_crp_all_topVar, neg_fillColors,"Output_files/VennDiagram/Remitters_crp/LFC0.5_allTVG/venn_aza_topVar_neg_lfc0.5.pdf")
