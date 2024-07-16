# EZE cohort: Therapy Signatures
  # DNA methylation-linked DEGs visualizations
  # Density plots: Figures 6D, 12D, 19B
  # Author: Graziella Credidio

rm(list = ls())

folder <- "Output_files/Methylation/antiTNF/GO"
if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading packages ----
library(tidyverse)
library(ggridges)

# Loading data ----
# Azathioprine
aza_dnam_degs_up <- read.table("Output_files/Methylation/aza/aza_upregulated_DNAm_DEG.txt", sep = '\t')
aza_dnam_degs_down <- read.table("Output_files/Methylation/aza/aza_downregulated_DNAm_DEG.txt", sep = '\t')

aza_go_up <- read.table("Output_files/Methylation/aza/GO/GO_aza_TVG_DNAm_up_DEGs.txt", sep = '\t')
aza_go_down <- read.table("Output_files/Methylation/aza/GO/GO_aza_TVG_DNAm_down_DEGs.txt", sep = '\t')

# Prednisolone
pred_dnam_degs_up <- read.table("Output_files/Methylation/pred/pred_upregulated_DNAm_DEG.txt", sep = '\t')
pred_dnam_degs_down <- read.table("Output_files/Methylation/pred/pred_downregulated_DNAm_DEG.txt", sep = '\t')

pred_go_up <- read.table("Output_files/Methylation/pred/GO/GO_pred_TVG_DNAm_up_DEGs.txt", sep = '\t')
pred_go_down <- read.table("Output_files/Methylation/pred/GO/GO_pred_TVG_DNAm_down_DEGs.txt", sep = '\t')

# Anti-TNF
antiTNF_dnam_degs_up <- read.table("Output_files/Methylation/antiTNF/antiTNF_upregulated_DNAm_DEG.txt", sep = '\t')
antiTNF_dnam_degs_down <- read.table("Output_files/Methylation/antiTNF/antiTNF_downregulated_DNAm_DEG.txt", sep = '\t')

antiTNF_go_up <- read.table("Output_files/Methylation/antiTNF/GO/GO_antiTNF_DNAm_up_DEGs.txt", sep = '\t')
antiTNF_go_down <- read.table("Output_files/Methylation/antiTNF/GO/GO_antiTNF_DNAm_down_DEGs.txt", sep = '\t')

# Selected GO terms ----
aza_terms <- c("natural killer cell mediated cytotoxicity", "positive regulation of lymphocyte mediated immunity",
               "G protein-coupled receptor signaling pathway", "positive regulation of cytokine production involved in immune response",
               "cell surface receptor signaling pathway", "positive regulation of programmed cell death",
               "positive regulation of innate immune response", "immune response", "response to type II interferon",
               "cell killing", "cellular response to tumor necrosis factor", "regulation of natural killer cell chemotaxis", 
               "positive regulation of T-helper 1 cell cytokine production","T cell mediated cytotoxicity", "monocyte chemotaxis", 
               "positive regulation of cell population proliferation","angiogenesis", "inorganic ion homeostasis", "tissue development", 
               "regulation of anatomical structure morphogenesis","phosphorylation", "monoatomic cation homeostasis", 
               "erythrocyte differentiation", "ammonium transmembrane transport", "gas transport")

pred_terms <- c("positive regulation of inflammatory response", "defense response to bacterium", "acute inflammatory response", 
                "positive regulation of cytokine production", "regulation of macrophage activation", "inflammatory response",
                "positive regulation of NF-kappaB transcription factor activity","neutrophil extravasation", 
                "regulation of toll-like receptor signaling pathway", "regulation of lipid storage", 
                "cell surface receptor signaling pathway", "B cell receptor signaling pathway",
                "positive regulation of nitric-oxide synthase biosynthetic process", "negative regulation of cytokine production",
                "interleukin-18-mediated signaling pathway", "production of molecular mediator involved in inflammatory response",
                "tumor necrosis factor production", "cytokine production involved in inflammatory response", "lipid homeostasis", 
                "regulation of long-chain fatty acid import across plasma membrane","B cell proliferation", "B cell differentiation", 
                "cell adhesion", "purine deoxyribonucleotide biosynthetic process")

antiTNF_terms <- c("immunoglobulin production", "adaptive immune response", "cartilage development", "regulation of cell cycle G1/S phase transition",
                   "interleukin-1-mediated signaling pathway", "positive regulation of chromosome condensation",
                   "DNA replication-dependent chromatin assembly", "response to osmotic stress", "positive regulation of signal transduction", 
                   "regulation of mRNA stability")


# Creating dataframes of selected GO terms
go_selected_df <- function(go_data, dir, terms){
  go_selected <- go_data %>% 
    mutate(plot = ifelse(Term %in% terms, "Yes", "No")) %>% 
    filter(plot == "Yes") %>% 
    mutate(direction = dir)
  go_selected$Term <- Term(as.vector(go_selected$GO.ID))
  go_selected$Term <- factor(go_selected$Term, levels = rev(go_selected$Term))
  go_selected$Genes <- sub('.', '', go_selected$Genes)
  
  return(go_selected)
}

aza_go_selected_up <- go_selected_df(aza_go_up, "upregulated", aza_terms)
aza_go_selected_down <- go_selected_df(aza_go_down, "downregulated", aza_terms)

pred_go_selected_up <- go_selected_df(pred_go_up, "upregulated", pred_terms)
pred_go_selected_down <- go_selected_df(pred_go_down, "downregulated", pred_terms)

antiTNF_go_selected_up <- go_selected_df(antiTNF_go_up, "upregulated", antiTNF_terms)
antiTNF_go_selected_down <- go_selected_df(antiTNF_go_down, "downregulated", antiTNF_terms)

# Creating data frames to be plotted ----
transform_data <- function(go_selected, dnam_degs, reg) {
  GO_plot_df <- data.frame(Gene = character(),
                           Term = character(),
                           Rho = character(),
                           FDR = character(),
                           Regulation = character(),
                           stringsAsFactors = FALSE)
  cor_genes <- dnam_degs
  rownum <- 1
  for (i in 1:nrow(go_selected)) {
    individual_genes <- go_selected[i, 10]
    individual_genes_vector <- unlist(strsplit(individual_genes, ","))
    
    for (j in 1:length(individual_genes_vector)) {
      Gene <- individual_genes_vector[j]
      Term <- as.character(go_selected[i, 2])
      Rho <- cor_genes[cor_genes$Gene_name %in% Gene, 4]
      FDR <- cor_genes[cor_genes$Gene_name %in% Gene, 5]
      Regulation <- reg
      
      if(Gene %in% cor_genes$Gene_name) { #exclude genes that are not in cor meth data
        Rho <- cor_genes[cor_genes$Gene_name %in% Gene, 4]
        FDR <- cor_genes[cor_genes$Gene_name %in% Gene, 5]
        Regulation <- reg
        
        if(length(Rho) > 1) {
          for(m in 1:length(Rho)) {
            GO_plot_df[rownum, ] <- c(as.character(Gene), Term, Rho[m], FDR[m], Regulation)
            rownum <- rownum + 1
          }
        } else {
          GO_plot_df[rownum, ] <- c(as.character(Gene), Term, Rho, FDR, Regulation)
          rownum <- rownum + 1
        }
      }
    }
  }
  return(GO_plot_df)
}


aza_plot_down_df <- transform_data(aza_go_selected_down, aza_dnam_degs_down, "downregulated")
aza_plot_up_df <- transform_data(aza_go_selected_up, aza_dnam_degs_up, "upregulated")

pred_plot_down_df <- transform_data(pred_go_selected_down, pred_dnam_degs_down, "downregulated")
pred_plot_up_df <- transform_data(pred_go_selected_up, pred_dnam_degs_up, "upregulated")

antiTNF_plot_down_df <- transform_data(antiTNF_go_selected_down, antiTNF_dnam_degs_down, "downregulated")
antiTNF_plot_up_df <- transform_data(antiTNF_go_selected_up, antiTNF_dnam_degs_up, "upregulated")

# Binding up and down dataframees
df_preprocessing <- function(df1, df2){
  plot_df <- rbind(df1, df2)
  plot_df$Rho <- as.numeric(plot_df$Rho)
  plot_df$FDR <- as.numeric(plot_df$FDR)
  plot_df$Term <- factor(plot_df$Term, levels = rev(unique(plot_df$Term))) 
  return(plot_df)
}

aza_plot_df <- df_preprocessing(aza_plot_down_df, aza_plot_up_df)
pred_plot_df <- df_preprocessing(pred_plot_down_df, pred_plot_up_df)
antiTNF_plot_df <- df_preprocessing(antiTNF_plot_down_df, antiTNF_plot_up_df)

# Ridgeplots: Aza and Pred ----
ridge_plot <- function(df, filePath){
  plot <- ggplot(pred_plot_df, aes(Rho, Term, fill = Regulation)) +
    scale_fill_manual("", values = c("upregulated" = "#d91f24", "downregulated" = "#365196")) +
    geom_density_ridges(alpha = .6) +
    scale_size_continuous(name="FDR", range = c(0.2,4), trans = 'reverse') + 
    xlim(-1,1) +
    theme_ridges() +
    ylab("") + xlab("Spearman's Rho") +
    theme(legend.position = "bottom",
          axis.text.y = element_text(angle = 0, size = 15, vjust = 0.5),
          legend.text = element_text(size = 15))
  ggsave(plot, file = filePath, height = 10, width = 13, units = "in")
  plot
  
}

ridge_plot(aza_plot_df, "Output_files/Methylation/aza/GO/GO_aza_ridgeplot_DNAm_DEGs.pdf")
ridge_plot(pred_plot_df, "Output_files/Methylation/pred/GO/GO_pred_ridgeplot_DNAm_DEGs.pdf")


# Density dot plot: Anti-TNF ----
density_dot_plot <- ggplot(antiTNF_plot_df, aes(Rho,Term, size = FDR, fill = Regulation)) +
  geom_point(pch = 21, colour = "black", stroke = 0.9, position = "jitter") +
  scale_size_continuous(name="FDR", range = c(0.2,4), trans = 'reverse') + 
  scale_fill_manual("", values = c("upregulated" = "#d91f24", "downregulated" = "#365196")) +
  xlim(-1,1) +
  geom_hline(yintercept = seq(1.5, 12.5, by = 1)) +
  theme_minimal() +
  ylab("") + xlab("Spearman's Rho") +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(angle = 0, size = 15, vjust = 0.5),
        legend.text = element_text(size = 15))
density_dot_plot
ggsave(density_dot_plot, file = "Output_files/Methylation/antiTNF/GO/GO_antiTNF_density_DNAm_DEGs.pdf", height = 6, width = 13, units = "in")