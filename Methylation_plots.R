# EZE cohort
# Methylation-transcriptome integration analysis
# Distance from TSSS and Rho value

graphics.off()
rm(list = ls())

# Methylation linked genes: direction and position plots ----
cor_genes_aza <- read.table("Output_files/Methylation/aza/aza_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = "\t")
cor_genes_pred <- read.table("Output_files/Methylation/pred/pred_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = "\t")
cor_genes_antiTNF <- read.table("Output_files/Methylation/antiTNF/antiTNF_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = "\t")

# Direction of correlation analysis
cor_genes_aza <- cor_genes_aza %>% 
  mutate(Direction = ifelse(Rho > 0, "Rho > 0", "Rho < 0")) %>% 
  mutate(Position = ifelse(Distance_from_TSS > 0, "Gene", "Promoter")) %>% 
  mutate(Analysis = "Aza x No syst")

cor_genes_pred <- cor_genes_pred %>% 
  mutate(Direction = ifelse(Rho > 0, "Rho > 0", "Rho < 0")) %>% 
  mutate(Position = ifelse(Distance_from_TSS > 0, "Gene", "Promoter")) %>% 
  mutate(Analysis = "Pred x No syst")

cor_genes_antiTNF <- cor_genes_antiTNF %>% 
  mutate(Direction = ifelse(Rho > 0, "Rho > 0", "Rho < 0")) %>% 
  mutate(Position = ifelse(Distance_from_TSS > 0, "Gene", "Promoter")) %>% 
  mutate(Analysis = "AntiTNF x No biologics")


# Plots
methylation_plot_single <- function(data, subtitle) {
  ggplot(data, aes(Position, fill = Direction)) +
    geom_bar(width = 0.3) +
    coord_flip() +
    geom_text(aes(label = after_stat(count)), 
              stat = "count", hjust = 0.5, position = position_stack(vjust = 0.8), color = "white", size = 3) +
    scale_fill_manual("", values = c ("Rho > 0" = "#d91f24", "Rho < 0" = "#365196")) + 
    theme_classic() +
    ylab("Number of sites") + xlab(" ") +
    labs(title = "DNA methylation-linked genes",
         subtitle = subtitle,
         fill = "") +
    theme(plot.title = element_text(angle = 0, size = 15, face = 'bold', vjust = 1, margin = margin(0,0,10,0)),
          plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(2,0,10,0)),
          axis.title.x = element_text(vjust = -1, size = 12),
          axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 12,),
          legend.position = "bottom")
}

methylation_plot_single(cor_genes_aza, "Azathioprine x No syst")
methylation_plot_single(cor_genes_pred, "Prednisolone x No syst")
methylation_plot_single(cor_genes_antiTNF, "Anti-TNF x No Biologics")




#### Methylation linked genes: position and direction pltos splitted by up and downregulated genes ----
# Methylation linked genes: direction and position plots 
aza_up_cor_genes <- read.table("Output_files/Methylation/aza/aza_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')
aza_up_cor_genes$Regulation <- "Upregulated"
aza_down_cor_genes <- read.table("Output_files/Methylation/aza/aza_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')
aza_down_cor_genes$Regulation <- "Downregulated"
cor_genes_aza <- rbind(aza_up_cor_genes, aza_down_cor_genes)

pred_up_cor_genes <- read.table("Output_files/Methylation/pred/pred_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')
pred_up_cor_genes$Regulation <- "Upregulated"
pred_down_cor_genes <- read.table("Output_files/Methylation/pred/pred_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')
pred_down_cor_genes$Regulation <- "Downregulated"
cor_genes_pred <- rbind(pred_up_cor_genes, pred_down_cor_genes)

antiTNF_up_cor_genes <- read.table("Output_files/Methylation/antiTNF/antiTNF_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')
antiTNF_up_cor_genes$Regulation <- "Upregulated"
antiTNF_down_cor_genes <- read.table("Output_files/Methylation/antiTNF/antiTNF_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')
antiTNF_down_cor_genes$Regulation <- "Downregulated"
cor_genes_antiTNF <- rbind(antiTNF_up_cor_genes, antiTNF_down_cor_genes)

# Direction of correlation analysis
cor_genes_aza <- cor_genes_aza %>% 
  mutate(Direction = ifelse(Rho > 0, "Rho > 0", "Rho < 0")) %>% 
  mutate(Position = ifelse(Distance_from_TSS > 0, "Gene", "Promoter")) %>% 
  mutate(Analysis = "Aza x No syst")

cor_genes_pred <- cor_genes_pred %>% 
  mutate(Direction = ifelse(Rho > 0, "Rho > 0", "Rho < 0")) %>% 
  mutate(Position = ifelse(Distance_from_TSS > 0, "Gene", "Promoter")) %>% 
  mutate(Analysis = "Pred x No syst")

cor_genes_antiTNF <- cor_genes_antiTNF %>% 
  mutate(Direction = ifelse(Rho > 0, "Rho > 0", "Rho < 0")) %>% 
  mutate(Position = ifelse(Distance_from_TSS > 0, "Gene", "Promoter")) %>% 
  mutate(Analysis = "AntiTNF x No biologics")


# Plots
methylation_plot_single <- function(data, subtitle) {
  ggplot(data, aes(Position, fill = Direction)) +
    geom_bar(width = 0.3) +
    coord_flip() +
    geom_text(aes(label = after_stat(count)), 
              stat = "count", hjust = 0.5, position = position_stack(vjust = 0.8), color = "white", size = 3) +
    scale_fill_manual("", values = c ("Rho > 0" = "#d91f24", "Rho < 0" = "#365196")) + 
    theme_classic() +
    ylab("Number of sites") + xlab(" ") +
    labs(title = "DNA methylation-linked genes",
         subtitle = subtitle,
         fill = "") +
    theme(plot.title = element_text(angle = 0, size = 15, face = 'bold', vjust = 1, margin = margin(0,0,10,0)),
          plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(2,0,10,0)),
          axis.title.x = element_text(vjust = -1, size = 12),
          axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 12,),
          legend.position = "bottom") +
    facet_grid(~Regulation)
}

methylation_plot_single(cor_genes_aza, "Azathioprine x No syst")
methylation_plot_single(cor_genes_pred, "Prednisolone x No syst")
methylation_plot_single(cor_genes_antiTNF, "Anti-TNF x No Biologics")







##### Biological processes, FDR and Rho plot ----
# AZA (DMLGs + topvar)
rm(list = ls())
aza_up_cor_genes <- read.table("Output_files/Methylation/aza/aza_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')
aza_down_cor_genes <- read.table("Output_files/Methylation/aza/aza_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')

aza_up_go <- read.table("Output_files/Methylation/aza/GO/correlated_topVar/aza_correlated_topVar_up_DEGs_GO.txt", sep = '\t')
aza_down_go <- read.table("Output_files/Methylation/aza/GO/correlated_topVar/aza_correlated_topVar_down_DEGs_GO.txt", sep = '\t')

# GO terms

terms <- c("natural killer cell mediated cytotoxicity", "positive regulation of lymphocyte mediated immunity",
               "G protein-coupled receptor signaling pathway", "positive regulation of cytokine production involved in immune response",
               "cell surface receptor signaling pathway", "positive regulation of programmed cell death",
               "positive regulation of innate immune response", "immune response", "response to type II interferon",
               "cell killing", "cellular response to tumor necrosis factor", # same down
           
           "regulation of natural killer cell chemotaxis", "positive regulation of T-helper 1 cell cytokine production",
           "T cell mediated cytotoxicity", "monocyte chemotaxis", "positive regulation of cell population proliferation", #new down
               
               "angiogenesis", "inorganic ion homeostasis", "tissue development", "regulation of anatomical structure morphogenesis",
               "phosphorylation", "monoatomic cation homeostasis", # same up
              
               "erythrocyte differentiation", "ammonium transmembrane transport", "gas transport" #new up
      )


plot_aza_up_go <- aza_up_go %>% 
  mutate(plot = ifelse(Term %in% terms, "Yes", "No")) %>% 
  filter(plot == "Yes") %>% 
  mutate(direction = "Upregulated")
plot_aza_up_go$Term <- Term(as.vector(plot_aza_up_go$GO.ID))
plot_aza_up_go$Term <- factor(plot_aza_up_go$Term, levels = rev(plot_aza_up_go$Term))
plot_aza_up_go$Genes <- sub('.', '', plot_aza_up_go$Genes)


plot_aza_down_go <- aza_down_go %>% 
  mutate(plot = ifelse(Term %in% terms, "Yes", "No")) %>% 
  filter(plot == "Yes") %>% 
  mutate(direction = "Downregulated")
plot_aza_down_go$Term <- Term(as.vector(plot_aza_down_go$GO.ID))
plot_aza_down_go$Term <- factor(plot_aza_down_go$Term, levels = rev(plot_aza_down_go$Term))
plot_aza_down_go$Genes <- sub('.', '', plot_aza_down_go$Genes)

# Creating data frame to be plotted
#UP
GO_plot_up <- data.frame(Gene = "",
                      Term = "",
                      Rho = "",
                      FDR = "",
                      Regulation = "")
cor_genes <- aza_up_cor_genes
rownum = 1

for (i in 1:nrow(plot_aza_up_go)){
  individual_genes <- plot_aza_up_go[i,10]
  individual_genes_vector <-  unlist(strsplit(individual_genes, ","))
  
  for (j in 1:length(individual_genes_vector)){
    Gene <- individual_genes_vector[j]
    Term <- as.character(plot_aza_up_go[i, 2])
    Rho <- cor_genes[cor_genes$Gene_name %in% Gene, 4]
    FDR <- cor_genes[cor_genes$Gene_name %in% Gene, 5]
    Regulation <- "upregulated"
    
    if(length(Rho) > 1){
      
      for(m in 1:length(Rho)){
        GO_plot_up[rownum, ] <- c(as.character(Gene), Term, Rho[m], FDR[m], Regulation)
        rownum <- rownum + 1
      }
      
    } else {
      GO_plot_up[rownum, ] <- c(as.character(Gene), Term, Rho, FDR, Regulation)
      rownum <- rownum + 1
  
    }
    
  }
}



# DOWN
GO_plot_down <- data.frame(Gene = "",
                         Term = "",
                         Rho = "",
                         FDR = "",
                         Regulation = "")
cor_genes <- aza_down_cor_genes
rownum = 1

for (i in 1:nrow(plot_aza_down_go)){
  individual_genes <- plot_aza_down_go[i,10]
  individual_genes_vector <-  unlist(strsplit(individual_genes, ","))
  
  for (j in 1:length(individual_genes_vector)){
    Gene <- individual_genes_vector[j]
    Term <- as.character(plot_aza_down_go[i, 2])
    Rho <- cor_genes[cor_genes$Gene_name %in% Gene, 4]
    FDR <- cor_genes[cor_genes$Gene_name %in% Gene, 5]
    Regulation <- "downregulated"
    
    if(length(Rho) > 1){
      
      for(m in 1:length(Rho)){
        GO_plot_down[rownum, ] <- c(as.character(Gene), Term, Rho[m], FDR[m], Regulation)
        rownum <- rownum + 1
      }
      
    } else {
      GO_plot_down[rownum, ] <- c(as.character(Gene), Term, Rho, FDR, Regulation)
      rownum <- rownum + 1
      
    }
    
  }
}

GO_plot <- rbind(GO_plot_up, GO_plot_down)
GO_plot$Rho <- as.numeric(GO_plot$Rho)
GO_plot$FDR <- as.numeric(GO_plot$FDR)
GO_plot$Term <- factor(GO_plot$Term, levels = rev(unique(GO_plot$Term))) # fixes order

write.table(GO_plot, "Output_files/Methylation/aza/GO/correlated_topVar/GO_plot_aza_DMLGs_topVar_sameTerms.txt", sep = "\t")
GO_plot <- read.table("Output_files/Methylation/aza/GO/correlated_topVar/GO_plot_aza_DMLGs_topVar_sameTerms.txt", sep = "\t")

# Plot 
pdf("Output_files/Methylation/aza/GO/correlated_topVar/dotplot_GO_Rho_aza_DMLGs_topVar_sameTerms.pdf", width = 13, height = 10)
dot_plot <- ggplot(GO_plot, aes(Rho,Term, size = FDR, fill = Regulation)) +
  geom_point(pch = 21, colour = "black", stroke = 0.9, position = "jitter") +
  scale_size_continuous(name="FDR", range = c(0.2,4), trans = 'reverse') + #reverse size of dots proportional to FDR)+
  scale_fill_manual("", values = c("upregulated" = "#d91f24", "downregulated" = "#365196")) +
  xlim(-1,1) +
  geom_hline(yintercept = seq(1.5, 14.5, by = 1)) +
  theme_minimal() +
  ylab("") + xlab("Spearman's Rho") +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        
        axis.text.y = element_text(angle = 0, size = 15, vjust = 0.5),
        legend.text = element_text(size = 15))

dot_plot
dev.off()

# Ridgeplot
pdf("Output_files/Methylation/aza/GO/correlated_topVar/ridgeplot_GO_Rho_aza_DMLGs_topVar_sameTerms.pdf", width = 13, height = 10)
ridge_plot <- ggplot(GO_plot, aes(Rho,Term, fill = Regulation)) +
  scale_fill_manual("", values = c("upregulated" = "#d91f24", "downregulated" = "#365196")) +
  geom_density_ridges(alpha = .6) +
  scale_size_continuous(name="FDR", range = c(0.2,4), trans = 'reverse') + #reverse size of dots proportional to FDR)+
  xlim(-1,1) +
  theme_ridges() +
  ylab("") + xlab("Spearman's Rho") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(angle = 0, size = 15, vjust = 0.5),
        legend.text = element_text(size = 15))
ridge_plot
dev.off()

 ##################################################
# PRED
rm(list = ls())

pred_up_cor_genes <- read.table("Output_files/Methylation/pred/pred_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')
pred_down_cor_genes <- read.table("Output_files/Methylation/pred/pred_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')

up_go_result <- read.csv("Output_files/Methylation/pred/GO/correlated_topvar/lfc0.5/pred_lfc0.5_correlated_topVar_up_DEGs_GO.txt", sep = '\t') #lfc 0.5 + topvar
down_go_result <- read.csv("Output_files/Methylation/pred/GO/correlated_topvar/lfc0.5/pred_lfc0.5_correlated_topVar_down_DEGs_GO.txt", sep = '\t') #lfc 0.5 + topvar

#new_terms + same terms
terms <- c("B cell receptor signaling pathway", #same neg
           "B cell proliferation", "B cell differentiation", "cell adhesion", "purine deoxyribonucleotide biosynthetic process", #new neg
           
           "cell surface receptor signaling pathway", #=
           
           "positive regulation of cytokine production", "regulation of macrophage activation", "inflammatory response",
           "positive regulation of NF-kappaB transcription factor activity",
           "neutrophil extravasation", "regulation of toll-like receptor signaling pathway", #same pos
           "positive regulation of nitric-oxide synthase biosynthetic process", "negative regulation of cytokine production",
           "interleukin-18-mediated signaling pathway", "production of molecular mediator involved in inflammatory response",
           "tumor necrosis factor production", "cytokine production involved in inflammatory response", "lipid homeostasis",
           "regulation of long-chain fatty acid import across plasma membrane")


plot_pred_up_go <- up_go_result %>% 
  mutate(plot = ifelse(Term %in% terms, "Yes", "No")) %>% 
  filter(plot == "Yes") %>% 
  mutate(direction = "Upregulated")
plot_pred_up_go$Term <- Term(as.vector(plot_pred_up_go$GO.ID))
plot_pred_up_go$Term <- factor(plot_pred_up_go$Term, levels = rev(plot_pred_up_go$Term))
plot_pred_up_go$Genes <- sub('.', '', plot_pred_up_go$Genes)


plot_pred_down_go <- down_go_result %>% 
  mutate(plot = ifelse(Term %in% terms, "Yes", "No")) %>% 
  filter(plot == "Yes") %>% 
  mutate(direction = "Downregulated")
plot_pred_down_go$Term <- Term(as.vector(plot_pred_down_go$GO.ID))
plot_pred_down_go$Term <- factor(plot_pred_down_go$Term, levels = rev(plot_pred_down_go$Term))
plot_pred_down_go$Genes <- sub('.', '', plot_pred_down_go$Genes)

# Creating data frame to be plotted
#UP
GO_plot_up <- data.frame(Gene = "",
                         Term = "",
                         Rho = "",
                         FDR = "",
                         Regulation = "")
cor_genes <- pred_up_cor_genes
rownum = 1

for (i in 1:nrow(plot_pred_up_go)){
  individual_genes <- plot_pred_up_go[i,10]
  individual_genes_vector <-  unlist(strsplit(individual_genes, ","))
  
  for (j in 1:length(individual_genes_vector)){
    Gene <- individual_genes_vector[j]
    Term <- as.character(plot_pred_up_go[i, 2])
    
    if(Gene %in% cor_genes$Gene_name) { #exclude genes that are not in cor meth data
      Rho <- cor_genes[cor_genes$Gene_name %in% Gene, 4]
      FDR <- cor_genes[cor_genes$Gene_name %in% Gene, 5]
      Regulation <- "upregulated"
      
      if(length(Rho) > 1){
        
        for(m in 1:length(Rho)){
          GO_plot_up[rownum, ] <- c(as.character(Gene), Term, Rho[m], FDR[m], Regulation)
          rownum <- rownum + 1
        }
        
      } else {
        GO_plot_up[rownum, ] <- c(as.character(Gene), Term, Rho, FDR, Regulation)
        rownum <- rownum + 1
        
      }
    }
    
    
  }
}



# DOWN
GO_plot_down <- data.frame(Gene = "",
                           Term = "",
                           Rho = "",
                           FDR = "",
                           Regulation = "")
cor_genes <- pred_down_cor_genes
rownum = 1

for (i in 1:nrow(plot_pred_down_go)){
  individual_genes <- plot_pred_down_go[i,10]
  individual_genes_vector <-  unlist(strsplit(individual_genes, ","))
  
  for (j in 1:length(individual_genes_vector)){
    Gene <- individual_genes_vector[j]
    Term <- as.character(plot_pred_down_go[i, 2])
    
    
    if(Gene %in% cor_genes$Gene_name) { #exclude genes that are not in cor meth data
      Rho <- cor_genes[cor_genes$Gene_name %in% Gene, 4]
      FDR <- cor_genes[cor_genes$Gene_name %in% Gene, 5]
      Regulation <- "downregulated"

    
    if(length(Rho) > 1){
      
      for(m in 1:length(Rho)){
        GO_plot_down[rownum, ] <- c(as.character(Gene), Term, Rho[m], FDR[m], Regulation)
        rownum <- rownum + 1
      }
      
    } else {
      GO_plot_down[rownum, ] <- c(as.character(Gene), Term, Rho, FDR, Regulation)
      rownum <- rownum + 1
      
    }
   }
  }
}

GO_plot <- rbind(GO_plot_up, GO_plot_down)
GO_plot$Rho <- as.numeric(GO_plot$Rho)
GO_plot$FDR <- as.numeric(GO_plot$FDR)
GO_plot$Term <- factor(GO_plot$Term, levels = rev(unique(GO_plot$Term))) # fixes order

write.table(GO_plot, "Output_files/Methylation/pred/GO/correlated_topvar/lfc0.5/GO_plot_pred_lfc0.5_topVar.txt", sep = "\t")

# # Plot 
# ggplot(GO_plot, aes(Rho,Term, size = FDR, fill = Regulation)) +
#   geom_point(pch = 21, colour = "black", stroke = 0.9, position = "jitter") +
#   scale_size_continuous(name="FDR", range = c(0.2,4), trans = 'reverse') + #reverse size of dots proportional to FDR)+
#   scale_fill_manual("", values = c("upregulated" = "#d91f24", "downregulated" = "#365196")) +
#   xlim(-1,1) +
#   geom_hline(yintercept = seq(1.5, 17.5, by = 1)) +
#   theme_minimal() +
#   ylab("") + xlab("Spearman's Rho") +
#   theme(legend.position = "bottom",
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         
#         axis.text.y = element_text(angle = 0, size = 15, vjust = 0.5),
#         legend.text = element_text(size = 15))


# Ridgeplot
pdf("Output_files/Methylation/pred/GO/correlated_topvar/lfc0.5/ridgeplot_GO_pred_lfc0.5_correlated_topVar_same_newTerms.pdf", width = 12, height = 10)
ridge_plot <- ggplot(GO_plot, aes(Rho,Term, fill = Regulation)) +
  scale_fill_manual("", values = c("upregulated" = "#d91f24", "downregulated" = "#365196")) +
  geom_density_ridges(alpha = .6) +
 # scale_size_continuous(name="FDR", range = c(0.2,4), trans = 'reverse') + #reverse size of dots proportional to FDR)+
  xlim(-1,1) +
  theme_ridges() +
  ylab("") + xlab("Spearman's Rho") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(angle = 0, size = 15, vjust = 0.5),
        legend.text = element_text(size = 15)+
          scale_y_discrete(limits = c("B cell receptor signaling pathway", #same neg
                                      "B cell proliferation", "B cell differentiation", "cell adhesion", "purine deoxyribonucleotide biosynthetic process", #new neg
                                      
                                      "cell surface receptor signaling pathway", #=
                                      
                                      "positive regulation of cytokine production", "regulation of macrophage activation", "inflammatory response",
                                      "positive regulation of NF-kappaB transcription factor activity",
                                      "neutrophil extravasation", "regulation of toll-like receptor signaling pathway", #same pos
                                      "positive regulation of nitric-oxide synthase biosynthetic process", "negative regulation of cytokine production",
                                      "interleukin-18-mediated signaling pathway", "production of molecular mediator involved in inflammatory response",
                                      "tumor necrosis factor production", "cytokine production involved in inflammatory response", "lipid homeostasis",
                                      "regulation of long-chain fatty acid import across plasma membrane"))) #new pos #new pos) 

ridge_plot
dev.off()

##################################################
# AntiTNF
rm(list = ls())

antiTNF_up_cor_genes <- read.table("Output_files/Methylation/antiTNF/antiTNF_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')
antiTNF_down_cor_genes <- read.table("Output_files/Methylation/antiTNF/antiTNF_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')

antiTNF_up_go <- read.csv("Output_files/Methylation/antiTNF/GO/antiTNF_correlated_up_DEGs_GO.txt", sep = '\t')
antiTNF_down_go <- read.table("Output_files/Methylation/antiTNF/GO/antiTNF_correlated_down_DEGs_GO.txt", sep = '\t')

# GO terms
terms <- c("immunoglobulin production", "adaptive immune response", "cartilage development", 
           "regulation of cell cycle G1/S phase transition",
           "interleukin-1-mediated signaling pathway", "positive regulation of chromosome condensation",
           
           "DNA replication-dependent chromatin assembly", "response to osmotic stress", 
           "positive regulation of signal transduction", "regulation of mRNA stability")


plot_antiTNF_up_go <- antiTNF_up_go %>% 
  mutate(plot = ifelse(Term %in% terms, "Yes", "No")) %>% 
  filter(plot == "Yes") %>% 
  mutate(direction = "Upregulated")
plot_antiTNF_up_go$Term <- Term(as.vector(plot_antiTNF_up_go$GO.ID))
plot_antiTNF_up_go$Term <- factor(plot_antiTNF_up_go$Term, levels = rev(plot_antiTNF_up_go$Term))
plot_antiTNF_up_go$Genes <- sub('.', '', plot_antiTNF_up_go$Genes)


plot_antiTNF_down_go <- antiTNF_down_go %>% 
  mutate(plot = ifelse(Term %in% terms, "Yes", "No")) %>% 
  filter(plot == "Yes") %>% 
  mutate(direction = "Downregulated")
plot_antiTNF_down_go$Term <- Term(as.vector(plot_antiTNF_down_go$GO.ID))
plot_antiTNF_down_go$Term <- factor(plot_antiTNF_down_go$Term, levels = rev(plot_antiTNF_down_go$Term))
plot_antiTNF_down_go$Genes <- sub('.', '', plot_antiTNF_down_go$Genes)

# Creating data frame to be plotted
#UP
GO_plot_up <- data.frame(Gene = "",
                         Term = "",
                         Rho = "",
                         FDR = "",
                         Regulation = "")
cor_genes <- antiTNF_up_cor_genes
rownum = 1

for (i in 1:nrow(plot_antiTNF_up_go)){
  individual_genes <- plot_antiTNF_up_go[i,10]
  individual_genes_vector <-  unlist(strsplit(individual_genes, ","))
  
  for (j in 1:length(individual_genes_vector)){
    Gene <- individual_genes_vector[j]
    Term <- as.character(plot_antiTNF_up_go[i, 2])
    Rho <- cor_genes[cor_genes$Gene_name %in% Gene, 4]
    FDR <- cor_genes[cor_genes$Gene_name %in% Gene, 5]
    Regulation <- "upregulated"
    
    if(length(Rho) > 1){
      
      for(m in 1:length(Rho)){
        GO_plot_up[rownum, ] <- c(as.character(Gene), Term, Rho[m], FDR[m], Regulation)
        rownum <- rownum + 1
      }
      
    } else {
      GO_plot_up[rownum, ] <- c(as.character(Gene), Term, Rho, FDR, Regulation)
      rownum <- rownum + 1
      
    }
    
  }
}



# DOWN
GO_plot_down <- data.frame(Gene = "",
                           Term = "",
                           Rho = "",
                           FDR = "",
                           Regulation = "")
cor_genes <- antiTNF_down_cor_genes
rownum = 1

for (i in 1:nrow(plot_antiTNF_down_go)){
  individual_genes <- plot_antiTNF_down_go[i,10]
  individual_genes_vector <-  unlist(strsplit(individual_genes, ","))
  
  for (j in 1:length(individual_genes_vector)){
    Gene <- individual_genes_vector[j]
    Term <- as.character(plot_antiTNF_down_go[i, 2])
    Rho <- cor_genes[cor_genes$Gene_name %in% Gene, 4]
    FDR <- cor_genes[cor_genes$Gene_name %in% Gene, 5]
    Regulation <- "downregulated"
    
    if(length(Rho) > 1){
      
      for(m in 1:length(Rho)){
        GO_plot_down[rownum, ] <- c(as.character(Gene), Term, Rho[m], FDR[m], Regulation)
        rownum <- rownum + 1
      }
      
    } else {
      GO_plot_down[rownum, ] <- c(as.character(Gene), Term, Rho, FDR, Regulation)
      rownum <- rownum + 1
      
    }
    
  }
}

GO_plot <- rbind(GO_plot_up, GO_plot_down)
GO_plot$Rho <- as.numeric(GO_plot$Rho)
GO_plot$FDR <- as.numeric(GO_plot$FDR)
GO_plot$Term <- factor(GO_plot$Term, levels = rev(unique(GO_plot$Term))) # fixes order

write.table(GO_plot, "Output_files/Methylation/GO_plot_antiTNF.txt", sep = "\t")
GO_plot <- read.csv( "Output_files/Methylation/GO_plot_antiTNF.txt", sep = "\t")

# Plot 
ggplot(GO_plot, aes(Rho,Term, size = FDR, fill = Regulation)) +
  geom_point(pch = 21, colour = "black", stroke = 0.9, position = "jitter") +
  scale_size_continuous(name="FDR", range = c(0.2,4), trans = 'reverse') + #reverse size of dots proportional to FDR)+
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


# Ridge plot
library(ggridges)
# excluding terms with too few sites: response to osmotic stress, interleukin-1-mediated signaling pathway and cartilage development
GO_plot_ridge <- GO_plot %>% 
  filter(!(Term == "response to osmotic stress") &
           !(Term == "interleukin-1-mediated signaling pathway") &
           !(Term == "cartilage development") &
           !(Term == "positive regulation of chromosome condensation"))


ggplot(GO_plot_ridge, aes(Rho,Term, fill = Regulation)) +
  scale_fill_manual("", values = c("upregulated" = "#d91f24", "downregulated" = "#365196")) +
  geom_density_ridges(alpha = .6) +
  scale_size_continuous(name="FDR", range = c(0.2,4), trans = 'reverse') + #reverse size of dots proportional to FDR)+
  xlim(-1,1) +
  #geom_hline(yintercept = seq(1.5, 12.5, by = 1)) +
  theme_ridges() +
  ylab("") + xlab("Spearman's Rho") +
  theme(legend.position = "bottom",
        #panel.grid.major.y = element_blank(),
        #panel.grid.minor.y = element_blank(),
        
        axis.text.y = element_text(angle = 0, size = 15, vjust = 0.5),
        legend.text = element_text(size = 15))



# geom_density_ridges(aes(point_size = FDR, point_color= Regulation, point_fill = Regulation),
#                     alpha = .6, jittered_points= T, point_alpha = 1, )+
#   scale_point_color_hue(l = 40) +