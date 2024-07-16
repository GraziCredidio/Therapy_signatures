# EZE cohort: Therapy Signatures
  # DNA methylation-linked DEGs visualizations
  # Location of methylated sites and expression of DEGs
  # Figures 6B, 12B, 17B
  # Author: Graziella Credidio

rm(list = ls())

folders <- c("Output_files/Methylation/aza", "Output_files/Methylation/pred", "Output_files/Methylation/antiTNF")
for (i in folders){
  if (!dir.exists(i)) {
    dir.create(i)}
}


# Loading packages ----
library(tidyverse)

# Loading data ----
# Azathioprine
aza_up <- read.table("Output_files/Methylation/aza/aza_upregulated_DNAm_DEG.txt", sep = '\t')
aza_up$Regulation <- "Upregulated"

aza_down <- read.table("Output_files/Methylation/aza/aza_downregulated_DNAm_DEG.txt", sep = '\t')
aza_down$Regulation <- "Downregulated"

cor_genes_aza <- rbind(aza_up, aza_down)

# Prednisolone
pred_up <- read.table("Output_files/Methylation/pred/pred_upregulated_DNAm_DEG.txt", sep = '\t')
pred_up$Regulation <- "Upregulated"

pred_down <- read.table("Output_files/Methylation/pred/pred_downregulated_DNAm_DEG.txt", sep = '\t')
pred_down$Regulation <- "Downregulated"

cor_genes_pred <- rbind(pred_up, pred_down)

# Anti-TNF
antiTNF_up <- read.table("Output_files/Methylation/antiTNF/antiTNF_upregulated_DNAm_DEG.txt", sep = '\t')
antiTNF_up$Regulation <- "Upregulated"

antiTNF_down <- read.table("Output_files/Methylation/antiTNF/antiTNF_downregulated_DNAm_DEG.txt", sep = '\t')
antiTNF_down$Regulation <- "Downregulated"
cor_genes_antiTNF <- rbind(antiTNF_up, antiTNF_down)

# Data preprocessing: creating columns direction, position and analysis----
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


# Plots: position of methylated sites and direction of expression ----
meth_position_direction_plot <- function(data, subtitle, filePath) {
plot <-  ggplot(data, aes(Position, fill = Direction)) +
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

ggsave(plot, file = filePath, height = 8, width = 10, units = "in", dpi = 300)
  plot
}

meth_position_direction_plot(cor_genes_aza, "Azathioprine x No syst", "Output_files/Methylation/aza/aza_DNAm_DEG_direction_position.png")
meth_position_direction_plot(cor_genes_pred, "Prednisolone x No syst", "Output_files/Methylation/pred/pred_DNAm_DEG_direction_position.png")
meth_position_direction_plot(cor_genes_antiTNF, "Anti-TNF x No Biologics", "Output_files/Methylation/antiTNF/antiTNF_DNAm_DEG_direction_position.png")
