# EZE cohort: Therapy Signatures
  # Linear Mixed Model: analysis of the quantity of up and downregulated significant genes 
  # Figures 7A and 7B (Master's thesis) 
  # Author: Graziella Credidio

rm(list = ls())

# Loading files ----
library(tidyverse)

folder <- "Output_files/Maaslin2/plots"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading data----
pred_sig_genes <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_pred_vs_noSyst.txt", sep = "\t")
aza_sig_genes <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_aza_vs_noSyst.txt", sep = "\t")
antiTNF_sig_genes <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_antiTNF_vs_noBio.txt", sep = "\t")

# Creation of a single df with all significant genes
pred_sig_genes <- pred_sig_genes %>% 
  mutate(Coef_distribution = ifelse(coef > 0, "Upregulated", "Downregulated"))

aza_sig_genes <- aza_sig_genes %>% 
  mutate(Coef_distribution = ifelse(coef > 0, "Upregulated", "Downregulated"))

antiTNF_sig_genes <- antiTNF_sig_genes %>% 
  mutate(Coef_distribution = ifelse(coef > 0, "Upregulated", "Downregulated"))

all_sig_genes <- rbind(pred_sig_genes, aza_sig_genes, antiTNF_sig_genes)
all_sig_genes <- all_sig_genes %>% 
  mutate(metadata_plot = sub(".mod", "", metadata))
all_sig_genes$metadata_plot <- as.factor(all_sig_genes$metadata_plot)

# Creation of a single df with all significant genes with absLFC > 0.5
pred_sig_genes_lfc <- pred_sig_genes %>% 
  filter(abs(coef) > 0.5) 

aza_sig_genes_lfc <- aza_sig_genes %>% 
  filter(abs(coef) > 0.5) 

antiTNF_sig_genes_lfc <- antiTNF_sig_genes %>% 
  filter(abs(coef) > 0.5)

all_sig_genes_lfc <- rbind(pred_sig_genes_lfc, aza_sig_genes_lfc, antiTNF_sig_genes_lfc)
all_sig_genes_lfc <- all_sig_genes_lfc %>% 
  mutate(metadata_plot = sub(".mod", "", metadata))
all_sig_genes_lfc$metadata_plot <- as.factor(all_sig_genes_lfc$metadata_plot)

# Plots ----
# Figure 7A: all DEGs
all_degs_plot <- ggplot(all_sig_genes, aes(Coef_distribution, fill = Coef_distribution)) +
  geom_bar() +
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = 0.5, position = position_stack(vjust = 1.1), colour = "black", size = 15) +
  facet_grid(~metadata_plot) +
  scale_fill_manual("Direction", values = c("Downregulated" = "#365196", "Upregulated" = "#d91f24")) +
  theme_bw() +
  ylab("Number of DEGs") + xlab(" ") +
  theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
        text=element_text(size=40))
all_degs_plot
ggsave(all_degs_plot, file = 'Output_files/Maaslin2/plots/all_degs_coeff.png', height = 8, width = 10, units = "in", dpi = 300)

# Figure 7B: DEGs absLFC > 0.5
degs_lfc_plot <- ggplot(all_sig_genes_lfc, aes(Coef_distribution, fill = Coef_distribution)) +
  geom_bar() +
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = 0.5, position = position_stack(vjust = 1.1), colour = "black", size = 15) +
  facet_grid(~metadata_plot) +
  scale_fill_manual("Direction", values = c("Downregulated" = "#365196", "Upregulated" = "#d91f24")) +
  theme_bw() +
  ylab("Number of DEGs") + xlab(" ") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text=element_text(size=40))
degs_lfc_plot
ggsave(degs_lfc_plot, file = 'Output_files/Maaslin2/plots/degs_lfc_coeff.png', height = 8, width = 10, units = "in", dpi = 300)