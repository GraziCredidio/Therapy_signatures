# EZE cohort: Therapy Signatures
  # Volcano plots of DEGs 
  # Figures 3B, 8B and 13B  (Master's thesis) 
  # Author: Graziella Credidio

rm(list = ls())


folder <- "Output_files/Volcano_plots"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading packages ----
library(data.table)
library(tidyverse)
library(ggrepel)

# Loading and preparing data ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

result_aza <- read.table("Output_files/Maaslin2/maaslin2_results_aza_vs_noSyst.txt", sep = "\t")
plot_data_aza <- data.frame()
result_aza <- subset(result_aza, result_aza$qval < 0.95 & abs(result_aza$coef) > 0.01)
result_aza$gene <- unique_ensg2gene[result_aza$feature, ]$hgnc_symbol
plot_data_aza <- rbind(plot_data_aza, result_aza)

result_pred <- read.table("Output_files/Maaslin2/maaslin2_results_pred_vs_noSyst.txt", sep = "\t")
plot_data_pred <- data.frame()
result_pred <- subset(result_pred, result_pred$qval < 0.95 & abs(result_pred$coef) > 0.01)
result_pred$gene <- unique_ensg2gene[result_pred$feature, ]$hgnc_symbol
plot_data_pred <- rbind(plot_data_pred, result_pred)


result_antiTNF <- read.table("Output_files/Maaslin2/maaslin2_results_antiTNF_vs_noBio.txt", sep = "\t")
plot_data_antiTNF <- data.frame()
result_antiTNF <- subset(result_antiTNF, result_antiTNF$qval < 0.95 & abs(result_antiTNF$coef) > 0.01)
result_antiTNF$gene <- unique_ensg2gene[result_antiTNF$feature, ]$hgnc_symbol
plot_data_antiTNF <- rbind(plot_data_antiTNF, result_antiTNF)

# Volcano plots ----
volcano <- function(data_plot, lfc_to_label, title, filePath){
  genes_label <- data_plot$gene[which(abs(data_plot$coef) > lfc_to_label)] 
  
  data_plot$color <- ifelse(data_plot$coef > 0.5 & data_plot$qval < 0.05, "red",
                                    ifelse(data_plot$coef < -0.5 & data_plot$qval < 0.05, "blue", 
                                           ifelse(data_plot$qval < 0.05, "black", "grey")))
  
  volcano_plot <- ggplot(data_plot, aes(x = coef, y = -log10(qval), text = paste("Gene:" = gene), color = color)) +
    geom_point(show.legend = FALSE) + 
    geom_label_repel(data = subset(data_plot, data_plot$gene %in% genes_label & (data_plot$qval) < 0.05),
                     aes(x = coef, y = -log10(qval), label = gene), 
                     nudge_x = -0.5, min.segment.length = 0,
                     #max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                     size = 9, alpha = 0.8) 
  volcano_plot <- volcano_plot + geom_vline(xintercept = 0, lty = 3) + 
    geom_vline(xintercept = 0.5, lty = 4) + 
    geom_vline(xintercept = -0.5, lty = 4) + 
    geom_hline(yintercept = 1.30103, lty = 4) + 
    geom_hline(yintercept = 3, lty = 4)
  volcano_plot <- volcano_plot + scale_color_manual(values = c("#606060", "#004C99",  "#C0C0C0", "#990000"))
  volcano_plot <- volcano_plot + theme_bw() + 
    theme(axis.text = element_text(size = 12), 
          axis.title = element_text(size = 15),
          plot.title = element_text(hjust = 0.5,size = 18),
          legend.position ="none") +
    ggtitle(title) + 
    xlab("Log2FoldChange") + ylab("-log10(qvalue)")
  ggsave(volcano_plot, file = filePath, height = 8, width = 10, units = "in", dpi = 300)
  volcano_plot
}

volcano(plot_data_aza, 1.8, "Azathioprine x No Syst", "Output_files/Volcano_plots/volcano_aza.png")
volcano(plot_data_pred, 1.3, "Pred x No Syst", "Output_files/Volcano_plots/volcano_pred.png")
volcano(plot_data_antiTNF, 0.9, "Anti-TNF x No Biologics", "Output_files/Volcano_plots/volcano_antiTNF.png")