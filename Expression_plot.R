# EZE cohort: Therapy Signatures
  # Anti-TNF: Expression plot of selected genes with absLFC > 0.5
  # Figure 3 C 
  # Author: Graziella Credidio

rm(list = ls())

folder <- "Output_files/Expression_plot"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading packages ----
library(tidyverse)

# Loading files ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_antiTNF <- read.table("Cleaned_tables/models/antiTNF/EZECohort_coldata_antiTNF.txt", sep = "\t")
norm_counts_antiTNF <- read.table("Cleaned_tables/models/antiTNF/normalized_counts_antiTNF.txt", sep = "\t")
degs_antiTNF <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_antiTNF_vs_noBio.txt", sep = "\t")

# Data preprocessing ----
coldata_antiTNF <- coldata_antiTNF %>% 
  mutate(across(c(diagnosis_class, sex, antiTNF_vs_noBiologics), as.factor))

selected_degs_antiTNF <- degs_antiTNF %>% 
  filter(abs(coef) > 0.5)

norm_counts_antiTNF <- norm_counts_antiTNF[selected_degs_antiTNF$feature, ]
norm_counts_antiTNF <- data.frame(norm_counts_antiTNF) %>% 
  rownames_to_column(var = "ensgene")
norm_counts_antiTNF_plot <- gather(norm_counts_antiTNF, key = "sample_id", value = "norm_counts", 2:138)
norm_counts_antiTNF_plot <- inner_join(norm_counts_antiTNF_plot,
                                       coldata_antiTNF,
                                       by = "sample_id")
norm_counts_antiTNF_plot$genes <- unique_ensg2gene[norm_counts_antiTNF_plot$ensgene, ]$hgnc_symbol

norm_counts_antiTNF_plot <- norm_counts_antiTNF_plot %>% 
           filter(!grepl('IGL', genes)) %>% 
           filter(!grepl('IGHV3', genes)) %>% 
           filter(!grepl('ENSG00000257275', ensgene))
  
# Plots ----
exp_plot <- ggplot(norm_counts_antiTNF_plot) +
  geom_point(aes(x = genes, y = norm_counts, color = antiTNF_vs_noBiologics), alpha = .8) +
  scale_y_log10(norm_counts_antiTNF_plot$norm_counts + 1) +
  xlab("Genes") +
  ylab("log10(normalized counts)") +
  labs(title="Expression plot - anti-TNF vs no Biologics",
       subtitle = "DEGs with abs(L2FC > 0.5)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
exp_plot
ggsave(exp_plot,
       file = "Output_files/Expression_plot/exp_plot_antiTNF_degs_lfc_0.5.png",
       height = 10, width = 6, units = "in", dpi = 300)
