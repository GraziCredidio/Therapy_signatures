# EZE cohort 
# Expression plot - EZE cohort
  # Date last modification: 20.04

graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

library(tidyverse)

# Loading tables
res_sig <- read.table("Output_files\\DESeq2result_genenames_sig_padj_lfc.txt", sep = "\t")
coldata <- read.csv("Cleaned_tables\\EZECohort_ord.coldata_10.02.txt", sep = "\t")
counts <- read.csv("Cleaned_tables\\EZECohort_ord.counts_10.02.txt", sep = "\t")

# Data preprocessing ----
# NAs
coldata$smoking <- str_replace_na(coldata$smoking, "NA")
coldata$remission <- str_replace_na(coldata$remission, "NA")
coldata$bmi_class <- str_replace_na(coldata$bmi_class, "NA")

# As factor
cols_factor <- coldata %>% 
  select("diagnosis_class", "age_group", "sex", "remission", "bmi_class", 
         "smoking", "Azathioprin":"No_syst", "biologics" : "comorb_diabetes") %>% 
  names()

coldata <- coldata %>% 
  mutate_at(cols_factor, factor)

# Deseq object ----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = coldata,
                                     design = ~ diagnosis_class + sex + remission + age_group + 
                                       Azathioprin + MTX + Prednisolon + Leflunomid + Mycophenolat.Mofetil + No_syst + comorb_hypertension + comorb_diabetes +
                                       biologics)
dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized = TRUE)

# Expression plot ----
sig_norm_counts <- normalized_counts[res_sig$gene_ID, ] 
top_20 <- data.frame(sig_norm_counts)[1:20, ] %>%
  rownames_to_column(var = "ensgene")
top_20 <- gather(top_20, key = "sample_id", value = "normalized_counts", 2:648)
top_20 <- inner_join(top_20,
                     coldata,
                     by = "sample_id")

exp_plot <- ggplot(top_20) +
  geom_point(aes(x = ensgene, y = normalized_counts, color = biologics)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("Normalized Counts") +
  labs(title="Top 20 Significant DE Genes",
       subtitle = "padj < 0.05, LFC threshhold 0.1",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(exp_plot, file = "Output_files\\EZE_expressionPlot.png", height = 10, width = 10, units = "in", dpi = 300)


##### Expression plot - aza x no syst (R_crp) significant genes
rm(list = ls())

# Loading files ----
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

coldata_R_aza <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_aza_R_crp_04.05.txt", sep = "\t")
norm_counts_R_aza <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_normalized_aza_R_crp_04.05.txt",  sep = "\t")
sig_genes_R_aza <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_crp.txt", sep = "\t") 
topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t") #TVG all patients

sig_tvg <- intersect(rownames(topVarGenes), sig_genes_R_aza$feature)
sig_genes_tvgs <- sig_genes_R_aza[sig_genes_R_aza$feature %in% sig_tvg, ]
sig_genes_tvgs <- sig_genes_tvgs %>% 
  arrange(desc(abs(coef)))
sig_genes_tvgs <- sig_genes_tvgs[!(is.na(sig_genes_tvgs$genes) | sig_genes_tvgs$genes == ""), ] #11 pseudogenes

# Coldata prep ----
coldata_R_aza <- coldata_R_aza %>% 
  mutate(across(c(diagnosis_class, sex, Azathioprin, aza_vs_noSyst), as.factor))

# Counts prep ----
norm_counts_sig_tvgs <- norm_counts_R_aza[sig_genes_tvgs$feature, ]


top20_norm_counts_sig_tvgs <-  data.frame(norm_counts_sig_tvgs)[1:20, ] %>%
  rownames_to_column(var = "ensgene")

top20_norm_counts_sig_tvgs <- gather(top20_norm_counts_sig_tvgs, key = "sample_id", value = "norm_counts", 2:98)
top20_norm_counts_sig_tvgs <- inner_join(top20_norm_counts_sig_tvgs,
                         coldata_R_aza,
                         by = "sample_id")
top20_norm_counts_sig_tvgs$genes <- unique_ensg2gene[top20_norm_counts_sig_tvgs$ensgene, ]$hgnc_symbol


# Plots ----
exp_plot_vst_int <- ggplot(top20_norm_counts_sig_tvgs) +
  geom_point(aes(x = genes, y = norm_counts, color = aza_vs_noSyst)) +
  scale_y_log10(top20_norm_counts_sig_tvgs$norm_counts + 1) +
  xlab("Genes") +
  ylab("Normalized counts") +
  labs(title="Aza x No Syst - inactive disease",
       subtitle = "Highest aL2FC intersecting TVGs",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(exp_plot_vst_int, file = "Output_files/QC/expression_plot_significant_genes/expression_plot_top20_aza_R_crp_TVGs_aL2FC.png", height = 10, width = 10, units = "in", dpi = 300)




### Expression plot: anti TNF x no Bio (R crp) sig genes

rm(list = ls())

# Loading files ----
coldata_noBio <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_antiTNF_R_crp_04.05.txt", sep = "\t")
norm_counts_R_noBio <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_normalized_antiTNF_R_crp_04.05.txt",  sep = "\t")

sig_genes_noBio <- read.csv("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_R_crp.txt", sep = "\t")

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Coldata prep ----
coldata_noBio <- coldata_noBio %>% 
  mutate(across(c(diagnosis_class, sex, antiTNF_vs_noBiologics), as.factor))

# Counts prep ----
selected_genes <- sig_genes_noBio %>% 
  filter(abs(coef) > 0.5)

norm_counts_R_noBio <- norm_counts_R_noBio[selected_genes$feature, ]
norm_counts_R_noBio <- data.frame(norm_counts_R_noBio) %>% 
  rownames_to_column(var = "ensgene")

norm_counts_R_noBio_plot <- gather(norm_counts_R_noBio, key = "sample_id", value = "norm_counts", 2:138)

norm_counts_R_noBio_plot <- inner_join(norm_counts_R_noBio_plot,
                                       coldata_noBio,
                                       by = "sample_id")
norm_counts_R_noBio_plot$genes <- unique_ensg2gene[norm_counts_R_noBio_plot$ensgene, ]$hgnc_symbol


norm_counts_R_noBio_plot <- norm_counts_R_noBio_plot %>% 
           filter(!grepl('IGL', genes)) %>% 
           filter(!grepl('IGHV3', genes)) %>% 
           filter(!grepl('ENSG00000257275', ensgene))
  
# Plots ----
exp_plot_vst_int <- ggplot(norm_counts_R_noBio_plot) +
  geom_point(aes(x = genes, y = norm_counts, color = antiTNF_vs_noBiologics), alpha = .8) +
  scale_y_log10(norm_counts_R_noBio_plot$norm_counts + 1) +
  xlab("Genes") +
  ylab("Normalized count") +
  labs(title="Significant Genes",
       subtitle = "anti-TNF vs no Biologics - DEGs with abs(L2FC > 0.5)",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))
# saved in QC
