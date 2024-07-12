# EZE cohort 
# Gene ontology of significant genes from model (maaslin2)
# Azathioprin x No Systemic therapies
  # Date: 05.04


graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC
# Loading packages ----
if (!requireNamespace("topGO", quietly = TRUE)) BiocManager::install("topGO")
if (!requireNamespace("genefilter", quietly = TRUE)) BiocManager::install("genefilter")
if (!requireNamespace("geneplotter", quietly = TRUE)) BiocManager::install("geneplotter")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rgraphviz")

install.packages("data.table")

library(genefilter)
library(geneplotter)
library(topGO)
library(data.table)
library(tidyverse)
library(ggrepel)
library(ggthemes)

folder <- "Output_files/GO/All_patients/aza_noSyst"

if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading data ----
sig_genes <- read.csv("Output_files/Maaslin2/significant_results/maaslin2_significant_results_aza_vs_noSyst_vst_q0.05.txt", sep = "\t")

coldata <- read.csv("Output_files/Maaslin2/Tables/coldata_maaslin2_aza_vs_noSyst_16.03.txt", sep = "\t")
counts <- read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_28.02.txt", sep="\t")

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Data prep ----
counts <- counts[,colnames(counts) %in% coldata$sample_id]
idx <- match(colnames(counts), rownames(coldata)) 
ord.coldata <- coldata[idx, ]

dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = ord.coldata,
                                     design = ~ 1)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)

dds <- DESeq(dds_counts, betaPrior = FALSE)

res_deseq <- results(dds, alpha = 0.05, independentFiltering = TRUE)
res_deseq <- as.data.frame(res_deseq) %>%
  arrange(padj) 

res_deseq$gene <- unique_ensg2gene[rownames(res_deseq), ]$hgnc_symbol
res_deseq$gene_id <- rownames(res_deseq)

write.table(res_deseq,"Output_files/DESeq2/TopGO/DESeq2_baseMean_aza_noSyst_31.03.txt", sep="\t", row.names=TRUE)

# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/DESeq2_baseMean_aza_noSyst_31.03.txt", sep="\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) # Filtering out genes that were never detected 

sig_genes_id <- subset(sig_genes, sig_genes$coef > 0) #modify coef
sig_genes_id <- sig_genes_id$feature
all_genes_id <- rownames(res_deseq)


## Removing background
overallBaseMean <- as.matrix(res_deseq[, "baseMean", drop = FALSE])
sig_idx <- match(sig_genes_id, rownames(overallBaseMean))

backG <- c()
for (i in sig_idx) {
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
}

backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]
backG <- setdiff(backG, sig_genes_id)

multidensity(list(all = log2(res_deseq[, "baseMean"]), 
                  foreground = log2(res_deseq[sig_genes_id, "baseMean"]), 
                  background = log2(res_deseq[backG, "baseMean"])), 
             xlab = "log2 mean normalized counts",
             main = "Matching for enrichment analysis") 

## Selecting genes 
inUniverse <- all_genes_id %in% c(sig_genes_id, backG)
inSelection <- all_genes_id %in% sig_genes_id
alg <- factor(as.integer(inSelection[inUniverse]))
names(alg) <- all_genes_id[inUniverse]

## Preparing data
tgd <- new("topGOdata", ontology = "BP", allGenes = alg, nodeSize = 5, annot = annFUN.org,
           mapping = "org.Hs.eg.db", ID = "ensembl") # Change to MF and CC


## Running tests
resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher")

if (length(nodes(graph(tgd))) < 200) {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim", numChar=1000, topNodes = length(nodes(graph(tgd))))
} else {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim",numChar=1000, topNodes = 200)
}

printGraph(tgd, resultTopGO.elim, firstSigNodes = 5, 
           fn.prefix = "Output_files/GO/All_patients/aza_noSyst/aza_noSyst_elim_NoBackground_BP_pos",
           useInfo = "all", pdfSW = TRUE) 


## Finding the genes associated with the enriched GO terms
for (i in 1:nrow(tab)) {
  go_id <- as.vector(tab[i, 1])
  genes_in_cat <- intersect(unlist(genesInTerm(tgd, go_id)), sig_genes_id)
  gene_sym_in_cat <- unique_ensg2gene[as.character(genes_in_cat), "hgnc_symbol"]
  gene_sym_in_cat_str <- ""
  
  if (length(genes_in_cat) > 0) {
    for (j in 1:length(gene_sym_in_cat)) {
      gene_sym_in_cat_str <- paste(gene_sym_in_cat_str, gene_sym_in_cat[j],
                                   sep = ",")
    }
  }
  if (gene_sym_in_cat_str == "") {
    gene_sym_in_cat_str <- NA
  }
  
  tab$Genes[i] <- gene_sym_in_cat_str
  
}

tab_aza <- copy(tab)

## Writing tables
### all terms
tab_aza <- tab_aza %>% 
  mutate(geneRatio = tab_aza[,4]/tab_aza[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_aza,"Output_files/GO/All_patients/aza_noSyst/TopGO_aza_noSyst_BP_pos.txt", sep="\t",row.names=TRUE)

### significant terms
tab_aza_sig <- tab_aza %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_aza_sig$Term = with(tab_aza_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_aza_sig,"Output_files/GO/All_patients/aza_noSyst/TopGO_significantTerms_aza_noSyst_BP_pos.txt", sep="\t",row.names=TRUE)

## Barplot
tab_aza_sig_plot <- tab_aza_sig[1:30,]

gg <- ggplot(tab_aza_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = reorder(Term, -Fisher.elim)), position = "stack", 
                    stat = "identity", fill="darkorange", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("GO Term") + ggtitle("GO Azathioprin x no Systemic Therapies (coeff > 0)")
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/All_patients/aza_noSyst/TopGO_aza_noSyst_BP_pos.png", height = 8, width = 10, units = "in", dpi = 300)






###################################################### Remitters ######################################################################
graphics.off()
rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission/aza_noSyst"

if (!dir.exists(folder)) {
  dir.create(folder)
}



# Loading data ----
sig_genes <- read.csv("Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_vst_q0.05.txt", sep = "\t")

coldata <- read.csv("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_aza_vs_noSyst_R_17.03.txt", sep = "\t")
counts <- read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_28.02.txt", sep="\t")

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Data prep ----
counts <- counts[,colnames(counts) %in% coldata$sample_id]
idx <- match(colnames(counts), rownames(coldata)) 
ord.coldata <- coldata[idx, ]

dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = ord.coldata,
                                     design = ~ 1)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)

dds <- DESeq(dds_counts, betaPrior = FALSE)

res_deseq <- results(dds, alpha = 0.05, independentFiltering = TRUE)

res_deseq$gene <- unique_ensg2gene[rownames(res_deseq), ]$hgnc_symbol
res_deseq$gene_id <- rownames(res_deseq)

write.table(res_deseq,"Output_files/DESeq2/TopGO/Remission/DESeq2_baseMean_aza_noSyst.txt", sep="\t", row.names=TRUE)

# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission/DESeq2_baseMean_aza_noSyst.txt", sep="\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) # Filtering out genes that were never detected

sig_genes_id <- subset(sig_genes, sig_genes$coef < 0) #change
sig_genes_id <- sig_genes_id$feature
all_genes_id <- rownames(res_deseq)

## Removing background
overallBaseMean <- as.matrix(res_deseq[, "baseMean", drop = FALSE])
sig_idx <- match(sig_genes_id, rownames(overallBaseMean))

backG <- c()
for (i in sig_idx) {
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
}
backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]
backG <- setdiff(backG, sig_genes_id)

multidensity(list(all = log2(res_deseq[, "baseMean"]), 
                  foreground = log2(res_deseq[sig_genes_id, "baseMean"]), 
                  background = log2(res_deseq[backG, "baseMean"])), 
             xlab = "log2 mean normalized counts",
             main = "Matching for enrichment analysis") 

## Selecting genes 
inUniverse <- all_genes_id %in% c(sig_genes_id, backG)
inSelection <- all_genes_id %in% sig_genes_id
alg <- factor(as.integer(inSelection[inUniverse]))
names(alg) <- all_genes_id[inUniverse]

## Preparing data
tgd <- new("topGOdata", ontology = "BP", allGenes = alg, nodeSize = 5, annot = annFUN.org,
           mapping = "org.Hs.eg.db", ID = "ensembl") # Change to MF and CC


## Running tests
resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher")

if (length(nodes(graph(tgd))) < 200) {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim", numChar=1000, topNodes = length(nodes(graph(tgd))))
} else {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim",numChar=1000, topNodes = 200)
}

printGraph(tgd, resultTopGO.elim, firstSigNodes = 5, 
           fn.prefix = "Output_files/GO/Remission/aza_noSyst/aza_noSyst_elim_NoBackground_BP_neg",
           useInfo = "all", pdfSW = TRUE) 


## Finding the genes associated with the enriched GO terms
for (i in 1:nrow(tab)) {
  go_id <- as.vector(tab[i, 1])
  genes_in_cat <- intersect(unlist(genesInTerm(tgd, go_id)), sig_genes_id)
  gene_sym_in_cat <- unique_ensg2gene[as.character(genes_in_cat), "hgnc_symbol"]
  gene_sym_in_cat_str <- ""
  
  if (length(genes_in_cat) > 0) {
    for (j in 1:length(gene_sym_in_cat)) {
      gene_sym_in_cat_str <- paste(gene_sym_in_cat_str, gene_sym_in_cat[j],
                                   sep = ",")
    }
  }
  if (gene_sym_in_cat_str == "") {
    gene_sym_in_cat_str <- NA
  }
  
  tab$Genes[i] <- gene_sym_in_cat_str
  
}

tab_aza <- copy(tab)

## Writing tables
### all terms
tab_aza <- tab_aza %>% 
  mutate(geneRatio = tab_aza[,4]/tab_aza[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_aza,"Output_files/GO/Remission/aza_noSyst/TopGO_aza_noSyst_BP_neg.txt", sep="\t",row.names=TRUE)

### significant terms
tab_aza_sig <- tab_aza %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_aza_sig$Term = with(tab_aza_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_aza_sig,"Output_files/GO/Remission/aza_noSyst/TopGO_significantTerms_aza_noSyst_BP_neg.txt", sep="\t",row.names=TRUE)

## Barplot
tab_aza_sig_plot <- tab_aza_sig[1:30,]

gg <- ggplot(tab_aza_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = reorder(Term, -Fisher.elim)), position = "stack", 
                    stat = "identity", fill="skyblue4", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=12, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("GO Term") + ggtitle("GO Azathioprin x no Systemic Therapies - remitters (coeff < 0)")
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission/aza_noSyst/TopGO_aza_noSyst_BP_neg.png", height = 8, width = 10, units = "in", dpi = 300)



##################################################### Intersection genes (topVar + Remitters)
graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Intersections/aza_noSyst"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
list_sig_genes_all_remitters <- readRDS("Output_files/Maaslin2/significant_results/list_sig_genes_all_remitters_intersect_diff_06.04.rds")
int_R_aza_topvar <- list_sig_genes_all_remitters$intersect_R_aza_topVar

sig_genes <- read.csv("Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_vst_q0.05.txt", sep = "\t")

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Subsetting sig genes to have only intersecting genes 
int_sig_genes <- sig_genes[sig_genes$feature %in% int_R_aza_topvar,]


# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission/DESeq2_baseMean_aza_noSyst.txt", sep="\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) # Filtering out genes that were never detected

sig_genes_id <- subset(int_sig_genes, int_sig_genes$coef > 0) #modify coef
sig_genes_id <- sig_genes_id$feature
all_genes_id <- rownames(res_deseq)

## Removing background
overallBaseMean <- as.matrix(res_deseq[, "baseMean", drop = FALSE])
sig_idx <- match(sig_genes_id, rownames(overallBaseMean))

backG <- c()
for (i in sig_idx) {
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
}
backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]
backG <- setdiff(backG, sig_genes_id)

multidensity(list(all = log2(res_deseq[, "baseMean"]), 
                  foreground = log2(res_deseq[sig_genes_id, "baseMean"]), 
                  background = log2(res_deseq[backG, "baseMean"])), 
             xlab = "log2 mean normalized counts",
             main = "Matching for enrichment analysis") 

## Selecting genes 
inUniverse <- all_genes_id %in% c(sig_genes_id, backG)
inSelection <- all_genes_id %in% sig_genes_id
alg <- factor(as.integer(inSelection[inUniverse]))
names(alg) <- all_genes_id[inUniverse]

## Preparing data
tgd <- new("topGOdata", ontology = "BP", allGenes = alg, nodeSize = 5, annot = annFUN.org,
           mapping = "org.Hs.eg.db", ID = "ensembl") # Change to MF and CC


## Running tests
resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher")

if (length(nodes(graph(tgd))) < 200) {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim", numChar=1000, topNodes = length(nodes(graph(tgd))))
} else {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim",numChar=1000, topNodes = 200)
}

printGraph(tgd, resultTopGO.elim, firstSigNodes = 5, 
           fn.prefix = "Output_files/GO/Intersections/aza_noSyst/R_aza_noSyst_elim_NoBackground_BP_pos",
           useInfo = "all", pdfSW = TRUE) 


## Finding the genes associated with the enriched GO terms
for (i in 1:nrow(tab)) {
  go_id <- as.vector(tab[i, 1])
  genes_in_cat <- intersect(unlist(genesInTerm(tgd, go_id)), sig_genes_id)
  gene_sym_in_cat <- unique_ensg2gene[as.character(genes_in_cat), "hgnc_symbol"]
  gene_sym_in_cat_str <- ""
  
  if (length(genes_in_cat) > 0) {
    for (j in 1:length(gene_sym_in_cat)) {
      gene_sym_in_cat_str <- paste(gene_sym_in_cat_str, gene_sym_in_cat[j],
                                   sep = ",")
    }
  }
  if (gene_sym_in_cat_str == "") {
    gene_sym_in_cat_str <- NA
  }
  
  tab$Genes[i] <- gene_sym_in_cat_str
  
}

tab_aza <- copy(tab)

## Writing tables
### all terms
tab_aza <- tab_aza %>% 
  mutate(geneRatio = tab_aza[,4]/tab_aza[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_aza,"Output_files/GO/Intersections/aza_noSyst/TopGO_R_aza_noSyst_BP_pos.txt", sep="\t",row.names=TRUE)

### significant terms
tab_aza_sig <- tab_aza %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_aza_sig$Term = with(tab_aza_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_aza_sig,"Output_files/GO/Intersections/aza_noSyst/TopGO_R_significantTerms_aza_noSyst_BP_pos.txt", sep="\t",row.names=TRUE)

## Barplot
tab_aza_sig_plot <- tab_aza_sig[1:30,]
gg <- ggplot(tab_aza_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="darkorange", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("") +
  labs(
    title = "GO Azathioprin x no Systemic Therapies (coef > 0)",
    subtitle = "Intersection: Inactive + TopVar genes"
  ) +
  theme(
    plot.title = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 12, vjust = 1, hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = -15))
  )
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Intersections/aza_noSyst/TopGO_R_aza_noSyst_BP_pos.png", height = 8, width = 10, units = "in", dpi = 300)






#################################################  Intersection var part > 25% and top var genes (remitters) ###################################################
graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Intersections/aza_25varpart"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
top25percent_varPart <- read.table("Output_files/VariancePartition/Remission/sig_genes/aza/lab_values/genes_top25percent_varPart_intersectVenn_aza_R_lab.txt")
top25percent_varPart_genes <- top25percent_varPart$Genes

sig_genes <- read.csv("Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_vst_q0.05.txt", sep = "\t")

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Subsetting sig genes to have only intersecting genes 
top25percent_varPart_sig <- sig_genes[sig_genes$feature %in% top25percent_varPart_genes,]

# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission/DESeq2_baseMean_aza_noSyst.txt", sep="\t")

sig_genes_id <- subset(top25percent_varPart_sig, top25percent_varPart_sig$coef < 0) #modify coef
sig_genes_id <- sig_genes_id$feature
all_genes_id <- rownames(res_deseq)

## Removing background
overallBaseMean <- as.matrix(res_deseq[, "baseMean", drop = FALSE])
sig_idx <- match(sig_genes_id, rownames(overallBaseMean))

backG <- c()
for (i in sig_idx) {
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
}
backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]
backG <- setdiff(backG, sig_genes_id)

multidensity(list(all = log2(res_deseq[, "baseMean"]), 
                  foreground = log2(res_deseq[sig_genes_id, "baseMean"]), 
                  background = log2(res_deseq[backG, "baseMean"])), 
             xlab = "log2 mean normalized counts",
             main = "Matching for enrichment analysis") 

## Selecting genes 
inUniverse <- all_genes_id %in% c(sig_genes_id, backG)
inSelection <- all_genes_id %in% sig_genes_id
alg <- factor(as.integer(inSelection[inUniverse]))
names(alg) <- all_genes_id[inUniverse]

## Preparing data
tgd <- new("topGOdata", ontology = "BP", allGenes = alg, nodeSize = 5, annot = annFUN.org,
           mapping = "org.Hs.eg.db", ID = "ensembl") # Change to MF and CC


## Running tests
resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher")

if (length(nodes(graph(tgd))) < 200) {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim", numChar=1000, topNodes = length(nodes(graph(tgd))))
} else {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim",numChar=1000, topNodes = 200)
}

printGraph(tgd, resultTopGO.elim, firstSigNodes = 5, 
           fn.prefix = "Output_files/GO/Intersections/aza_25varPart/R_aza_noSyst_elim_NoBackground_BP_neg",
           useInfo = "all", pdfSW = TRUE) 


## Finding the genes associated with the enriched GO terms
for (i in 1:nrow(tab)) {
  go_id <- as.vector(tab[i, 1])
  genes_in_cat <- intersect(unlist(genesInTerm(tgd, go_id)), sig_genes_id)
  gene_sym_in_cat <- unique_ensg2gene[as.character(genes_in_cat), "hgnc_symbol"]
  gene_sym_in_cat_str <- ""
  
  if (length(genes_in_cat) > 0) {
    for (j in 1:length(gene_sym_in_cat)) {
      gene_sym_in_cat_str <- paste(gene_sym_in_cat_str, gene_sym_in_cat[j],
                                   sep = ",")
    }
  }
  if (gene_sym_in_cat_str == "") {
    gene_sym_in_cat_str <- NA
  }
  
  tab$Genes[i] <- gene_sym_in_cat_str
  
}

tab_aza <- copy(tab)

## Writing tables
### all terms
tab_aza <- tab_aza %>% 
  mutate(geneRatio = tab_aza[,4]/tab_aza[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_aza,"Output_files/GO/Intersections/aza_25varpart/TopGO_R_aza_varPart_BP_neg.txt", sep="\t",row.names=TRUE)

### significant terms
tab_aza_sig <- tab_aza %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_aza_sig$Term = with(tab_aza_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_aza_sig,"Output_files/GO/Intersections/aza_25varpart/TopGO_R_significantTerms_aza_varPart_BP_neg.txt", sep="\t",row.names=TRUE)

## Barplot
tab_aza_sig_plot <- tab_aza_sig[1:30,]
gg <- ggplot(tab_aza_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="skyblue4", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("") +
  labs(
    title = "GO Azathioprin x no Systemic Therapies (coef > 0)",
    subtitle = "Intersection: VarPart > 25% + top var"
  ) +
  theme(
    plot.title = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 12, vjust = 1, hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = -15))
  )
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Intersections/aza_25varpart/TopGO_R_aza_varPart_BP_neg.png", height = 8, width = 10, units = "in", dpi = 300)







###################################################### Remitters + crp ######################################################################
graphics.off()
rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission_crp/aza_noSyst"

if (!dir.exists(folder)) {
  dir.create(folder)
}



# Loading data ----
sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_crp.txt", sep = "\t")
coldata <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_aza_R_crp_04.05.txt", sep = "\t")
counts <- read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_28.02.txt", sep="\t")

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Data prep ----
counts <- counts[,colnames(counts) %in% coldata$sample_id]
idx <- match(colnames(counts), rownames(coldata)) 
ord.coldata <- coldata[idx, ]

dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = ord.coldata,
                                     design = ~ 1)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)

dds <- DESeq(dds_counts, betaPrior = FALSE)

res_deseq <- results(dds, alpha = 0.05, independentFiltering = TRUE)

res_deseq$gene <- unique_ensg2gene[rownames(res_deseq), ]$hgnc_symbol
res_deseq$gene_id <- rownames(res_deseq)

write.table(res_deseq,"Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_aza_R_crp.txt", sep="\t", row.names=TRUE)

# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_aza_R_crp.txt", sep="\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) # Filtering out genes that were never detected

sig_genes_id <- subset(sig_genes, sig_genes$coef > 0) #change
sig_genes_id <- sig_genes_id$feature
all_genes_id <- rownames(res_deseq)

## Removing background
overallBaseMean <- as.matrix(res_deseq[, "baseMean", drop = FALSE])
sig_idx <- match(sig_genes_id, rownames(overallBaseMean))

backG <- c()
for (i in sig_idx) {
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
}
backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]
backG <- setdiff(backG, sig_genes_id)

multidensity(list(all = log2(res_deseq[, "baseMean"]), 
                  foreground = log2(res_deseq[sig_genes_id, "baseMean"]), 
                  background = log2(res_deseq[backG, "baseMean"])), 
             xlab = "log2 mean normalized counts",
             main = "Matching for enrichment analysis") 

## Selecting genes 
inUniverse <- all_genes_id %in% c(sig_genes_id, backG)
inSelection <- all_genes_id %in% sig_genes_id
alg <- factor(as.integer(inSelection[inUniverse]))
names(alg) <- all_genes_id[inUniverse]

## Preparing data
tgd <- new("topGOdata", ontology = "BP", allGenes = alg, nodeSize = 5, annot = annFUN.org,
           mapping = "org.Hs.eg.db", ID = "ensembl") # Change to MF and CC


## Running tests
resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher")

if (length(nodes(graph(tgd))) < 200) {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim", numChar=1000, topNodes = length(nodes(graph(tgd))))
} else {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim",numChar=1000, topNodes = 200)
}

printGraph(tgd, resultTopGO.elim, firstSigNodes = 5, 
           fn.prefix = "Output_files/GO/Remission_crp/aza_noSyst/aza_noSyst_elim_NoBackground_BP_pos",
           useInfo = "all", pdfSW = TRUE) 


## Finding the genes associated with the enriched GO terms
for (i in 1:nrow(tab)) {
  go_id <- as.vector(tab[i, 1])
  genes_in_cat <- intersect(unlist(genesInTerm(tgd, go_id)), sig_genes_id)
  gene_sym_in_cat <- unique_ensg2gene[as.character(genes_in_cat), "hgnc_symbol"]
  gene_sym_in_cat_str <- ""
  
  if (length(genes_in_cat) > 0) {
    for (j in 1:length(gene_sym_in_cat)) {
      gene_sym_in_cat_str <- paste(gene_sym_in_cat_str, gene_sym_in_cat[j],
                                   sep = ",")
    }
  }
  if (gene_sym_in_cat_str == "") {
    gene_sym_in_cat_str <- NA
  }
  
  tab$Genes[i] <- gene_sym_in_cat_str
  
}

tab_aza <- copy(tab)

## Writing tables
### all terms
tab_aza <- tab_aza %>% 
  mutate(geneRatio = tab_aza[,4]/tab_aza[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_aza,"Output_files/GO/Remission_crp/aza_noSyst/TopGO_BP_aza_R_crp_pos.txt", sep="\t",row.names=TRUE)

### significant terms
tab_aza_sig <- tab_aza %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_aza_sig$Term = with(tab_aza_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_aza_sig,"Output_files/GO/Remission_crp/aza_noSyst/TopGO_BP_significantTerms_aza_R_crp_pos.txt", sep="\t",row.names=TRUE)

## Barplot
tab_aza_sig_plot <- tab_aza_sig[1:30,]

gg <- ggplot(tab_aza_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = reorder(Term, -Fisher.elim)), position = "stack", 
                    stat = "identity", fill="darkorange", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=12, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("GO Term") + ggtitle("GO BP Aza x No systemic therapies - remitters + crp (upregulated)")
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission_crp/aza_noSyst/TopGO_BP_aza_R_crp_pos.png", height = 8, width = 10, units = "in", dpi = 300)




################### Intersection genes (topVar and Remitters + crp) #########################
graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission_crp/aza_noSyst/aza_topVar"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_crp.txt", sep = "\t")
topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")

intersect_sig_topvar <- base::intersect(sig_genes$feature, rownames(topVar))

sig_topVar <- sig_genes[sig_genes$feature %in% intersect_sig_topvar,]

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_aza_R_crp.txt", sep="\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) # Filtering out genes that were never detected

sig_genes_id <- subset(sig_topVar, sig_topVar$coef < 0) #modify coef
sig_genes_id <- sig_genes_id$feature
all_genes_id <- rownames(res_deseq)

## Removing background
overallBaseMean <- as.matrix(res_deseq[, "baseMean", drop = FALSE])
sig_idx <- match(sig_genes_id, rownames(overallBaseMean))

backG <- c()
for (i in sig_idx) {
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
}
backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]
backG <- setdiff(backG, sig_genes_id)

multidensity(list(all = log2(res_deseq[, "baseMean"]), 
                  foreground = log2(res_deseq[sig_genes_id, "baseMean"]), 
                  background = log2(res_deseq[backG, "baseMean"])), 
             xlab = "log2 mean normalized counts",
             main = "Matching for enrichment analysis") 

## Selecting genes 
inUniverse <- all_genes_id %in% c(sig_genes_id, backG)
inSelection <- all_genes_id %in% sig_genes_id
alg <- factor(as.integer(inSelection[inUniverse]))
names(alg) <- all_genes_id[inUniverse]

## Preparing data
tgd <- new("topGOdata", ontology = "BP", allGenes = alg, nodeSize = 5, annot = annFUN.org,
           mapping = "org.Hs.eg.db", ID = "ensembl") # Change to MF and CC


## Running tests
resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher")

if (length(nodes(graph(tgd))) < 200) {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim", numChar=1000, topNodes = length(nodes(graph(tgd))))
} else {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim",numChar=1000, topNodes = 200)
}

printGraph(tgd, resultTopGO.elim, firstSigNodes = 5, 
           fn.prefix = "Output_files/GO/Remission_crp/aza_noSyst/aza_topvar/aza_elim_NoBackground_BP_neg",
           useInfo = "all", pdfSW = TRUE) 


## Finding the genes associated with the enriched GO terms
for (i in 1:nrow(tab)) {
  go_id <- as.vector(tab[i, 1])
  genes_in_cat <- intersect(unlist(genesInTerm(tgd, go_id)), sig_genes_id)
  gene_sym_in_cat <- unique_ensg2gene[as.character(genes_in_cat), "hgnc_symbol"]
  gene_sym_in_cat_str <- ""
  
  if (length(genes_in_cat) > 0) {
    for (j in 1:length(gene_sym_in_cat)) {
      gene_sym_in_cat_str <- paste(gene_sym_in_cat_str, gene_sym_in_cat[j],
                                   sep = ",")
    }
  }
  if (gene_sym_in_cat_str == "") {
    gene_sym_in_cat_str <- NA
  }
  
  tab$Genes[i] <- gene_sym_in_cat_str
  
}

tab_aza <- copy(tab)

## Writing tables
### all terms
tab_aza <- tab_aza %>% 
  mutate(geneRatio = tab_aza[,4]/tab_aza[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_aza,"Output_files/GO/Remission_crp/aza_noSyst/aza_topvar/TopGO_BP_intersect_aza_R_crp_topVar_neg.txt", sep="\t",row.names=TRUE)

### significant terms
tab_aza_sig <- tab_aza %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_aza_sig$Term = with(tab_aza_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_aza_sig,"Output_files/GO/Remission_crp/aza_noSyst/aza_topvar/TopGO_BP_sigTerms_intersect_aza_R_crp_topVar_neg.txt", sep="\t",row.names=TRUE)

## Barplot
tab_aza_sig_plot <- tab_aza_sig[1:30,]
gg <- ggplot(tab_aza_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="skyblue4", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("") +
  labs(
    title = "GO Azathioprin x no Systemic Therapies (downregulated)",
    subtitle = "Intersection: Inactive (crp) + TopVar genes"
  ) +
  theme(
    plot.title = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 12, vjust = 1, hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = -15))
  )
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission_crp/aza_noSyst/aza_topvar/TopGO_BP_intersect_aza_R_crp_topVar_neg.png", height = 8, width = 10, units = "in", dpi = 300)





###################################################### Remitters + crp: var part 25% + top var ######################################################################
graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission_crp/aza_25varpart_topvar"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
top25percent_varPart <- read.table("Output_files/VariancePartition/Remission_crp/aza/top25percent_varPart_aza_R_crp.txt")
top25percent_varPart_genes <- rownames(top25percent_varPart)

topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVar_genes <- rownames(topVar) 
intersect_topVar_varPart25 <- base::intersect(top25percent_varPart_genes, topVar_genes)

sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_aza_vs_noSyst_R_crp.txt", sep = "\t")

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Subsetting sig genes to have only intersecting genes 
intersect_topVar_varPart25_sig <- sig_genes[sig_genes$feature %in% intersect_topVar_varPart25,]

# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_aza_R_crp.txt", sep="\t")

sig_genes_id <- subset(intersect_topVar_varPart25_sig, intersect_topVar_varPart25_sig$coef > 0) #modify coef
sig_genes_id <- sig_genes_id$feature
all_genes_id <- rownames(res_deseq)

## Removing background
overallBaseMean <- as.matrix(res_deseq[, "baseMean", drop = FALSE])
sig_idx <- match(sig_genes_id, rownames(overallBaseMean))

backG <- c()
for (i in sig_idx) {
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
}
backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]
backG <- setdiff(backG, sig_genes_id)

multidensity(list(all = log2(res_deseq[, "baseMean"]), 
                  foreground = log2(res_deseq[sig_genes_id, "baseMean"]), 
                  background = log2(res_deseq[backG, "baseMean"])), 
             xlab = "log2 mean normalized counts",
             main = "Matching for enrichment analysis") 

## Selecting genes 
inUniverse <- all_genes_id %in% c(sig_genes_id, backG)
inSelection <- all_genes_id %in% sig_genes_id
alg <- factor(as.integer(inSelection[inUniverse]))
names(alg) <- all_genes_id[inUniverse]

## Preparing data
tgd <- new("topGOdata", ontology = "BP", allGenes = alg, nodeSize = 5, annot = annFUN.org,
           mapping = "org.Hs.eg.db", ID = "ensembl") # Change to MF and CC


## Running tests
resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher")

if (length(nodes(graph(tgd))) < 200) {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim", numChar=1000, topNodes = length(nodes(graph(tgd))))
} else {
  tab <- GenTable(tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,
                  orderBy = "Fisher.elim",numChar=1000, topNodes = 200)
}

printGraph(tgd, resultTopGO.elim, firstSigNodes = 5, 
           fn.prefix = "Output_files/GO/Remission_crp/aza_25varPart_topvar/aza_R_crp_elim_NoBackground_BP_pos",
           useInfo = "all", pdfSW = TRUE) 


## Finding the genes associated with the enriched GO terms
for (i in 1:nrow(tab)) {
  go_id <- as.vector(tab[i, 1])
  genes_in_cat <- intersect(unlist(genesInTerm(tgd, go_id)), sig_genes_id)
  gene_sym_in_cat <- unique_ensg2gene[as.character(genes_in_cat), "hgnc_symbol"]
  gene_sym_in_cat_str <- ""
  
  if (length(genes_in_cat) > 0) {
    for (j in 1:length(gene_sym_in_cat)) {
      gene_sym_in_cat_str <- paste(gene_sym_in_cat_str, gene_sym_in_cat[j],
                                   sep = ",")
    }
  }
  if (gene_sym_in_cat_str == "") {
    gene_sym_in_cat_str <- NA
  }
  
  tab$Genes[i] <- gene_sym_in_cat_str
  
}

tab_aza <- copy(tab)

## Writing tables
tab_aza <- tab_aza %>% 
  mutate(geneRatio = tab_aza[,4]/tab_aza[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_aza,"Output_files/GO/Remission_crp/aza_25varPart_topvar/TopGO_BP_aza_R_crp_varPart_topvar_pos.txt", sep="\t",row.names=TRUE)

### significant terms
tab_aza_sig <- tab_aza %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_aza_sig$Term = with(tab_aza_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_aza_sig,"Output_files/GO/Remission_crp/aza_25varPart_topvar/TopGO_BP_sigTerms_aza_R_crp_varPart_topvar_pos.txt", sep="\t",row.names=TRUE)

## Barplot
tab_aza_sig_plot <- tab_aza_sig#[1:30,]
gg <- ggplot(tab_aza_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="darkorange", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("") +
  labs(
    title = "GO BP Aza x No systemic therapies (upregulated)",
    subtitle = "Intersection: VarPart > 25% + top var"
  ) +
  theme(
    plot.title = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 12, vjust = 1, hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = -15))
  )
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission_crp/aza_25varPart_topvar/TopGO_BP_aza_R_crp_varPart_topvar_pos.png", height = 8, width = 10, units = "in", dpi = 300)







