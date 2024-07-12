# EZE cohort 
# Gene ontology of significant genes from model (maaslin2)
# Prednisolon x No Systemic therapies
  # Date: 10.06


graphics.off()
rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC
# Loading packages ----
if (!requireNamespace("topGO", quietly = TRUE)) BiocManager::install("topGO")
if (!requireNamespace("genefilter", quietly = TRUE)) BiocManager::install("genefilter")
if (!requireNamespace("geneplotter", quietly = TRUE)) BiocManager::install("geneplotter")

library(genefilter)
library(geneplotter)
library(topGO)
library(data.table)
library(tidyverse)
library(ggrepel)
library(ggthemes)
library(DESeq2)

folder <- "Output_files/GO/All_patients/pred_noSyst"

if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading data ----
sig_genes <- read.csv("Output_files/Maaslin2/significant_results/maaslin2_significant_results_pred_vs_noSyst_vst_q0.05.txt",
                      sep = "\t")

coldata <- read.csv("Output_files/Maaslin2/Tables/coldata_maaslin2_pred_vs_noSyst_16.03.txt", sep = "\t")
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
res_deseq <- as.data.frame(res_deseq)

res_deseq$gene <- unique_ensg2gene[rownames(res_deseq), ]$hgnc_symbol
res_deseq$gene_id <- rownames(res_deseq)

write.table(res_deseq,"Output_files/DESeq2/TopGO/DESeq2_baseMean_pred_noSyst_31.03.txt", sep="\t", row.names=TRUE)

# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/DESeq2_baseMean_pred_noSyst_31.03.txt", sep="\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) # Filtering out genes that were never detected 

sig_genes_id <- subset(sig_genes, sig_genes$coef < 0) #modify coef
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
           fn.prefix = "Output_files/GO/All_patients/pred_noSyst/pred_noSyst_elim_NoBackground_BP_neg",
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

tab_pred <- copy(tab)

## Writing tables
### all terms
tab_pred <- tab_pred %>% 
  mutate(geneRatio = tab_pred[,4]/tab_pred[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_pred,"Output_files/GO/All_patients/pred_noSyst/TopGO_pred_noSyst_BP_neg.txt", sep="\t",row.names=TRUE)

### significant terms
tab_pred_sig <- tab_pred %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_pred_sig$Term = with(tab_pred_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_pred_sig,"Output_files/GO/All_patients/pred_noSyst/TopGO_significantTerms_pred_noSyst_BP_neg.txt", sep="\t",row.names=TRUE)

## Barplot
tab_pred_sig_plot <- tab_pred_sig[1:30,]
gg <- ggplot(tab_pred_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term),
                    position = "stack", stat = "identity", fill="skyblue4", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("GO Term") + ggtitle("GO Prednisolon x no Systemic Therapies (coeff < 0)")
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/All_patients/pred_noSyst/TopGO_pred_noSyst_BP_neg.png", height = 8, width = 10, units = "in", dpi = 300)





###################################################### Remitters ######################################################################
graphics.off()
rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission/pred_noSyst"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
sig_genes <- read.csv("Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_vst_q0.05.txt",
                      sep = "\t")

coldata <- read.csv("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_pred_vs_noSyst_R_17.03.txt", sep = "\t")
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
res_deseq <- as.data.frame(res_deseq)

res_deseq$gene <- unique_ensg2gene[rownames(res_deseq), ]$hgnc_symbol
res_deseq$gene_id <- rownames(res_deseq)

write.table(res_deseq,"Output_files/DESeq2/TopGO/Remission/DESeq2_results_pred_noSyst_R_31.03.txt", sep="\t", row.names=TRUE)

# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission/DESeq2_results_pred_noSyst_R_31.03.txt", sep="\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) # Filtering out genes that were never detected 

sig_genes_id <- subset(sig_genes, sig_genes$coef < 0) #modify coef
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
           fn.prefix = "Output_files/GO/Remission/pred_noSyst/pred_noSyst_elim_NoBackground_BP_neg",
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

tab_pred <- copy(tab)

## Writing tables
### all terms
tab_pred <- tab_pred %>% 
  mutate(geneRatio = tab_pred[,4]/tab_pred[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_pred,"Output_files/GO/Remission/pred_noSyst/TopGO_pred_noSyst_BP_neg.txt", sep="\t",row.names=TRUE)

### significant terms
tab_pred_sig <- tab_pred %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_pred_sig$Term = with(tab_pred_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_pred_sig,"Output_files/GO/Remission/pred_noSyst/TopGO_significantTerms_pred_noSyst_BP_neg.txt", sep="\t",row.names=TRUE)

## Barplot
tab_pred_sig_plot <- tab_pred_sig[1:30,]
gg <- ggplot(tab_pred_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term),
                    position = "stack", stat = "identity", fill="skyblue4", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("GO Term") + ggtitle("GO Prednisolon x no Systemic Therapies - remitters (coeff < 0)")
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission/pred_noSyst/TopGO_pred_noSyst_BP_neg.png", height = 8, width = 10, units = "in", dpi = 300)




##################################################### Intersection genes (topVar + Remitters)
graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Intersections/pred_noSyst"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
list_sig_genes_all_remitters <- readRDS("Output_files/Maaslin2/significant_results/list_sig_genes_all_remitters_intersect_diff_06.04.rds")
int_R_pred_topvar <- list_sig_genes_all_remitters$intersect_R_pred_topVar

sig_genes <- read.csv("Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_vst_q0.05.txt",
                      sep = "\t")

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Subsetting sig genes to have only intersecting genes 
int_sig_genes <- sig_genes[sig_genes$feature %in% int_R_pred_topvar,]


# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission/DESeq2_results_pred_noSyst_R_31.03.txt", sep="\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) # Filtering out genes that were never detected 

sig_genes_id <- subset(int_sig_genes, int_sig_genes$coef < 0) #modify coef
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
           fn.prefix = "Output_files/GO/Intersections/pred_noSyst/R_pred_noSyst_elim_NoBackground_BP_neg",
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

tab_pred <- copy(tab)

## Writing tables
### all terms
tab_pred <- tab_pred %>% 
  mutate(geneRatio = tab_pred[,4]/tab_pred[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_pred,"Output_files/GO/Intersections/pred_noSyst/TopGO_R_pred_noSyst_BP_neg.txt", sep="\t",row.names=TRUE)

### significant terms
tab_pred_sig <- tab_pred %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_pred_sig$Term = with(tab_pred_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_pred_sig,"Output_files/GO/Intersections/pred_noSyst/TopGO_R_significantTerms_pred_noSyst_BP_neg.txt", sep="\t",row.names=TRUE)

## Barplot
tab_pred_sig_plot <- tab_pred_sig[1:30,]
gg <- ggplot(tab_pred_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="skyblue4", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("") +
  labs(
    title = "GO Prednisolon x no Systemic Therapies (coef < 0)",
    subtitle = "Intersection: Inactive + TopVar genes"
  ) +
  theme(
    plot.title = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 12, vjust = 1, hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = -15))
  )
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Intersections/pred_noSyst/TopGO_R_pred_noSyst_BP_neg.png", height = 8, width = 10, units = "in", dpi = 300)








################# Remitters + crp ###################

rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission_crp/pred_noSyst"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
sig_genes <- read.table( "Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt", sep = "\t")
coldata <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_pred_R_crp_04.05.txt", sep = "\t")
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
res_deseq <- as.data.frame(res_deseq)

res_deseq$gene <- unique_ensg2gene[rownames(res_deseq), ]$hgnc_symbol
res_deseq$gene_id <- rownames(res_deseq)

write.table(res_deseq,"Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_pred_R_crp.txt", sep="\t", row.names=TRUE)

# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_pred_R_crp.txt", sep="\t")
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
           fn.prefix = "Output_files/GO/Remission_crp/pred_noSyst/pred_noSyst_elim_NoBackground_BP_pos",
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

tab_pred <- copy(tab)

## Writing tables
tab_pred <- tab_pred %>% 
  mutate(geneRatio = tab_pred[,4]/tab_pred[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_pred,"Output_files/GO/Remission_crp/pred_noSyst/TopGO_BP_pred_R_crp_pos.txt", sep="\t",row.names=TRUE)

### significant terms
tab_pred_sig <- tab_pred %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_pred_sig$Term = with(tab_pred_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_pred_sig,"Output_files/GO/Remission_crp/pred_noSyst/TopGO_BP_sigTerms_pred_R_crp_pos.txt", sep="\t",row.names=TRUE)

## Barplot
tab_pred_sig_plot <- tab_pred_sig[1:30,]
gg <- ggplot(tab_pred_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term),
                    position = "stack", stat = "identity", fill="darkorange", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("GO Term") + ggtitle("GO Prednisolon x no Systemic Therapies - remitters + crp (upregulated")
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission_crp/pred_noSyst/TopGO_BP_pred_R_crp_pos.png", height = 8, width = 10, units = "in", dpi = 300)





################# Intersection: topVar + (Remitters + crp) ###################
graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission_crp/pred_noSyst/pred_topVar"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
sig_genes <- read.table( "Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt", sep = "\t")
topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")

intersect_sig_topvar <- base::intersect(sig_genes$feature, rownames(topVar))

sig_topVar <- sig_genes[sig_genes$feature %in% intersect_sig_topvar,]

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_pred_R_crp.txt", sep="\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) # Filtering out genes that were never detected 

sig_genes_id <- subset(sig_topVar, sig_topVar$coef > 0) #modify coef
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
           fn.prefix = "Output_files/GO/Remission_crp/pred_noSyst/pred_topVar/R_pred_noSyst_elim_NoBackground_BP_pos",
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

tab_pred <- copy(tab)

## Writing tables
### all terms
tab_pred <- tab_pred %>% 
  mutate(geneRatio = tab_pred[,4]/tab_pred[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_pred,"Output_files/GO/Remission_crp/pred_noSyst/pred_topVar/TopGO_BP_intersect_pred_R_crp_topVar_pos.txt", sep="\t",row.names=TRUE)

### significant terms
tab_pred_sig <- tab_pred %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_pred_sig$Term = with(tab_pred_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_pred_sig,"Output_files/GO/Remission_crp/pred_noSyst/pred_topVar/TopGO_BP_sigTerms_intersect_pred_R_crp_topVar_pos.txt", sep="\t",row.names=TRUE)

## Barplot
tab_pred_sig_plot <- tab_pred_sig[1:30,]
gg <- ggplot(tab_pred_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="darkorange", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("") +
  labs(
    title = "GO Pred x No systemic therapies (upregulated)",
    subtitle = "Intersection: Inactive (crp) + TopVar genes"
  ) +
  theme(
    plot.title = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 12, vjust = 1, hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = -15))
  )
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission_crp/pred_noSyst/pred_topVar/TopGO_BP_intersect_pred_R_crp_topVar_pos.png", height = 8, width = 10, units = "in", dpi = 300)






################# Remitters + crp (LFC 0.5) ###################
graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
sig_genes <- read.table( "Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt", sep = "\t")

sig_genes_lfc0.5 <- sig_genes %>% 
  filter(abs(coef) > 0.5)

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_pred_R_crp.txt", sep="\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) # Filtering out genes that were never detected 

sig_genes_id <- subset(sig_genes_lfc0.5, sig_genes_lfc0.5$coef > 0) #modify coef

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
           fn.prefix = "Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5/pred_lfc0.5_elim_NoBackground_BP_pos",
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

tab_pred <- copy(tab)

## Writing tables
### all terms
tab_pred <- tab_pred %>% 
  mutate(geneRatio = tab_pred[,4]/tab_pred[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_pred,"Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5/TopGO_BP_pred_lfc0.5_up.txt", sep="\t",row.names=TRUE)

### significant terms
tab_pred_sig <- tab_pred %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_pred_sig$Term = with(tab_pred_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_pred_sig,"Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5/TopGO_BP_sig_pred_lfc0.5_up.txt", sep="\t",row.names=TRUE)

## Barplot
tab_pred_sig_plot <- tab_pred_sig[1:30,]
gg <- ggplot(tab_pred_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="darkorange", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("") +
  labs(
    title = "GO Pred x No systemic therapies (upregulated)",
    subtitle = "Inactive (crp) lfc 0.5 "
  ) +
  theme(
    plot.title = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 12, vjust = 1, hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = -15))
  )
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5/TopGO_BP_pred_lfc0.5_up.png", height = 8, width = 10, units = "in", dpi = 300)




################# Intersection: topVar + (Remitters + crp LFC 0.5) ###################
graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5_topVar"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
sig_genes <- read.table( "Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt", sep = "\t")
topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")

sig_genes_lfc0.5 <- sig_genes %>% 
  filter(abs(coef) > 0.5)

sig_genes_TVG <- intersect(rownames(topVar),sig_genes_lfc0.5$feature )
sig_genes_lfc0.5_TVG <- sig_genes_lfc0.5[sig_genes_lfc0.5$feature %in% sig_genes_TVG,]

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_pred_R_crp.txt", sep="\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) # Filtering out genes that were never detected 

sig_genes_id <- subset(sig_genes_lfc0.5_TVG, sig_genes_lfc0.5_TVG$coef < 0) #modify coef

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
           fn.prefix = "Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5_topVar/pred_lfc0.5_elim_NoBackground_BP_neg",
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

tab_pred <- copy(tab)

## Writing tables
### all terms
tab_pred <- tab_pred %>% 
  mutate(geneRatio = tab_pred[,4]/tab_pred[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_pred,"Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5_topVar/TopGO_BP_pred_lfc0.5_topVar_down.txt", sep="\t",row.names=TRUE)

### significant terms
tab_pred_sig <- tab_pred %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_pred_sig$Term = with(tab_pred_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_pred_sig,"Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5_topVar/TopGO_BP_sig_pred_lfc0.5_topVar_down.txt", sep="\t",row.names=TRUE)

## Barplot
tab_pred_sig_plot <- tab_pred_sig[1:30,]
gg <- ggplot(tab_pred_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="skyblue4", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("") +
  labs(
    title = "GO Pred x No systemic therapies (downregulated)",
    subtitle = "Intersection: Inactive (crp) lfc 0.5 + TopVar genes"
  ) +
  theme(
    plot.title = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 12, vjust = 1, hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = -15))
  )
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5_topVar/TopGO_BP_pred_lfc0.5_topVar_down.png", height = 8, width = 10, units = "in", dpi = 300)


## Dotplot up and down
go_pred_up <- read.table("Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5_topVar/TopGO_BP_sig_pred_lfc0.5_topVar_up.txt", sep="\t")
go_pred_down <- read.table("Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5_topVar/TopGO_BP_sig_pred_lfc0.5_topVar_down.txt", sep="\t")

terms <- c("positive regulation of inflammatory response", "defense response to bacterium", "acute inflammatory response", 
           "positive regulation of cytokine production", "regulation of macrophage activation", "inflammatory response",
           "response to cAMP", "regulation of macrophage activation", "positive regulation of NF-kappaB transcription factor activity",
           "neutrophil extravasation", "regulation of toll-like receptor signaling pathway", "regulation of lipid storage",
           
           "adaptive immune response", "cell surface receptor signaling pathway", "B cell receptor signaling pathway",
           "T cell activation", "cell population proliferation", "leukocyte differentiation", "leukocyte proliferation", 
           "cell population proliferation")


plot_up_go_result <- go_pred_up %>% 
  mutate(plot = ifelse(Term %in% terms, "Yes", "No")) %>% 
  filter(plot == "Yes") %>% 
  mutate(direction = "Upregulated")

plot_down_go_result <- go_pred_down %>% 
  mutate(plot = ifelse(Term %in% terms, "Yes", "No")) %>% 
  filter(plot == "Yes") %>% 
  mutate(direction = "Downregulated")

selected_GO <- rbind(plot_up_go_result, plot_down_go_result)
selected_GO$Term <- Term(as.vector(selected_GO$GO.ID))
selected_GO$Term <- factor(selected_GO$Term, levels = rev(selected_GO$Term)) # fixes order

pdf("Output_files/GO/Remission_crp/pred_noSyst/pred_LFC0.5_topVar/dotplot_GO_BP_pred_lfc0.5_topVar.pdf", width = 7, height = 8)
dot_plot <- ggplot(selected_GO, aes(x = direction, y= Term, size = geneRatio, color = Fisher.elim)) +
  geom_point() +
  scale_colour_gradient(low = "#990000", high = "#FF9999") +
  #scale_size_continuous(range=c(3,10)) +
  xlab('') + ylab('') +
  labs(
    title = "GO BP Pred x No Syst",
    subtitle = "Sig genes lfc 0.5 + topVar",
    caption = paste0("Produced on ", Sys.time()),
    color = "p-value",
    size = "Gene ratio") +
  theme_bw() + #base_size = 20
  theme(
    axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
    axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1),
    axis.title = element_text(size = 20)
  ) +
  scale_y_discrete(limits = c("adaptive immune response",  "B cell receptor signaling pathway",
                              "T cell activation", "cell population proliferation", "leukocyte differentiation", "leukocyte proliferation",
                              
                              "defense response to bacterium","cell surface receptor signaling pathway",
                              
                              "positive regulation of inflammatory response", "acute inflammatory response", "regulation of lipid storage",
                              "positive regulation of cytokine production", "regulation of macrophage activation", "inflammatory response",
                              "response to cAMP", "regulation of macrophage activation", "positive regulation of NF-kappaB transcription factor activity",
                              "neutrophil extravasation", "regulation of toll-like receptor signaling pathway", "cell population proliferation"))
dot_plot
dev.off()
