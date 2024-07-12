# EZE cohort 
# Gene ontology of significant genes from model (maaslin2)
# Anti TNF vs No Biologics
  # Date 05.04


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

folder <- "Output_files/GO/All_patients/antiTNF_noBio"

if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading data ----
sig_genes <- read.csv("Output_files/Maaslin2/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_vst_q0.05.txt",
                      sep = "\t")

coldata <- read.csv("Output_files/Maaslin2/Tables/coldata_maaslin2_antiTNF_vs_noBiologics_diag_20.03.txt", sep = "\t")
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

write.table(res_deseq,"Output_files/DESeq2/TopGO/DESeq2_baseMean_antiTNF_noBio_31.03.txt", sep="\t", row.names=TRUE)


# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/DESeq2_baseMean_antiTNF_noBio_31.03.txt", sep="\t")

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
           fn.prefix = "Output_files/GO/All_patients/antiTNF_noBio/antiTNF_noBio_elim_NoBackground_BP_neg",
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

tab_tnf <- copy(tab)


## Writing tables
### all terms
tab_tnf <- tab_tnf %>% 
  mutate(geneRatio = tab_tnf[,4]/tab_tnf[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_tnf,"Output_files/GO/All_patients/antiTNF_noBio/TopGO_antiTNF_noBio_BP_neg.txt", sep="\t",row.names=TRUE)

### significant terms
tab_tnf_sig <- tab_tnf %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_tnf_sig$Term = with(tab_tnf_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_tnf_sig,"Output_files/GO/All_patients/antiTNF_noBio/TopGO_significantTerms_antiTNF_noBio_BP_neg.txt", sep="\t",row.names=TRUE)

## Barplot
tab_tng_sig_plot <- tab_tnf_sig[1:30,]
gg <- ggplot(tab_tng_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="skyblue4", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("GO Term") + ggtitle("GO Anti TNF x no Biologics BP (coeff < 0)")
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/All_patients/antiTNF_noBio/TopGO_antiTNF_noBio_BP_neg.png", height = 8, width = 10, units = "in", dpi = 300)






###################################################### Remitters ######################################################################
graphics.off()
rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission/antiTNF_noBio"

if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading data ----
sig_genes <- read.csv("Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_R_vst_q0.05.txt",
                      sep = "\t")

coldata <- read.csv("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_antiTNF_vs_noBio_diag_R_20.03.txt", sep = "\t")
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

write.table(res_deseq,"Output_files/DESeq2/TopGO/Remission/DESeq2_baseMean_antiTNF_noBio_R_31.03.txt", sep="\t", row.names=TRUE)


# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission/DESeq2_baseMean_antiTNF_noBio_R_31.03.txt", sep="\t")

res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) 
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
           fn.prefix = "Output_files/GO/Remission/antiTNF_noBio/antiTNF_noBio_elim_NoBackground_BP_pos",
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

tab_tnf <- copy(tab)


## Writing tables
### all terms
tab_tnf <- tab_tnf %>% 
  mutate(geneRatio = tab_tnf[,4]/tab_tnf[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_tnf,"Output_files/GO/Remission/antiTNF_noBio/TopGO_antiTNF_noBio_BP_pos.txt", sep="\t",row.names=TRUE)

### significant terms
tab_tnf_sig <- tab_tnf %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_tnf_sig$Term = with(tab_tnf_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_tnf_sig,"Output_files/GO/Remission/antiTNF_noBio/TopGO_significantTerms_antiTNF_noBio_BP_pos.txt", sep="\t",row.names=TRUE)

## Barplot
tab_tng_sig_plot <- tab_tnf_sig#[1:30,]
gg <- ggplot(tab_tng_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="darkorange", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("GO Term") + ggtitle("GO Anti TNF x no Biologics BP - remitters (coef > 0)")
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission/antiTNF_noBio/TopGO_antiTNF_noBio_BP_pos.png", height = 8, width = 10, units = "in", dpi = 300)






###################################################### Intersection genes (topVar + Remitters)
graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Intersections/antiTNF_noBio"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
list_sig_genes_all_remitters <- readRDS("Output_files/Maaslin2/significant_results/list_sig_genes_all_remitters_intersect_diff_06.04.rds")
int_R_noBio_topvar <- list_sig_genes_all_remitters$intersect_R_noBio_topVar

sig_genes <- read.csv("Output_files/Maaslin2/Results/Remission/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_R_vst_q0.05.txt",
                      sep = "\t")

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Subsetting sig genes to have only intersecting genes 
int_sig_genes <- sig_genes[sig_genes$feature %in% int_R_noBio_topvar,]

# Top GO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission/DESeq2_baseMean_antiTNF_noBio_R_31.03.txt", sep="\t")

res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) 
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
           fn.prefix = "Output_files/GO/Intersections/antiTNF_noBio/R_antiTNF_noBio_elim_NoBackground_BP_neg",
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

tab_tnf <- copy(tab)


## Writing tables
### all terms
tab_tnf <- tab_tnf %>% 
  mutate(geneRatio = tab_tnf[,4]/tab_tnf[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_tnf,"Output_files/GO/Intersections/antiTNF_noBio/TopGO_R_antiTNF_noBio_BP_neg.txt", sep="\t",row.names=TRUE)

### significant terms
tab_tnf_sig <- tab_tnf %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_tnf_sig$Term = with(tab_tnf_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_tnf_sig,"Output_files/GO/Intersections/antiTNF_noBio/TopGO_R_significantTerms_antiTNF_noBio_BP_neg.txt", sep="\t",row.names=TRUE)

## Barplot
tab_tnf_sig_plot <- tab_tnf_sig#[1:30,]
gg <- ggplot(tab_tnf_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="skyblue4", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("") +
  labs(
    title = "GO Anti TNF x no Biologics BP (coef < 0)",
    subtitle = "Intersection: Inactive + TopVar genes"
  ) +
  theme(
    plot.title = element_text(angle = 0, size = 20, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 15, vjust = 1, hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = -15))
  )
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Intersections/antiTNF_noBio/TopGO_R_antiTNF_noBio_BP_neg.png", height = 8, width = 10, units = "in", dpi = 300)






######################### Remitters + crp ###############################
graphics.off()
rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission_crp/antiTNF_noBio"

if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading data ----
sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_R_crp.txt", sep = "\t")

coldata <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_antiTNF_R_crp_04.05.txt", sep = "\t")
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

write.table(res_deseq,"Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_antiTNF_R_crp.txt", sep="\t", row.names=TRUE)


# TopGO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_antiTNF_R_crp.txt", sep= "\t")

res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) 
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
           fn.prefix = "Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_elim_NoBackground_BP_neg",
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

tab_tnf <- copy(tab)


## Writing tables
tab_tnf <- tab_tnf %>% 
  mutate(geneRatio = tab_tnf[,4]/tab_tnf[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_tnf,"Output_files/GO/Remission_crp/antiTNF_noBio/TopGO_BP_antiTNF_R_crp_neg.txt", sep="\t",row.names=TRUE)

### significant terms
tab_tnf_sig <- tab_tnf %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_tnf_sig$Term = with(tab_tnf_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_tnf_sig,"Output_files/GO/Remission_crp/antiTNF_noBio/TopGO_BP_significantTerms_antiTNF_R_crp_neg.txt", sep="\t",row.names=TRUE)

## Barplot
tab_tng_sig_plot <- tab_tnf_sig#[1:30,]
gg <- ggplot(tab_tng_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="skyblue4", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("GO Term") + ggtitle("GO BP Anti TNF x no Biologics - remitters + crp (downregulated)")
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission_crp/antiTNF_noBio/TopGO_BP_antiTNF_R_crp_neg.png", height = 8, width = 10, units = "in", dpi = 300)




###################################################### Intersection genes (topVar + Remitters)
graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_topVar"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_R_crp.txt", sep = "\t")
topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")

intersect_sig_topvar <- base::intersect(sig_genes$feature, rownames(topVar))
sig_topVar <- sig_genes[sig_genes$feature %in% intersect_sig_topvar,]

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Top GO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_antiTNF_R_crp.txt", sep= "\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) 

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
           fn.prefix = "Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_topVar/intersect_antiTNF_R_crp_topVar_elim_NoBackground_BP_pos",
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

tab_tnf <- copy(tab)


## Writing tables
### all terms
tab_tnf <- tab_tnf %>% 
  mutate(geneRatio = tab_tnf[,4]/tab_tnf[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_tnf,"Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_topVar/TopGO_BP_intersection_antiTNF_R_crp_pos.txt", sep="\t",row.names=TRUE)

### significant terms
tab_tnf_sig <- tab_tnf %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_tnf_sig$Term = with(tab_tnf_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_tnf_sig,"Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_topVar/TopGO_BP_sigTerms_intersection_antiTNF_R_crp_pos.txt", sep="\t",row.names=TRUE)

## Barplot
tab_tnf_sig_plot <- tab_tnf_sig#[1:30,]
gg <- ggplot(tab_tnf_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="darkorange", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("") +
  labs(
    title = "GO Anti TNF x no Biologics BP (upregulated)",
    subtitle = "Intersection: Inactive (crp) + TopVar genes"
  ) +
  theme(
    plot.title = element_text(angle = 0, size = 20, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 15, vjust = 1, hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = -15))
  )
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_topVar/TopGO_BP_intersect_antiTNF_R_crp_pos.png", height = 8, width = 10, units = "in", dpi = 300)


# other dot plot
go_tnf_down <- read.table("Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_topVar/TopGO_BP_sigTerms_intersection_antiTNF_R_crp_neg.txt", sep="\t")
go_tnf_down$Analysis <- "Downregulated"
go_tnf_up <- read.table("Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_topVar/TopGO_BP_sigTerms_intersection_antiTNF_R_crp_pos.txt", sep="\t")
go_tnf_up$Analysis <- "Upregulated"


all_terms_tnf <- rbind(go_tnf_up[c(1,2,4),], go_tnf_down[1:6,])
all_terms_tnf  <- all_terms_tnf %>%
  group_by(Analysis) %>% 
  mutate(Term = fct_reorder2(Term, geneRatio,-geneRatio)) %>%
  ungroup()

pdf("Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_topVar/reddotplot_GO_BP_antiTNF_topVar_selectedTerms.pdf", width = 7, height = 7)
dot_plot <- ggplot(all_terms_tnf, aes(x = Analysis, y= Term, size = geneRatio, color = Fisher.elim)) +
  geom_point() +
  scale_colour_gradient(low = "#990000", high = "#FF9999") +
  scale_size_continuous(range=c(3,8)) +
  xlab('') + ylab('') +
  labs(
    title = "GO BP Anti-TNF x No Bio",
    subtitle = "Inactive disease (crp) sig. genes in TVGs",
    caption = paste0("Produced on ", Sys.time()),
    color = "p-value",
    size = "Gene ratio") +
  theme_bw() + #base_size = 20
  theme(
    axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
    axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1),
    axis.title = element_text(size = 20)
  )

dot_plot
dev.off()





###################################################### Unique genes (DEGS not part of TVGs)
graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

folder <- "Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_not_topVar"

if (!dir.exists(folder)) {
  dir.create(folder)
}


# Loading data ----
sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_R_crp.txt", sep = "\t")
topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")

unique_sig <- setdiff(sig_genes$feature, rownames(topVar))

unique_sig <- sig_genes[sig_genes$feature %in% unique_sig,]

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# Top GO ----
res_deseq <- read.table("Output_files/DESeq2/TopGO/Remission_crp/DESeq2_baseMean_antiTNF_R_crp.txt", sep= "\t")
res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) 

sig_genes_id <- subset(unique_sig, unique_sig$coef < 0) #modify coef
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
           fn.prefix = "Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_not_topVar/unique_antiTNF_R_crp_topVar_elim_NoBackground_BP_neg",
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

tab_tnf <- tab


## Writing tables
### all terms
tab_tnf <- tab_tnf %>% 
  mutate(geneRatio = tab_tnf[,4]/tab_tnf[,3]) %>% 
  relocate(geneRatio, .after = Expected)

write.table(tab_tnf,"Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_not_topVar/TopGO_BP_unique_antiTNF_R_crp_neg.txt", sep="\t",row.names=TRUE)

### significant terms
tab_tnf_sig <- tab_tnf %>% 
  mutate(across(Fisher.elim, as.numeric)) %>% 
  filter(Fisher.elim < 0.05)

tab_tnf_sig$Term = with(tab_tnf_sig, reorder(Term, Fisher.elim, decreasing = TRUE))

write.table(tab_tnf_sig,"Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_not_topVar/TopGO_BP_sigTerms_unique_antiTNF_R_crp_neg.txt", sep="\t",row.names=TRUE)

## Barplot
tab_tnf_sig_plot <- tab_tnf_sig#[1:30,]
gg <- ggplot(tab_tnf_sig_plot)
gg <- gg + geom_bar(aes(x = -log10(Fisher.elim), y = Term), position = "stack", stat = "identity", 
                    fill="skyblue4", alpha=.6, width=.4) # "darkorange" or "skyblue4"
gg <- gg + theme_bw() + theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"))
gg <- gg + theme(axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 1)))
gg <- gg + xlab("-log10 (p-value)") + ylab("") +
  labs(
    title = "GO Anti TNF x no Biologics BP (downregulated)",
    subtitle = "Intersection: Inactive (crp) not in TVGs"
  ) +
  theme(
    plot.title = element_text(angle = 0, size = 20, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 15, vjust = 1, hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = -15))
  )
gg <- gg + geom_text(aes(-log10(Fisher.elim), y = Term, label = round(geneRatio,2), hjust = -0.2), size = 4)
gg
ggsave(gg, file = "Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_not_topVar/TopGO_BP_unique_antiTNF_R_crp_neg.png", height = 8, width = 10, units = "in", dpi = 300)


# other dot plot
go_tnf_down <- read.table("Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_not_topVar/TopGO_BP_sigTerms_unique_antiTNF_R_crp_neg.txt", sep="\t")
go_tnf_down$Analysis <- "Downregulated"
go_tnf_up <- read.table("Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_not_topVar/TopGO_BP_sigTerms_unique_antiTNF_R_crp_pos.txt", sep="\t")
go_tnf_up$Analysis <- "Upregulated"

  
all_terms_tnf <- rbind(go_tnf_up, go_tnf_down)
all_terms_tnf  <- all_terms_tnf %>%
  group_by(Analysis) %>% 
  mutate(Term = fct_reorder2(Term, geneRatio,-geneRatio)) %>%
  ungroup()

pdf("Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_not_topVar/dotplot_GO_BP_unique_antiTNF.pdf", width = 6, height = 7)
dot_plot <- ggplot(all_terms_tnf, aes(x = Analysis, y= Term, size = geneRatio, color = Fisher.elim)) +
  geom_point() +
  scale_colour_gradient(low = "#990000", high = "#FF9999") +
  scale_size_continuous(range=c(3,8)) +
  xlab('') + ylab('') +
  labs(
    title = "GO BP Anti-TNF x No Bio",
    subtitle = "Inactive disease (crp) sig. genes not in TVGs",
    caption = paste0("Produced on ", Sys.time()),
    color = "p-value",
    size = "Gene ratio") +
  theme_bw() + #base_size = 20
  theme(
    axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
    axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1),
    axis.title = element_text(size = 20)
  )

dot_plot
dev.off()



### Same dot plot for shared (with TVGs) and unique DEGs

go_tnf_down <- read.table("Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_topVar/TopGO_BP_sigTerms_intersection_antiTNF_R_crp_neg.txt", sep="\t")
go_tnf_down$Direction <- "Downregulated"
go_tnf_down$Analysis <- "Sig genes + TVGs"

go_tnf_up <- read.table("Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_topVar/TopGO_BP_sigTerms_intersection_antiTNF_R_crp_pos.txt", sep="\t")
go_tnf_up$Direction <- "Upregulated"
go_tnf_up$Analysis <- "Sig genes + TVGs"


unique_go_tnf_down <- read.table("Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_not_topVar/TopGO_BP_sigTerms_unique_antiTNF_R_crp_neg.txt", sep="\t")
unique_go_tnf_down$Direction <- "Downregulated"
unique_go_tnf_down$Analysis <- "Unique sig genes"

unique_go_tnf_up <- read.table("Output_files/GO/Remission_crp/antiTNF_noBio/antiTNF_not_topVar/TopGO_BP_sigTerms_unique_antiTNF_R_crp_pos.txt", sep="\t")
unique_go_tnf_up$Direction <- "Upregulated"
unique_go_tnf_up$Analysis <- "Unique sig genes"


all_terms_tnf <- rbind(go_tnf_up[c(1,2,4),], go_tnf_down[1:6,], unique_go_tnf_down, unique_go_tnf_up)
all_terms_tnf  <- all_terms_tnf %>%
  group_by(Analysis) %>% 
  mutate(Term = fct_reorder2(Term, geneRatio,-geneRatio)) %>%
  ungroup()

pdf("Output_files/GO/Remission_crp/antiTNF_noBio/dotplot_GO_BP_unique_shared_antiTNF.pdf", width = 8, height = 7)
dot_plot <- ggplot(all_terms_tnf, aes(x = Direction, y= Term, size = geneRatio, color = Fisher.elim)) +
  geom_point() +
  scale_colour_gradient(low = "#990000", high = "#FF9999") +
  scale_size_continuous(range=c(3,8)) +
  xlab('') + ylab('') +
  labs(
    title = "GO BP Anti-TNF x No Bio",
    subtitle = "Unique and shared DEGs",
    caption = paste0("Produced on ", Sys.time()),
    color = "p-value",
    size = "Gene ratio") +
  theme_bw() + #base_size = 20
  theme(
    axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
    axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1),
    axis.title = element_text(size = 20)
  ) + facet_wrap(~Analysis)

dot_plot
dev.off()
