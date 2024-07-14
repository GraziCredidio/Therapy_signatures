# EZE cohort: Therapy Signatures
  # GO analysis of DEGs

  # Anti-TNF: DEGs in TVG and DEGs not in topVar. Figure 5B
  # Pred: DEGs with absLFC > 0.5 in topVar. Figure 11
  # Aza: DEGs with varPart > 25% in topVar. Figure 16

  # Author: Graziella Credidio


rm(list = ls())

# Loading packages ----
library(DESeq2)
library(genefilter)
library(tidyverse)
library(topGO)
library(data.table)

# Crating folders ----
folders <- c("Output_files/GO/aza", "Output_files/GO/pred", "Output_files/GO/antiTNF")

for (i in folders){
  if (!dir.exists(i)) {
    dir.create(i)}
}
  
# Loading data ----
topVar <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVar_genes <- rownames(topVar)

ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

counts <- read.csv("Cleaned_tables/EZECohort_counts_ord.txt", sep="\t")

# Aza
aza_coldata <- read.table("Cleaned_tables/models/aza/EZECohort_coldata_aza.txt", sep = "\t")
aza_top25percent_varPart <- read.table("Output_files/Variance_partition/aza/top25percent_varPart_aza.txt")
aza_degs <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_aza_vs_noSyst.txt", sep = '\t')

# Pred
pred_coldata <- read.table("Cleaned_tables/models/pred/EZECohort_coldata_pred.txt", sep = "\t")
pred_degs <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_pred_vs_noSyst.txt", sep = '\t')

# Anti-TNF
antiTNF_coldata <- read.table("Cleaned_tables/models/antiTNF/EZECohort_coldata_antiTNF.txt", sep = "\t")
antiTNF_degs <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_antiTNF_vs_noBio.txt", sep = '\t')

# Getting set of genes and splitting into up and downregulated genes ----
# Aza
aza_degs_topVar_varPart25 <- aza_degs[aza_degs$feature %in% base::intersect(rownames(aza_top25percent_varPart), topVar_genes),]

aza_up_degs_topVar_varPart25 <- aza_degs_topVar_varPart25 %>% 
  filter(coef > 0) %>% 
  pull(feature)

aza_down_degs_topVar_varPart25 <- aza_degs_topVar_varPart25 %>% 
  filter(coef < 0) %>% 
  pull(feature)

# Pred
pred_degs_lfc0.5 <- pred_degs %>% 
  filter(abs(coef) > 0.5)

pred_degs_topVar <- pred_degs_lfc0.5[pred_degs_lfc0.5$feature %in% base::intersect(pred_degs_lfc0.5$feature, rownames(topVar)),]

pred_up_degs_topVar <- pred_degs_topVar %>% 
  filter(coef > 0) %>% 
  pull(feature)

pred_down_degs_topVar <- pred_degs_topVar %>% 
  filter(coef < 0) %>% 
  pull(feature)

# Anti-TNF
## Shared with topVar
antiTNF_degs_topVar <- antiTNF_degs[antiTNF_degs$feature %in% base::intersect(antiTNF_degs$feature, rownames(topVar)),]
antiTNF_up_degs_topVar <- antiTNF_degs_topVar %>% 
  filter(coef > 0) %>% 
  pull(feature)

antiTNF_down_degs_topVar <- antiTNF_degs_topVar %>% 
  filter(coef < 0) %>% 
  pull(feature)

## Unique DEGs (not shared with topVar)
antiTNF_unique_degs <- antiTNF_degs[antiTNF_degs$feature %in% setdiff(antiTNF_degs$feature, rownames(topVar)),]
antiTNF_up_unique_degs <- antiTNF_unique_degs %>% 
  filter(coef > 0) %>% 
  pull(feature)

antiTNF_down_unique_degs <- antiTNF_unique_degs %>% 
  filter(coef < 0) %>% 
  pull(feature)

# Computing base mean expression of genes ----
basemean_go <- function(coldata){
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
  res_deseq <- subset(res_deseq, res_deseq$baseMean > 0) # Filtering out genes that were never detected
  
  return(res_deseq)
}

aza_res_deseq <- basemean_go(aza_coldata)
pred_res_deseq <- basemean_go(pred_coldata)
antiTNF_res_deseq <- basemean_go(antiTNF_coldata)

# GO ----
GO_analysis <- function(degs_id, res_deseq, regulation, top_go_path, go_table_path){
  # Universe genes
  all_genes_id <- rownames(res_deseq)
  
  overallBaseMean <- as.matrix(res_deseq[, "baseMean", drop = FALSE])
  sig_idx <- match(degs_id, rownames(overallBaseMean))
  
  backG <- c()
  for (i in sig_idx) {
    ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
    backG <- c(backG, ind)
  }
  backG <- unique(backG)
  backG <- rownames(overallBaseMean)[backG]
  backG <- setdiff(backG, degs_id)
  
  multidensity(list(all = log2(res_deseq[, "baseMean"]), 
                    foreground = log2(res_deseq[degs_id, "baseMean"]), 
                    background = log2(res_deseq[backG, "baseMean"])), 
               xlab = "log2 mean normalized counts",
               main = paste("Matching", regulation, "for enrichment analysis", sep = ' ')) 
  
  # Selecting genes 
  inUniverse <- all_genes_id %in% c(degs_id, backG)
  inSelection <- all_genes_id %in% degs_id
  alg <- factor(as.integer(inSelection[inUniverse]))
  names(alg) <- all_genes_id[inUniverse]
  
  # Preparing data
  tgd <- new("topGOdata", ontology = "BP", allGenes = alg, nodeSize = 5, annot = annFUN.org,
             mapping = "org.Hs.eg.db", ID = "ensembl") 
  
  # Running tests
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
             fn.prefix = top_go_path,
             useInfo = "all", pdfSW = TRUE) 
  
  # Finding the genes associated with the enriched GO terms
  for (i in 1:nrow(tab)) {
    go_id <- as.vector(tab[i, 1])
    genes_in_cat <- intersect(unlist(genesInTerm(tgd, go_id)), degs_id)
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
  
  tab <- tab %>% 
    mutate(geneRatio = tab[,4]/tab[,3]) %>% 
    relocate(geneRatio, .after = Expected) %>% 
    mutate(across(Fisher.elim, as.numeric)) %>% 
    filter(Fisher.elim < 0.05)
  
  tab$Term = with(tab, reorder(Term, Fisher.elim, decreasing = TRUE))
  
  write.table(tab, go_table_path, quote = FALSE, sep = '\t')
  
  return(tab)
}

# Aza
go_aza_up <- GO_analysis(aza_up_degs_topVar_varPart25, aza_res_deseq, "aza upregulated genes", 
                         "Output_files/GO/aza/aza_elim_NoBackground_up",
                         "Output_files/GO/aza/aza_up_go.txt")
go_aza_down <- GO_analysis(aza_down_degs_topVar_varPart25, aza_res_deseq, "aza downregulated genes", 
                           "Output_files/GO/aza/aza_elim_NoBackground_down",
                           "Output_files/GO/aza/aza_down_go.txt")

# Pred
go_pred_up<- GO_analysis(pred_up_degs_topVar, pred_res_deseq, "pred upregulated genes", 
                         "Output_files/GO/pred/pred_elim_NoBackground_up",
                         "Output_files/GO/pred/pred_up_go.txt")
go_pred_down <- GO_analysis(pred_down_degs_topVar, pred_res_deseq, "pred downregulated genes", 
                            "Output_files/GO/pred/pred_elim_NoBackground_down",
                            "Output_files/GO/pred/pred_down_go.txt")

# Anti-TNF
go_antiTNF_up_shared <- GO_analysis(antiTNF_up_degs_topVar, antiTNF_res_deseq, "antiTNF shared upregulated genes", 
                                    "Output_files/GO/antiTNF/antiTNF_shared_elim_NoBackground_up",
                                    "Output_files/GO/antiTNF/antiTNF_shared_up_go.txt")
go_antiTNF_down_shared <- GO_analysis(antiTNF_down_degs_topVar, antiTNF_res_deseq, "antiTNF shared downregulated genes", 
                                      "Output_files/GO/antiTNF/antiTNF_shared_elim_NoBackground_down",
                                      "Output_files/GO/antiTNF/antiTNF_shared_down_go.txt")

go_antiTNF_up_unique <- GO_analysis(antiTNF_up_unique_degs, antiTNF_res_deseq, "antiTNF unique upregulated genes", 
                                    "Output_files/GO/antiTNF/antiTNF_unique_elim_NoBackground_up",
                                    "Output_files/GO/antiTNF/antiTNF_unique_up_go.txt")
go_antiTNF_down_unique <- GO_analysis(antiTNF_down_unique_degs, antiTNF_res_deseq, "antiTNF unique upregulated genes", 
                                      "Output_files/GO/antiTNF/antiTNF_unique_elim_NoBackground_up",
                                      "Output_files/GO/antiTNF/antiTNF_unique_down_go.txt")

# Dot plots ----
# Aza
go_terms_aza_up <- c("angiogenesis", "inorganic ion homeostasis", "tissue development", "regulation of anatomical structure morphogenesis",
                  "phosphorylation", "monoatomic cation homeostasis")
go_terms_aza_down <- c("natural killer cell mediated cytotoxicity", "positive regulation of lymphocyte mediated immunity",
               "G protein-coupled receptor signaling pathway", "positive regulation of cytokine production involved in immune response",
               "cell surface receptor signaling pathway", "positive regulation of programmed cell death",
               "positive regulation of innate immune response", "immune response", "response to type II interferon",
               "cell killing", "cellular response to tumor necrosis factor")

plot_aza_up <- go_aza_up %>% 
  mutate(Dotplot = ifelse(Term %in% go_terms_aza_up, "YES", "NO")) %>% 
  filter(Dotplot == "YES") %>% 
  mutate(Direction = "Upregulated")

plot_aza_down <- go_aza_down %>% 
  mutate(Dotplot = ifelse(Term %in% go_terms_aza_down, "YES", "NO")) %>% 
  filter(Dotplot == "YES") %>% 
  mutate(Direction = "Downregulated")

all_terms_aza <- rbind(plot_aza_up, plot_aza_down)
all_terms_aza  <- all_terms_aza %>%
  group_by(Analysis) %>% 
  mutate(Term = fct_reorder2(Term, geneRatio,-geneRatio)) %>%
  ungroup()

pdf("Output_files/GO/aza/aza_dotplot_GO.pdf", width = 12, height = 8)
dot_plot_aza <- ggplot(all_terms_aza, aes(x = Analysis, y= Term, size = geneRatio, color = Fisher.elim)) +
  geom_point() +
  scale_colour_gradient(low = "#990000", high = "#FF9999") +
  xlab('') + ylab('') +
  labs(
    title = "GO Azathioprine x No Syst",
    color = "Adj p-value",
    size = "Gene ratio") +
  theme_bw() +
  theme(
    axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
    axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1),
    axis.title = element_text(size = 20),
    legend.position = "bottom"
  )
dot_plot_aza
dev.off()

# Pred
pred_terms <- c("regulation of toll-like receptor signaling pathway", "neutrophil extravasation", 
                "positive regulation of NF-kappaB transcription factor activity", "response to cAMP",
                "inflammatory response", "regulation of macrophage activation", "positive regulation of cytokine production",
                "regulation of lipid storage", "acute inflammatory response", "positive regulation of inflammatory response",
                "cell surface receptor signaling pathway", "defense response to bacterium", "leukocyte proliferation", 
                "leukocyte differentiation", "cell population proliferation", "T cell activation", "B cell receptor signaling pathway",
                "adaptive immune response")
                
plot_pred_up <- go_pred_up %>% 
  mutate(plot = ifelse(Term %in% pred_terms, "Yes", "No")) %>% 
  filter(plot == "Yes") %>% 
  mutate(Direction = "Upregulated")

plot_pred_down <- go_pred_down %>% 
  mutate(plot = ifelse(Term %in% pred_terms, "Yes", "No")) %>% 
  filter(plot == "Yes") %>% 
  mutate(Direction = "Downregulated")

all_terms_pred <- rbind(plot_pred_up, plot_pred_down)
all_terms_pred$Term <- Term(as.vector(all_terms_pred$GO.ID))
all_terms_pred$Term <- factor(all_terms_pred$Term, levels = rev(unique(all_terms_pred$Term))) 

pdf("Output_files/GO/pred/pred_dotplot_GO.pdf", width = 12, height = 8)
dot_plot_pred <- ggplot(all_terms_pred, aes(x = Direction, y= Term, size = geneRatio, color = Fisher.elim)) +
  geom_point() +
  scale_colour_gradient(low = "#990000", high = "#FF9999") +
  xlab('') + ylab('') +
  labs(
    title = "GO BP Pred x No Syst",
    color = "Adj p-value",
    size = "Gene ratio") +
  theme_bw() + 
  theme(
    axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
    axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1),
    axis.title = element_text(size = 20),
    legend.position = "bottom"
  ) +
    scale_y_discrete(limits = c("adaptive immune response",  "B cell receptor signaling pathway",
                                "T cell activation", "cell population proliferation", "leukocyte differentiation", "leukocyte proliferation",
   
                                "defense response to bacterium","cell surface receptor signaling pathway",
   
                                "positive regulation of inflammatory response", "acute inflammatory response", "regulation of lipid storage",
                                "positive regulation of cytokine production", "regulation of macrophage activation", "inflammatory response",
                                "response to cAMP", "regulation of macrophage activation", "positive regulation of NF-kappaB transcription factor activity",
                                "neutrophil extravasation", "regulation of toll-like receptor signaling pathway", "cell population proliferation"))
dot_plot_pred
dev.off()

# Anti-TNF
go_antiTNF_up_shared$Direction <- "Upregulated"
go_antiTNF_up_shared$Analysis <- "Shared"

go_antiTNF_down_shared$Direction <- "Downregulated"
go_antiTNF_down_shared$Analysis <- "Shared"

go_antiTNF_up_unique$Direction <- "Upregulated" 
go_antiTNF_up_unique$Analysis <- "Unique"

go_antiTNF_down_unique$Direction <- "Downregulated" 
go_antiTNF_down_unique$Analysis <- "Unique"

all_terms_tnf <- rbind(go_antiTNF_up_shared[c(1,2,4),], go_antiTNF_down_shared[1:6,], go_antiTNF_down_unique, go_antiTNF_up_unique)
all_terms_tnf  <- all_terms_tnf %>%
  group_by(Analysis) %>% 
  mutate(Term = fct_reorder2(Term, geneRatio,-geneRatio)) %>%
  ungroup()

pdf("Output_files/GO/antiTNF/antiTNF_dotplot_GO_unique_shared.pdf", width = 12, height = 7)
dot_plot_antiTNF <- ggplot(all_terms_tnf, aes(x = Direction, y= Term, size = geneRatio, color = Fisher.elim)) +
  geom_point() +
  scale_colour_gradient(low = "#990000", high = "#FF9999") +
  scale_size_continuous(range=c(3,8)) +
  xlab('') + ylab('') +
  labs(
    title = "GO Anti-TNF x No Bio",
    color = "Adj p-value",
    size = "Gene ratio") +
  theme_bw() + 
  theme(
    axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
    axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1),
    axis.title = element_text(size = 20),
    legend.position = 'bottom'
  ) + facet_wrap(~Analysis)
dot_plot_antiTNF
dev.off()
