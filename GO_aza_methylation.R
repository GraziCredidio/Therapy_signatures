# EZE cohort: Therapy Signatures
  # DNA methylation-linked genes shared with TVGs: Gene ontology analysis
  # Azathioprine x No Systemic Therapy 
  # Author: Graziella Credidio

rm(list = ls())

folder <- "Output_files/Methylation/aza/GO"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading packages ----
library(topGO)
library(tidyverse)
library(kableExtra)

# Loading data ----
sig_gene_meth_site_correlation <- read.table("Output_files/Methylation/aza/aza_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')
transcriptome_data <- read.table("Cleaned_tables/models/aza/normalized_counts_aza.txt", sep = "\t")
degs <- read.table("Output_files/Maaslin2/significant_results/maaslin2_significant_results_aza_vs_noSyst.txt", sep = '\t')
topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVarGenes_names <- rownames(topVarGenes)

# Split into down and upregulated DNAm linked DEGs in common with TVGs ----
up_sig_genes <- degs %>% 
  filter(coef > 0) %>% 
  pull(feature)
up_dnam_degs <- sig_gene_meth_site_correlation[sig_gene_meth_site_correlation$Gene %in% up_sig_genes,]
upgenes <- unique(up_dnam_degs$Gene)
upgenes_topVar <- intersect(upgenes, topVarGenes_names)
write.table(up_dnam_degs, file = "Output_files/Methylation/aza/aza_upregulated_DNAm_DEG.txt", quote = FALSE, sep = '\t')

down_sig_genes <- degs %>% 
  filter(coef < 0) %>% 
  pull(feature)
down_dnam_degs <- sig_gene_meth_site_correlation[sig_gene_meth_site_correlation$Gene %in% down_sig_genes,]
downgenes <- unique(down_dnam_degs$Gene)
downgenes_topvar <- intersect(downgenes, topVarGenes_names)
write.table(down_dnam_degs, file = "Output_files/Methylation/aza/aza_downregulated_DNAm_DEG.txt", quote = FALSE, sep = '\t')

# GO analysis for upregulated and downregulated genes ----
gene_set <- list(upgenes_topVar, downgenes_topvar)
direction <- 1 #1 = upregulated, 2= downregulated

# Universe gene set
all_genes <- rownames(transcriptome_data)

for (g in 1:length(gene_set)){
  #Define background set and test set
  inUniverse <- all_genes %in% all_genes
  inSelection <- all_genes %in% as.vector(gene_set[[g]])
  alg <- factor( as.integer( inSelection[inUniverse] ) )
  names(alg) <- all_genes[inUniverse]
  
  #Run topGO analysis
  tgd <- new("topGOdata", ontology="BP", allGenes = alg, nodeSize=5,
             annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )
  
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  
  go_results <- GenTable(tgd, Fisher.elim = resultTopGO.elim, 
                         Fisher.classic = resultTopGO.classic,
                         orderBy = "Fisher.elim", numChar=1000, topNodes = length(nodes(graph(tgd))))
  
  #Extract significant results and add genes contributing to each GO term
  go_results <- subset(go_results, go_results$Fisher.elim < 0.05)
  for(i in 1:nrow(go_results)){
    go_id <- as.vector(go_results[i,1])
    alleges <- get(go_id, org.Hs.egGO2ALLEGS)
    genes <- unlist(mget(alleges, org.Hs.egENSEMBL))
    genes_in_cat <- intersect(as.vector(genes), as.vector(gene_set[[g]])) 
    gene_sym_in_cat <- as.vector(unlist(mget(unlist(mget(genes_in_cat, org.Hs.egENSEMBL2EG)), org.Hs.egSYMBOL)))
    gene_sym_in_cat_str <- ""
    
    if(length(genes_in_cat) > 0){
      for(j in 1:length(gene_sym_in_cat)){
        gene_sym_in_cat_str <- paste(gene_sym_in_cat_str, gene_sym_in_cat[j], sep = ',')
      }
    }
    if (gene_sym_in_cat_str == ""){
      gene_sym_in_cat_str <- NA
    }
    go_results$Genes[i] <- gene_sym_in_cat_str
  }
  
  # Calculating gene ratio and adding it as a column of go_results
  go_results <- go_results %>% 
    mutate(geneRatio = go_results[,4]/go_results[,3]) %>% 
    relocate(geneRatio, .after = Expected)
  
  # Saving files
  if (direction > 1){
    write.table(go_results, "Output_files/Methylation/aza/GO/GO_aza_TVG_DNAm_down_DEGs.txt", sep = '\t', quote = FALSE)
  } else{
    write.table(go_results, "Output_files/Methylation/aza/GO/GO_aza_TVG_DNAm_up_DEGs.txt",sep = '\t', quote = FALSE)
  }
  
  direction <- direction + 1
}

# Plotting GO analysis ----
# Reading GO for upregulated and downregulated DNAm-DEGs 
up_go_result <- read.table("Output_files/Methylation/aza/GO/GO_aza_TVG_DNAm_up_DEGs.txt", sep = '\t')
down_go_result <- read.table("Output_files/Methylation/aza/GO/GO_aza_TVG_DNAm_down_DEGs.txt", sep = '\t')


# Selected GO enriched terms 
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


# Data preprocessing
plot_up_go_result <- up_go_result %>% 
  mutate(plot = ifelse(Term %in% terms, "Yes", "No")) %>% 
  filter(plot == "Yes") %>% 
  mutate(direction = "Upregulated")

plot_down_go_result <- down_go_result %>% 
  mutate(plot = ifelse(Term %in% terms, "Yes", "No")) %>% 
  filter(plot == "Yes") %>% 
  mutate(direction = "Downregulated")

selected_GO <- rbind(plot_up_go_result, plot_down_go_result)
selected_GO$Term <- Term(as.vector(selected_GO$GO.ID))
selected_GO$Term <- factor(selected_GO$Term, levels = rev(selected_GO$Term)) # fixes order
selected_GO$Genes <- sub('.', '', selected_GO$Genes)

# Table with information
selected_GO %>%
  dplyr::select(GO.ID, Term, Annotated, Significant, geneRatio, Fisher.elim, Genes, direction) %>% 
  kbl() %>%
  kable_styling()

# Plot (Figure 19A)
pdf("Output_files/Methylation/aza/GO/GO_dotplot_BP_DNAm_TVG_aza.pdf", width = 10, height = 8)
dot_plot <- ggplot(selected_GO, aes(x = direction, y = Term, size = geneRatio, color = Fisher.elim)) +
  geom_point() +
  scale_colour_gradient(low = "#990000", high = "#FF9999") +
  xlab('') + ylab('') +
  labs(
    title = "GO: DNam-DEGs shared with TVG",
    subtitle = "Azathioprine x No Systemic Therapy",
    color = "Adj p-value",
    size = "Gene ratio") +
  theme_bw() + 
  theme(
    axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
    axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1),
    axis.title = element_text(size = 20),
    legend.position="bottom"
  )
dot_plot
dev.off()
