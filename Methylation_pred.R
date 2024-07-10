# EZE Cohort - therapy signatures
# Methylation data integration with transcriptome
# Pred x No syst - remitters + crp significant genes

rm(list = ls())


library(topGO)
library(org.Hs.eg.db)
library(GO.db)
library(ComplexHeatmap)
library(RColorBrewer)
library(plyr)
library(tibble)
library(dplyr)
library(ggplot2)



#Correlation analysis between transcription and methylation--------
#Read and filter methylation data for sites
methylation_data <- read.csv("Methylation_data/m_values_hg38_round.txt", header = TRUE, row.names = NULL, 
                             sep = '\t')

which(is.na(methylation_data$row.names)) #10359 = NA
methylation_data <- methylation_data[-10359,]

meth_sample_info <- read.csv("Methylation_data/DZHK_EZE_epicarray_metadata_with_remission_status_090223.csv",
                             header = TRUE, sep = ',')


methylation_data_filtered <- methylation_data[, colnames(methylation_data) %in% paste('X', as.character(meth_sample_info$Sample_ID), sep = '')]
rownames(methylation_data_filtered) <- methylation_data$row.names
names(methylation_data_filtered) <- gsub(pattern= "X", replacement = "", x=names(methylation_data_filtered))

#Matching order of sample info and methylation data filtered
meth_sample_info_filterd <- meth_sample_info[meth_sample_info$Sample_ID %in% colnames(methylation_data_filtered),]
idx <- match(colnames(methylation_data_filtered), meth_sample_info_filterd$Sample_ID) 
ord.meth_sample_info_filtered <- meth_sample_info_filterd[idx, ]

all(colnames(methylation_data_filtered) %in% ord.meth_sample_info_filtered$Sample_ID)
all(colnames(methylation_data_filtered) == ord.meth_sample_info_filtered$Sample_ID)

# replacing sample_id by Study ID
colnames(methylation_data_filtered) <- ord.meth_sample_info_filtered$Study.ID


#Read and filter transcriptome data for genes
transcriptome_data <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_normalized_pred_R_crp_04.05.txt", sep = "\t")
transcriptome_info_file <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_pred_R_crp_04.05.txt", sep = "\t") #coldata

transcriptome_data <- transcriptome_data[, as.character(transcriptome_info_file$sample_id)] #normalized counts filtered by pred patients

#Matching order of sample info and transcriptome sig
idx2 <- match(colnames(transcriptome_data), transcriptome_info_file$sample_id) 
ord.transcriptome_info_file <- transcriptome_info_file[idx2, ]

all(colnames(transcriptome_data) %in% ord.transcriptome_info_file$sample_id)
all(colnames(transcriptome_data) == ord.transcriptome_info_file$sample_id)

colnames(transcriptome_data) <- ord.transcriptome_info_file$study_id

common_samples <- intersect(colnames(transcriptome_data), colnames(methylation_data_filtered))
transcriptome_data <- transcriptome_data[, common_samples]
methylation_data_filtered <- methylation_data_filtered[, common_samples]

#Read list of differentially expressed genes 
deg_gene_list <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt",
                            sep = '\t')
deg_gene_list <- deg_gene_list$feature

#Read and filter gene and linked methylation sites info, remove sites on sex chromosomes
gene_meth_sites <- read.csv("Methylation_data/gene_meth_sites_5000bp.txt", sep = '\t')
gene_meth_sites <- gene_meth_sites[,-9] #removing column X

gene_meth_sites <- subset(gene_meth_sites, gene_meth_sites$Chr != 'chrX' & gene_meth_sites$Chr != 'chrY' & 
                            gene_meth_sites$Chr != 'chrM')
rownames(gene_meth_sites) <- gene_meth_sites$Gene_id

deg_gene_list <- intersect(deg_gene_list, rownames(gene_meth_sites))

deg_gene_meth_sites <- gene_meth_sites[as.character(deg_gene_list), ]
deg_gene_meth_sites <- subset(deg_gene_meth_sites, deg_gene_meth_sites$no_of_meth_sites > 0)


#Calculate correlation coefficient between the gene expression of each gene and the methylation intensity of each of its link methylated sites
#Calculate FDR using a permutation approach
deg_gene_meth_sites$meth_sites <- as.character(deg_gene_meth_sites$meth_sites)
deg_gene_meth_sites$no_of_meth_sites <- as.numeric(deg_gene_meth_sites$no_of_meth_sites)
gene_meth_site_correlation <- matrix(nrow = sum(deg_gene_meth_sites$no_of_meth_sites), ncol = 5)
rownum = 1


for (i in 1:nrow(deg_gene_meth_sites)) {
  meth_sites <- strsplit(deg_gene_meth_sites[i,8], split = ';')
  print(meth_sites)
  for (j in 1:length(meth_sites[[1]])) {
    site_id_info <- strsplit(meth_sites[[1]][j], split = ':')
    site_id <- site_id_info[[1]][1]
    distance <- site_id_info[[1]][3]
    print(site_id)
    
    if(site_id %in% rownames(methylation_data_filtered)){ #exclude sites that are not in methylation data
      corr_data_frame <- data.frame(meth=t(methylation_data_filtered[site_id,]), 
                                    expr=t(transcriptome_data[as.character(deg_gene_meth_sites[i,1]),]))
    }
    
    corr_data_frame <- corr_data_frame[complete.cases(corr_data_frame),] 
    rho <- cor(corr_data_frame[,1], corr_data_frame[,2], method="spearman")
    fdr <- 0
    for (k in 1:1000) {
      x <- sample(corr_data_frame[,1])
      y <- sample(corr_data_frame[,2])
      rand_rho <- cor(x, y, method="spearman")
      if(abs(rand_rho) >= abs(rho)){
        fdr <- fdr + 1
      }
    }
    fdr <- fdr/1000
    gene_meth_site_correlation[rownum, ] <- c(as.character(deg_gene_meth_sites[i,1]), site_id, distance, rho, fdr)
    rownum <- rownum + 1
  }
}

gene_meth_site_correlation <- as.data.frame(gene_meth_site_correlation)
colnames(gene_meth_site_correlation) <- c("Gene", "Site", "Distance_from_TSS", "Rho", "FDR")
gene_meth_site_correlation$FDR <- as.numeric(as.character(gene_meth_site_correlation$FDR)) 

write.table(gene_meth_site_correlation, file = "Output_files/Methylation/pred/pred_DEG_meth_site_correlation_5000bp_1000rep.txt", quote = FALSE,
            sep = '\t')


#Extract significant DEGs  ----
#gene_meth_site_correlation <- read.table("Output_files/Methylation/pred/pred_DEG_meth_site_correlation_5000bp_1000rep.txt", sep = '\t')

sig_gene_meth_site_correlation <- gene_meth_site_correlation %>% #27045
  filter(FDR < 0.05)

# Add gene names
ensg2gene <- read.table("Cleaned_tables/ensg2gene_EZE.txt", header = TRUE, sep = ",")
unique_ensg2gene <- subset(ensg2gene, duplicated(ensg2gene$ensembl_gene_id) == FALSE)
rownames(unique_ensg2gene) <- unique_ensg2gene$ensembl_gene_id

# save
sig_gene_meth_site_correlation$Gene_name <- unique_ensg2gene[as.character(sig_gene_meth_site_correlation$Gene), "hgnc_symbol"]
write.table(sig_gene_meth_site_correlation, file = "Output_files/Methylation/pred/pred_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", quote = FALSE,
            sep = '\t')


# split into up and downregulated significant genes from correlation
sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt",
                        sep = '\t')

up_sig_genes <- sig_genes %>% 
  filter(coef > 0) %>% 
  pull(feature)

up_sig_gene_meth_site_correlation <- sig_gene_meth_site_correlation[sig_gene_meth_site_correlation$Gene %in% up_sig_genes,]
write.table(up_sig_gene_meth_site_correlation, file = "Output_files/Methylation/pred/pred_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", quote = FALSE,
            sep = '\t')


down_sig_genes <- sig_genes %>% 
  filter(coef < 0) %>% 
  pull(feature)

down_sig_gene_meth_site_correlation <- sig_gene_meth_site_correlation[sig_gene_meth_site_correlation$Gene %in% down_sig_genes,]
write.table(down_sig_gene_meth_site_correlation, file = "Output_files/Methylation/pred/pred_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", quote = FALSE,
            sep = '\t')


#GO term enrichment analysis: all correlated genes --------
#loading files
sig_gene_meth_site_correlation <- read.table("Output_files/Methylation/pred/pred_significant_DEG_meth_site_correlation_5000bp_1000rep.txt",
                                             sep = '\t')
up_sig_gene_meth_site_correlation <- read.table("Output_files/Methylation/pred/pred_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt",
                                                sep = '\t')
down_sig_gene_meth_site_correlation <- read.table("Output_files/Methylation/pred/pred_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", 
                                                  sep = '\t')
transcriptome_data <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_normalized_pred_R_crp_04.05.txt", 
                                 sep = "\t")

# setting genes to GO analysis 
upgenes <- unique(up_sig_gene_meth_site_correlation$Gene)
downgenes <- unique(down_sig_gene_meth_site_correlation$Gene)
all_genes <- rownames(transcriptome_data)

#Perform GO term analysis for upregulated and downregulated genes
direction <- c("up")
gene_set <- upgenes

direction <- c("down")
gene_set <- downgenes

for (i in 1:length(direction)) {
  #Define background set and test set
  inUniverse <- all_genes %in% all_genes
  inSelection <- all_genes %in% as.vector(gene_set)
  alg <- factor( as.integer( inSelection[inUniverse] ) )
  names(alg) <- all_genes[inUniverse]
  
  #Run topGO analysis
  tgd <- new( "topGOdata", ontology="BP", allGenes = alg, nodeSize=5,
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
    genes_in_cat <- intersect(as.vector(genes), as.vector(gene_set)) #changed
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
  
  go_results <- go_results %>% 
    mutate(geneRatio = go_results[,4]/go_results[,3]) %>% 
    relocate(geneRatio, .after = Expected)
  
  
  write.table(go_results, paste("Output_files/Methylation/pred/GO/pred_correlated_", direction, "_DEGs_GO.txt", sep = ''), 
              sep = '\t', quote = FALSE)
  
}


up_go_result <- read.csv("Output_files/Methylation/pred/GO/pred_correlated_up_DEGs_GO.txt", sep = '\t')
down_go_result <- read.csv("Output_files/Methylation/pred/GO/pred_correlated_down_DEGs_GO.txt", sep = '\t')

#Plot selected GO terms enriched in DEGs 
terms <- c("protein autophosphorylation", "innate immune response in mucosa", "NIK/NF-kappaB signaling", "negative regulation of cell population proliferation",
           "positive regulation of cell population proliferation", "cell population proliferation", "toll-like receptor signaling pathway",
           "leukocyte chemotaxis", "innate immune response in mucosa", "G protein-coupled receptor signaling pathway", "response to cAMP", 
           "negative regulation of inflammatory response","positive regulation of cytokine production involved in inflammatory response",
           
           
           "B cell proliferation", "pyrimidine nucleoside diphosphate biosynthetic process", "peptide biosynthetic process", "ribosomal large subunit biogenesis",
           "ribosomal small subunit biogenesis", "rRNA modification")


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

pdf("Output_files/Methylation/pred/GO/dotplot_GO_BP_pred_R_crp_data_integration.pdf", width = 9, height = 6)
dot_plot <- ggplot(selected_GO, aes(x = direction, y= Term, size = geneRatio, color = Fisher.elim)) +
  geom_point() +
  scale_colour_gradient(low = "#990000", high = "#FF9999") +
  #scale_size_continuous(range=c(3,10)) +
  xlab('') + ylab('') +
  labs(
    title = "GO BP Pred x No Syst",
    subtitle = "Transcriptome-methylome integration",
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







#GO term enrichment analysis: correlated genes + topvar --------
#loading files
topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVarGenes_names <- rownames(topVarGenes)

sig_gene_meth_site_correlation <- read.table("Output_files/Methylation/pred/pred_significant_DEG_meth_site_correlation_5000bp_1000rep.txt",
                                             sep = '\t')
up_sig_gene_meth_site_correlation <- read.table("Output_files/Methylation/pred/pred_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt",
                                                sep = '\t')
down_sig_gene_meth_site_correlation <- read.table("Output_files/Methylation/pred/pred_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", 
                                                  sep = '\t')
transcriptome_data <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_normalized_pred_R_crp_04.05.txt", 
                                 sep = "\t")

# setting genes to GO analysis 
upgenes <- unique(up_sig_gene_meth_site_correlation$Gene)
upgenes_topVar <- intersect(upgenes, topVarGenes_names)

downgenes <- unique(down_sig_gene_meth_site_correlation$Gene)
downgenes_topvar <- intersect(downgenes, topVarGenes_names)

all_genes <- rownames(transcriptome_data)

#Perform GO term analysis for upregulated and downregulated genes
direction <- c("up")
gene_set <- upgenes_topVar

direction <- c("down")
gene_set <- downgenes_topvar

for (i in 1:length(direction)) {
  #Define background set and test set
  inUniverse <- all_genes %in% all_genes
  inSelection <- all_genes %in% as.vector(gene_set)
  alg <- factor( as.integer( inSelection[inUniverse] ) )
  names(alg) <- all_genes[inUniverse]
  
  #Run topGO analysis
  tgd <- new( "topGOdata", ontology="BP", allGenes = alg, nodeSize=5,
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
    genes_in_cat <- intersect(as.vector(genes), as.vector(gene_set)) #changed
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
  
  go_results <- go_results %>% 
    mutate(geneRatio = go_results[,4]/go_results[,3]) %>% 
    relocate(geneRatio, .after = Expected)
  
  
  write.table(go_results, paste("Output_files/Methylation/pred/GO/correlated_topvar/pred_correlated_topVar_", direction, "_DEGs_GO.txt", sep = ''), 
              sep = '\t', quote = FALSE)
  
}


up_go_result <- read.csv("Output_files/Methylation/pred/GO/correlated_topvar/pred_correlated_topVar_up_DEGs_GO.txt", sep = '\t')
down_go_result <- read.csv("Output_files/Methylation/pred/GO/correlated_topvar/pred_correlated_topVar_down_DEGs_GO.txt", sep = '\t')

#Plot selected GO terms enriched in DEGs 
terms <- c("inflammatory response", "response to cAMP", "positive regulation of epidermal growth factor receptor signaling pathway",
           "acute inflammatory response", "macrophage activation", "detection of molecule of bacterial origin",
           "response to type II interferon", "neutrophil activation", "viral entry into host cell", "regulation of cholesterol storage",
           "regulation of toll-like receptor signaling pathway",

           "B cell receptor signaling pathway", "B cell proliferation", "B cell differentiation", "regulation of natural killer cell mediated cytotoxicity",
           "regulation of humoral immune response mediated by circulating immunoglobulin", "purine deoxyribonucleotide biosynthetic process")



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

pdf("Output_files/Methylation/pred/GO/correlated_topvar/dotplot_GO_BP_pred_R_crp_data_integration_topVar.pdf", width = 9, height = 6)
dot_plot <- ggplot(selected_GO, aes(x = direction, y= Term, size = geneRatio, color = Fisher.elim)) +
  geom_point() +
  scale_colour_gradient(low = "#990000", high = "#FF9999") +
  #scale_size_continuous(range=c(3,10)) +
  xlab('') + ylab('') +
  labs(
    title = "GO BP Pred x No Syst",
    subtitle = "Transcriptome-methylome integration + topVar",
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







#GO term enrichment analysis: correlated genes (lfc 0.5) + topvar --------
#loading files
topVarGenes <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")
topVarGenes_names <- rownames(topVarGenes)

sig_genes <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_pred_vs_noSyst_R_crp.txt", sep = "\t")
sig_genes_lfc0.5 <- sig_genes %>% 
  filter(abs(coef) > 0.5)

# sig_gene_meth_site_correlation <- read.table("Output_files/Methylation/pred/pred_significant_DEG_meth_site_correlation_5000bp_1000rep.txt",
#                                              sep = '\t')

up_sig_gene_meth_site_correlation <- read.table("Output_files/Methylation/pred/pred_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt",
                                                sep = '\t')
up_sig_lfc0.5_correlation <- intersect(up_sig_gene_meth_site_correlation$Gene,sig_genes_lfc0.5$feature )
up_sig_gene_meth_site_correlation_lfc0.5 <- up_sig_gene_meth_site_correlation[up_sig_gene_meth_site_correlation$Gene %in% up_sig_lfc0.5_correlation,]


down_sig_gene_meth_site_correlation <- read.table("Output_files/Methylation/pred/pred_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", 
                                                  sep = '\t')
down_sig_lfc0.5_correlation <- intersect(down_sig_gene_meth_site_correlation$Gene,sig_genes_lfc0.5$feature )
down_sig_gene_meth_site_correlation_lfc0.5 <- down_sig_gene_meth_site_correlation[down_sig_gene_meth_site_correlation$Gene %in% down_sig_lfc0.5_correlation,]


transcriptome_data <- read.table("Output_files/DESeq2/Maaslin2/Remission_crp/DESeq2_normalized_pred_R_crp_04.05.txt", 
                                 sep = "\t")

# setting genes to GO analysis 
upgenes <- unique(up_sig_gene_meth_site_correlation_lfc0.5$Gene)
upgenes_topVar <- intersect(upgenes, topVarGenes_names) #292

downgenes <- unique(down_sig_gene_meth_site_correlation$Gene)
downgenes_topvar <- intersect(downgenes, topVarGenes_names) #114

all_genes <- rownames(transcriptome_data)

#Perform GO term analysis for upregulated and downregulated genes
direction <- c("up")
gene_set <- upgenes_topVar

direction <- c("down")
gene_set <- downgenes_topvar

for (i in 1:length(direction)) {
  #Define background set and test set
  inUniverse <- all_genes %in% all_genes
  inSelection <- all_genes %in% as.vector(gene_set)
  alg <- factor( as.integer( inSelection[inUniverse] ) )
  names(alg) <- all_genes[inUniverse]
  
  #Run topGO analysis
  tgd <- new( "topGOdata", ontology="BP", allGenes = alg, nodeSize=5,
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
    genes_in_cat <- intersect(as.vector(genes), as.vector(gene_set)) #changed
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
  
  go_results <- go_results %>% 
    mutate(geneRatio = go_results[,4]/go_results[,3]) %>% 
    relocate(geneRatio, .after = Expected)
  
  
  write.table(go_results, paste("Output_files/Methylation/pred/GO/correlated_topvar/lfc0.5/pred_lfc0.5_correlated_topVar_", direction, "_DEGs_GO.txt", sep = ''), 
              sep = '\t', quote = FALSE)
  
}

# comparing terms of DMLGs lfc0.5 + topvar with DEG lfc 0.5 + topvar
up_go_result <- read.csv("Output_files/Methylation/pred/GO/correlated_topvar/lfc0.5/pred_lfc0.5_correlated_topVar_up_DEGs_GO.txt", sep = '\t')
down_go_result <- read.csv("Output_files/Methylation/pred/GO/correlated_topvar/lfc0.5/pred_lfc0.5_correlated_topVar_down_DEGs_GO.txt", sep = '\t')

#Plot selected GO terms enriched in DEGs 
#same terms as GO with DEGS lfc 0.5 + topvar
same_terms <- c("positive regulation of inflammatory response", "defense response to bacterium", "acute inflammatory response", 
           "positive regulation of cytokine production", "inflammatory response",
           "regulation of macrophage activation", "positive regulation of NF-kappaB transcription factor activity",
           "neutrophil extravasation", "regulation of toll-like receptor signaling pathway", "regulation of lipid storage",
            "cell surface receptor signaling pathway", "B cell receptor signaling pathway")


#new_terms + same terms
terms <- c("positive regulation of inflammatory response", "defense response to bacterium", "acute inflammatory response", 
                   "positive regulation of cytokine production", "regulation of macrophage activation", "inflammatory response",
                    "positive regulation of NF-kappaB transcription factor activity",
                   "neutrophil extravasation", "regulation of toll-like receptor signaling pathway", "regulation of lipid storage", 
                    "cell surface receptor signaling pathway", "B cell receptor signaling pathway", #same terms
           
          "positive regulation of nitric-oxide synthase biosynthetic process", "negative regulation of cytokine production",
          "interleukin-18-mediated signaling pathway", "production of molecular mediator involved in inflammatory response",
          "tumor necrosis factor production", "cytokine production involved in inflammatory response", "lipid homeostasis", 
          "regulation of long-chain fatty acid import across plasma membrane",
          
          "B cell proliferation", "B cell differentiation", "cell adhesion", "purine deoxyribonucleotide biosynthetic process")

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


selected_GO$Genes <- sub('.', '', selected_GO$Genes)

selected_GO %>%
  dplyr::select(GO.ID, Term, Annotated, Significant, geneRatio, Fisher.elim, Genes, direction) %>% 
  kbl() %>%
  kable_styling()


selected_GO <- selected_GO %>% #highlight terms that were the same as in the DEGs lfc 0.5+topvar GOA
  mutate(same_terms = ifelse(selected_GO$Term %in% same_terms, "same", "new"))# %>% 
 # mutate(color_term = ifelse(same_terms == "same", "grey60", "black"))


pdf("Output_files/Methylation/pred/GO/correlated_topvar/lfc0.5/GO_pred_lfc0.5_correlated_topVar_same_newTerms.pdf", width = 9, height = 10)
dot_plot <- ggplot(selected_GO, aes(x = direction, y= Term, size = geneRatio, color = Fisher.elim)) +
  geom_point() +
  scale_colour_gradient(low = "#990000", high = "#FF9999") +
  #scale_size_continuous(range=c(3,10)) +
  xlab('') + ylab('') +
  labs(
    title = "GO BP Pred x No Syst",
    subtitle = "Transcriptome-methylome integration lfc 0.5 + topVar",
    caption = paste0("Produced on ", Sys.time()),
    color = "p-value",
    size = "Gene ratio") +
  theme_bw() + #base_size = 20
  theme(
    axis.text.y = element_text(hjust = 1, size=12), # color = selected_GO$color_term
    axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1),
    axis.title = element_text(size = 20))+ 
  scale_y_discrete(limits = c("B cell receptor signaling pathway", #same neg
                              "B cell proliferation", "B cell differentiation", "cell adhesion", "purine deoxyribonucleotide biosynthetic process", #new neg

                              "cell surface receptor signaling pathway", #=

                              "positive regulation of cytokine production", "regulation of macrophage activation", "inflammatory response",
                               "positive regulation of NF-kappaB transcription factor activity",
                               "neutrophil extravasation", "regulation of toll-like receptor signaling pathway", #same pos
                              "positive regulation of nitric-oxide synthase biosynthetic process", "negative regulation of cytokine production",
                              "interleukin-18-mediated signaling pathway", "production of molecular mediator involved in inflammatory response",
                              "tumor necrosis factor production", "cytokine production involved in inflammatory response", "lipid homeostasis",
                              "regulation of long-chain fatty acid import across plasma membrane")) #new pos
                              


  dot_plot
dev.off()



