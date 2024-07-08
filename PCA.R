# PCA - EZE cohort

graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("D:\\Documentos/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

# Loading packages ----
library(data.table)
library(tidyverse)
library(ggrepel)
library(DESeq2)
library(RColorBrewer)
library(Rgraphviz)
library(plotly)

folder <- "Output_files/PCA"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading files ----
coldata <- read.csv("Cleaned_tables/EZECohort_ord.coldata_maaslin_16.03.txt", sep = "\t")
counts <- read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_16.03.txt", sep = "\t")

# Data preprocessing ----
# NAs
coldata$smoking <- str_replace_na(coldata$smoking, "NA")
coldata$remission <- str_replace_na(coldata$remission, "NA")
coldata$bmi_class <- str_replace_na(coldata$bmi_class, "NA")

coldata <- coldata %>% 
  mutate(Syst_therapy = recode(No_syst,
                               "0" = "Yes",
                               "1" = "No"))

# As factor
coldata <- coldata%>%
  mutate(across(c(diagnosis_class, age_group, age_group2, sex, remission, bmi_class, Syst_therapy,
                  smoking, biologics, biologics_TNF), as.factor))

# DESeq2 object----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = coldata,
                                     design = ~ diagnosis_class + sex + age_group + 
                                       bmi_class + biologics)
dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

# PCA ----
pc <- prcomp(t(vst_counts_norm))
pc_df <- as.data.frame(pc$x)
summary(pc)

pc_df <- cbind(pc_df, coldata)
pc_df$depth <- colSums(counts)/1e6

pc.var <- pc$sdev^2 
pc.per <- round(pc.var/sum(pc.var)*100, 1)

pc_df <- pc_df%>%
  mutate(across(c(diagnosis_class, age_group, age_group2, sex, remission, bmi_class, Syst_therapy, No_syst,
                  smoking, biologics, biologics_TNF, Azathioprin, Prednisolon, MTX), as.factor))
# PCA plots ----
# Simple PCA plot
simple_plot <- 
  ggplot(pc_df, aes(pc_df$PC1, pc_df$PC2, color = age_group2))+ #, label = study_id)) +
 # geom_text(hjust = 1, nudge_x = 2, nudge_y = 3, cex = 3) +
  geom_vline(xintercept = (mean(pc_df$PC1) - 3*sd(pc_df$PC1))) +
  geom_vline(xintercept = (mean(pc_df$PC1) + 3*sd(pc_df$PC1))) +
  geom_hline(yintercept = (mean(pc_df$PC2) - 3*sd(pc_df$PC2))) +
  geom_hline(yintercept = (mean(pc_df$PC2) + 3*sd(pc_df$PC2))) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() + theme(axis.text = element_text(size = 15),
                      axis.title = element_text(size = 15),
                     legend.title = element_text(size = 10),
                     legend.text = element_text(size = 10))+
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  labs(color= "systemic therapy")
simple_plot

simple_plot <- simple_plot +  stat_ellipse() 
simple_plot <- simple_plot + scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(7)[1:9])
simple_plot


ggsave(simple_plot, file = 'Output_files/PCA/No_outliers/PCA_syst_therapy_ellipse.png', height = 8, width = 10, units = "in", dpi = 300)


# More complete PCA plot
PCA_plot <- function(data, x, y, plotColor, plotLabel= NULL, legendTitleColor, legendTitleShape, filePath, plotShape) {
  if (is.null(plotLabel)) {
    aes <- aes(x, y, color = plotColor, shape = plotShape)
    labs(color=legendTitleColor)
  } else {
    aes <- aes(x, y, color = plotColor, label = plotLabel, shape = plotShape)
  }
  
  plot <- ggplot(data, aes) +
    #geom_text(hjust = 0, nudge_x = 2, nudge_y = 1, cex = 3) +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() + theme(axis.text = element_text(size = 15),
                       axis.title = element_text(size = 15),
                       legend.title = element_text(size = 10),
                       legend.text = element_text(size = 10))+
    xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
    #geom_vline(xintercept = (mean(pc_df$PC1) - 3*sd(pc_df$PC1))) +
    #geom_vline(xintercept = (mean(pc_df$PC1) + 3*sd(pc_df$PC1))) +
    #geom_hline(yintercept = (mean(pc_df$PC2) - 3*sd(pc_df$PC2))) +
    #geom_hline(yintercept = (mean(pc_df$PC2) + 3*sd(pc_df$PC2))) +
    ylab(paste0("PC2 (",pc.per[2],"%",")")) +
    coord_fixed() +
    labs(color= legendTitleColor, shape = legendTitleShape) #+ 
  # scale_color_manual(values = colors) #+
  #stat_ellipse() 
  ggsave(plot, file = filePath, height = 8, width = 10, units = "in", dpi = 300)
  plot
}

PCA_plot(data = pc_df, 
         x = pc_df$PC1, 
         y = pc_df$PC2, 
         plotColor = pc_df$Syst_therapy, 
         #plotLabel = pc_df$biologics,
         legendTitleColor = "Use of Systemic Therapies", 
         legendTitleShape = "Sex",
         filePath = 'Output_files/PCA/01.06/PCA_noSyst_sex.png',
         plotShape = pc_df$sex)



range(pc_df$crp)
# plotPCA() function (top 500 genes) ----
DESeq2::plotPCA(vst, intgroup="age_group") #500 top genes, selected by highest row variance
# # calculate the variance for each gene
# rv <- rowVars(assay(object))
# 
# # select the ntop genes by variance
# select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
# 
# # perform a PCA on the data in assay(x) for the selected genes
# pca <- prcomp(t(assay(object)[select,]))
# 
# # the contribution to the total variance for each component
# percentVar <- pca$sdev^2 / sum( pca$sdev^2 )


# Correlation ----
cor1 <- cor.test(pc_df$PC1, pc_df$crp,  method = "spearman", exact = FALSE)
cor1
round(cor1$estimate,1)


cor1 <- cor.test(pc_df$PC3, pc_df$crp,  method = "spearman", exact = FALSE)
cor1

cor2 <- cor.test(pc_df$PC1, pc_df$depth, method = "spearman", exact = FALSE)
cor2

cor3 <- cor.test(pc_df$PC1, pc_df$age, method = "spearman", exact = FALSE)
round(cor3$p.value,1)

cor4 <- cor.test(pc_df$PC1, pc_df$bmi, method = "spearman", exact = FALSE)
cor4
# Genes contributions to PC3 ----
eigenvalues3  <- data.frame(sort(abs(pc$rotation[,"PC3"]), decreasing=TRUE)[1:500]) 


####################################### PCA for each comparison###########################################################
# Loading files ----
counts <-  read.csv("Cleaned_tables/EZECohort_ord.counts_maaslin_28.02.txt", sep="\t")

coldata_noBio <- read.table("Output_files/Maaslin2/Tables/coldata_maaslin2_antiTNF_vs_noBiologics_28.02.txt", sep = "\t")
coldata_nonTNF <- read.table("Output_files/Maaslin2/Tables/coldata_maaslin2_antiTNF_vs_nonAntiTNF_28.02.txt", sep = "\t")

coldata_noPred <- read.table("Output_files/Maaslin2/Tables/coldata_maaslin2_pred_vs_noPred_28.02.txt", sep = "\t")
coldata_pred_noSyst <- read.table("Output_files/Maaslin2/Tables/coldata_maaslin2_pred_vs_noSyst_28.02.txt", sep = "\t")

coldata_noAza <- read.table("Output_files/Maaslin2/Tables/coldata_maaslin2_aza_vs_noAza_IBD_28.02.txt", sep = "\t")
coldata_aza_noSyst <- read.table("Output_files/Maaslin2/Tables/coldata_maaslin2_aza_vs_noSyst_IBD_28.02.txt", sep = "\t")

coldata_noMtx <- read.table("Output_files/Maaslin2/Tables/coldata_maaslin2_mtx_vs_noMtx_RA-PsA_28.02.txt", sep = "\t")
coldata_mtx_noSyst <- read.table("Output_files/Maaslin2/Tables/coldata_maaslin2_mtx_vs_noSyst_RA-PsA_28.02.txt", sep = "\t")

#PCA ----
## PCA tnf x no Bio
counts_noBio <- counts[,colnames(counts) %in% coldata_noBio$sample_id]
all(rownames(coldata_noBio) == colnames(counts_noBio))

dds_counts <- DESeqDataSetFromMatrix(countData = counts_noBio,
                                     colData = coldata_noBio,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group2 +
                                       antiTNF_vs_noBiologics)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

pc <- prcomp(t(vst_counts_norm))
pc_df <- as.data.frame(pc$x)
pc_df$depth <- colSums(counts_noBio)/1e6
pc.var <- pc$sdev^2 
pc.per <- round(pc.var/sum(pc.var)*100, 1)
pc_df <- cbind(pc_df, coldata_noBio)

simple_plot <- 
  ggplot(pc_df, aes(PC1, PC2, color = bmi_class)) + #, label = study_id)) +
  # geom_text(hjust = 1, nudge_x = 2, nudge_y = 3, cex = 3) +
  geom_vline(xintercept = (mean(pc_df$PC1) - 3*sd(pc_df$PC1))) +
  geom_vline(xintercept = (mean(pc_df$PC1) + 3*sd(pc_df$PC1))) +
  geom_hline(yintercept = (mean(pc_df$PC2) - 3*sd(pc_df$PC2))) +
  geom_hline(yintercept = (mean(pc_df$PC2) + 3*sd(pc_df$PC2))) +
  geom_point(size = 5, alpha = 0.8) +
  theme_bw() + theme(axis.text = element_text(size = 15),
                     axis.title = element_text(size = 20),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 15))+
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title = "PCA anti TNF x no Biologics",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  labs(color= "BMI")
simple_plot
simple_plot <- simple_plot + scale_color_brewer(palette = "Oranges")
simple_plot
ggsave(simple_plot, file = 'Output_files/PCA/antiTNF/PCA_antiTNF_noBio_bmi.png', 
       height = 8, width = 10, units = "in", dpi = 300)

## PCA pred x no syst
counts_pred_noSyst <- counts[,colnames(counts) %in% coldata_pred_noSyst$sample_id]
all(rownames(coldata_pred_noSyst) == colnames(counts_pred_noSyst))

dds_counts <- DESeqDataSetFromMatrix(countData = counts_pred_noSyst,
                                     colData = coldata_pred_noSyst,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group2 +
                                       pred_vs_noSyst)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

pc <- prcomp(t(vst_counts_norm))
pc_df <- as.data.frame(pc$x)
pc_df$depth <- colSums(counts_pred_noSyst)/1e6
pc.var <- pc$sdev^2 
pc.per <- round(pc.var/sum(pc.var)*100, 1)
pc_df <- cbind(pc_df, coldata_pred_noSyst)

simple_plot <- 
  ggplot(pc_df, aes(PC1, PC2, color = bmi_class)) + #, label = study_id)) +
  # geom_text(hjust = 1, nudge_x = 2, nudge_y = 3, cex = 3) +
  geom_vline(xintercept = (mean(pc_df$PC1) - 3*sd(pc_df$PC1))) +
  geom_vline(xintercept = (mean(pc_df$PC1) + 3*sd(pc_df$PC1))) +
  geom_hline(yintercept = (mean(pc_df$PC2) - 3*sd(pc_df$PC2))) +
  geom_hline(yintercept = (mean(pc_df$PC2) + 3*sd(pc_df$PC2))) +
  geom_point(size = 5, alpha = 0.8) +
  theme_bw() + theme(axis.text = element_text(size = 15),
                     axis.title = element_text(size = 20),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 15))+
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title = "PCA pred x no systemic therapies",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  labs(color= "BMI")
simple_plot
simple_plot <- simple_plot + scale_color_brewer(palette = "Oranges")
simple_plot
ggsave(simple_plot, file = 'Output_files/PCA/Pred/PCA_pred_noSyst_age.png', 
       height = 8, width = 10, units = "in", dpi = 300)



## PCA mtx x no syst
counts_mtx_noSyst <- counts[,colnames(counts) %in% coldata_mtx_noSyst$sample_id]
all(rownames(coldata_mtx_noSyst) == colnames(counts_mtx_noSyst))

dds_counts <- DESeqDataSetFromMatrix(countData = counts_mtx_noSyst,
                                     colData = coldata_mtx_noSyst,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group2 +
                                       mtx_vs_noSyst)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

pc <- prcomp(t(vst_counts_norm))
pc_df <- as.data.frame(pc$x)
pc_df$depth <- colSums(counts_mtx_noSyst)/1e6
pc.var <- pc$sdev^2 
pc.per <- round(pc.var/sum(pc.var)*100, 1)
pc_df <- cbind(pc_df, coldata_mtx_noSyst)

simple_plot <- 
  ggplot(pc_df, aes(PC1, PC3, color = diagnosis_class, shape = mtx_vs_noSyst)) + #, label = study_id)) +
  # geom_text(hjust = 1, nudge_x = 2, nudge_y = 3, cex = 3) +
  geom_vline(xintercept = (mean(pc_df$PC1) - 3*sd(pc_df$PC1))) +
  geom_vline(xintercept = (mean(pc_df$PC1) + 3*sd(pc_df$PC1))) +
  geom_hline(yintercept = (mean(pc_df$PC2) - 3*sd(pc_df$PC2))) +
  geom_hline(yintercept = (mean(pc_df$PC2) + 3*sd(pc_df$PC2))) +
  geom_point(size = 5, alpha = 0.8) +
  theme_bw() + theme(axis.text = element_text(size = 15),
                     axis.title = element_text(size = 20),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 15))+
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[3],"%",")")) +
  labs(title = "PCA mtx x no systemic therapies",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  labs(color= "Diagnosis")
simple_plot

simple_plot <- simple_plot + scale_color_brewer(palette = "Oranges")
simple_plot

ggsave(simple_plot, file = 'Output_files/PCA/MTX/PCA3_mtx_noSyst_diagnosis.png', 
       height = 8, width = 10, units = "in", dpi = 300)


## Aza x no Syst
counts_aza_noSyst <- counts[,colnames(counts) %in% coldata_aza_noSyst$sample_id]
all(rownames(coldata_aza_noSyst) == colnames(counts_aza_noSyst))

dds_counts <- DESeqDataSetFromMatrix(countData = counts_aza_noSyst,
                                     colData = coldata_aza_noSyst,
                                     design = ~ diagnosis_class + sex + bmi_class + age_group2 +
                                       aza_vs_noSyst)

dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

pc <- prcomp(t(vst_counts_norm))
pc_df <- as.data.frame(pc$x)
pc_df$depth <- colSums(counts_aza_noSyst)/1e6
pc.var <- pc$sdev^2 
pc.per <- round(pc.var/sum(pc.var)*100, 1)
pc_df <- cbind(pc_df, coldata_aza_noSyst)

simple_plot <- 
  ggplot(pc_df, aes(PC1, PC2, color = diagnosis_class)) + #, label = study_id)) +
  # geom_text(hjust = 1, nudge_x = 2, nudge_y = 3, cex = 3) +
  geom_vline(xintercept = (mean(pc_df$PC1) - 3*sd(pc_df$PC1))) +
  geom_vline(xintercept = (mean(pc_df$PC1) + 3*sd(pc_df$PC1))) +
  geom_hline(yintercept = (mean(pc_df$PC2) - 3*sd(pc_df$PC2))) +
  geom_hline(yintercept = (mean(pc_df$PC2) + 3*sd(pc_df$PC2))) +
  geom_point(size = 5, alpha = 0.8) +
  theme_bw() + theme(axis.text = element_text(size = 15),
                     axis.title = element_text(size = 20),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 15))+
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title = "PCA aza x no systemic therapies",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  labs(color= "Diagnosis")
simple_plot

simple_plot <- simple_plot + scale_color_brewer(palette = "Oranges")
simple_plot

ggsave(simple_plot, file = 'Output_files/PCA/Aza/PCA_aza_noSyst_crp.png', 
       height = 8, width = 10, units = "in", dpi = 300)
