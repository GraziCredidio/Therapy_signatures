# EZE cohort: Therapy Signatures
  # Principal Component Analysis
  # Author: Graziella Credidio

rm(list = ls())

folder <- "Output_files/PCA"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading packages ----
library(data.table)
library(tidyverse)
library(ggrepel)
library(DESeq2)
library(RColorBrewer)

# Loading files ----
counts <- read.table("Cleaned_tables/EZECohort_counts_ord.txt", sep = "\t")
coldata <- read.table("Cleaned_tables/EZECohort_coldata_clean_ord.txt", sep = "\t")

# Data preprocessing to be plotted ----
coldata <- coldata %>% 
  mutate(Syst_therapy = recode(No_syst,
                               "0" = "Yes",
                               "1" = "No"))

coldata <- coldata %>% 
  mutate(Prednisolon = recode(Prednisolon,
                               "0" = "Yes",
                               "1" = "No"))

# Transforming categorical columns as factor
coldata <- coldata%>%
  mutate(across(c(diagnosis_class, age_group, sex, remission, bmi_class, Syst_therapy,
                  smoking, biologics, biologics_TNF), as.factor))

# Counts normalization: DESeq2 ----
dds_counts <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = coldata,
                                     design = ~ diagnosis_class + sex + age_group + bmi_class + biologics)
dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(counts), ]
dds_counts <- estimateSizeFactors(dds_counts)

vst <- vst(dds_counts, blind=TRUE) 
vst_counts_norm <- assay(vst)

# PCA ----
pc <- prcomp(t(vst_counts_norm))
pc_df <- as.data.frame(pc$x)

# PCs dataframe
pc_df <- cbind(pc_df, coldata)

# Sequencing depth and PC percentage
pc_df$depth <- colSums(counts)/1e6 

pc.var <- pc$sdev^2 
pc.per <- round(pc.var/sum(pc.var)*100, 1)

pc_df <- pc_df%>%
  mutate(across(c(diagnosis_class, age_group, sex, remission, bmi_class, Syst_therapy, No_syst,
                  smoking, biologics, biologics_TNF, Azathioprin, Prednisolon, MTX), as.factor))

# PCA plots ----
PCA_plot <- function(data, x, y, plotColor, plotLabel= NULL, legendTitleColor, legendTitleShape, filePath, plotShape) {
  if (is.null(plotLabel)) {
    aes <- aes(x, y, color = plotColor, shape = plotShape)
    labs(color=legendTitleColor)
  } else {
    aes <- aes(x, y, color = plotColor, label = plotLabel, shape = plotShape)
  }
  
  plot <- ggplot(data, aes) +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() + theme(axis.text = element_text(size = 15),
                       axis.title = element_text(size = 15),
                       legend.title = element_text(size = 10),
                       legend.text = element_text(size = 10))+
    xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
    
    # outliers boundaries:
    #geom_vline(xintercept = (mean(pc_df$PC1) - 3*sd(pc_df$PC1))) +
    #geom_vline(xintercept = (mean(pc_df$PC1) + 3*sd(pc_df$PC1))) +
    #geom_hline(yintercept = (mean(pc_df$PC2) - 3*sd(pc_df$PC2))) +
    #geom_hline(yintercept = (mean(pc_df$PC2) + 3*sd(pc_df$PC2))) +
    
    ylab(paste0("PC2 (",pc.per[2],"%",")")) +
    coord_fixed() +
    labs(color= legendTitleColor, shape = legendTitleShape)  
  ggsave(plot, file = filePath, height = 8, width = 10, units = "in", dpi = 300)
  plot
}

# Diagnosis and sex
PCA_plot(data = pc_df, 
         x = pc_df$PC1, 
         y = pc_df$PC2, 
         plotColor = pc_df$diagnosis_class, 
         legendTitleColor = "Diagnosis", 
         legendTitleShape = "Sex",
         filePath = 'Output_files/PCA/PCA_diag_sex.png',
         plotShape = pc_df$sex)

# BMI and sex
PCA_plot(data = pc_df, 
         x = pc_df$PC1, 
         y = pc_df$PC2, 
         plotColor = pc_df$bmi, 
         legendTitleColor = "BMI", 
         legendTitleShape = "Sex",
         filePath = 'Output_files/PCA/PCA_bmi_sex.png',
         plotShape = pc_df$sex)

# Age and sex
PCA_plot(data = pc_df, 
         x = pc_df$PC1, 
         y = pc_df$PC2, 
         plotColor = pc_df$age, 
         legendTitleColor = "Age (years)", 
         legendTitleShape = "Sex",
         filePath = 'Output_files/PCA/PCA_age_sex.png',
         plotShape = pc_df$sex)

# log(crp) and sex
PCA_plot(data = pc_df, 
         x = pc_df$PC1, 
         y = pc_df$PC2, 
         plotColor = pc_df$crp_log, 
         legendTitleColor = "log2(CRP)", 
         legendTitleShape = "Sex",
         filePath = 'Output_files/PCA/PCA_crp_sex.png',
         plotShape = pc_df$sex)

# Prednisolone and sex
PCA_plot(data = pc_df, 
         x = pc_df$PC1, 
         y = pc_df$PC2, 
         plotColor = pc_df$Prednisolon, 
         legendTitleColor = "Prednisolone", 
         legendTitleShape = "Sex",
         filePath = 'Output_files/PCA/PCA_pred_sex.png',
         plotShape = pc_df$sex)

# Systemic therapies and sex
PCA_plot(data = pc_df, 
         x = pc_df$PC1, 
         y = pc_df$PC2, 
         plotColor = pc_df$Syst_therapy, 
         legendTitleColor = "Systemic therapies", 
         legendTitleShape = "Sex",
         filePath = 'Output_files/PCA/PCA_systTherapies_sex.png',
         plotShape = pc_df$sex)

# PC1 x PC3 by sex
PC1_PC3_plot <- 
  ggplot(pc_df, aes(PC2, PC3, color = sex)) + 
  geom_point(size = 5, alpha = 0.8) +
  theme_bw() + theme(axis.text = element_text(size = 15),
                     axis.title = element_text(size = 20),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 15))+
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[3],"%",")")) +
  labs(color= "Sex") +
  coord_fixed() 
PC1_PC3_plot
ggsave(PC1_PC3_plot, file = 'Output_files/PCA/PCA_PC1_PC3_sex.png', height = 8, width = 10, units = "in", dpi = 300)


