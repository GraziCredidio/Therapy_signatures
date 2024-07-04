# EZE Cohort - therapy signatures
# Anti TNF patients: time of treatment

coldata_R_noBio <- read.table("Output_files/Maaslin2/Tables/Remission_crp/coldata_maaslin2_antiTNF_R_crp_04.05.txt", sep = "\t")

coldata_antiTNF <- coldata_R_noBio %>% 
  filter(antiTNF_vs_noBiologics == "anti_tnf") 

anti_tnf_syst <- coldata_antiTNF %>% 
  filter(No_syst == 0)


EZE <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")
EZE_antiTNF <- EZE[EZE$study_id %in% coldata_antiTNF$study_id, ]


date_EZE_antiTNF <- EZE_antiTNF[,c(1, 26, 275, 276)]
visit <- as.Date(date_EZE_antiTNF$date_visit)
first_therapy <- as.Date(date_EZE_antiTNF$therapy_date_first)
date_EZE_antiTNF$biologics_duration <- difftime(visit, first_therapy)

date_EZE_antiTNF <- EZE_antiTNF[,c(1, 26, 275, 276)]
date_EZE_antiTNF[10,3] <- "2009-11-01"

date_EZE_antiTNF$biologics_duration_years <- date_EZE_antiTNF$biologics_duration / 365
date_EZE_antiTNF$biologics_duration_weeks <- round(date_EZE_antiTNF$biologics_duration / 7, 2)
date_EZE_antiTNF$biologics_duration_months <- round(date_EZE_antiTNF$biologics_duration / 30, 2)

date_EZE_antiTNF <- date_EZE_antiTNF %>% 
  filter(!(biologics_duration == "NA"))
min(date_EZE_antiTNF$biologics_duration)
median(date_EZE_antiTNF$biologics_duration)

ggplot(date_EZE_antiTNF, aes(x = biologics_duration_weeks)) +
  geom_histogram (
    breaks = seq(0, 400, by = 14), 
    binwidth=1, center= 0,
    aes(fill = ..count..)) + 
  scale_x_continuous(breaks = seq(0, 400, by=14)) +
  scale_y_continuous(breaks = seq(0,10, by =1)) +
  theme_bw() +
  xlab("Weeks since first biologics use") +
  ylab("# of patients") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=45, hjust=1)) 


ggplot(date_EZE_antiTNF, aes(x = biologics_duration_months)) +
  geom_histogram (
    breaks = seq(0, 95, by = 1), 
    binwidth=1, center= 0,
    aes(fill = ..count..)) + 
  scale_x_continuous(breaks = seq(0, 95, by=1)) +
  theme_bw() +
  xlab("Months since first biologics use") +
  ylab("# of patients") +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=90, hjust=1)) 



## IG genes exploration
sig_genes_antiTNF <- read.table("Output_files/Maaslin2/Results/Remission_crp/significant_results/maaslin2_significant_results_antiTNF_vs_noBio_R_crp.txt", sep = "\t")

IG_genes <- sig_genes_antiTNF %>% 
  filter(grepl('IG', genes)) 

IGK <- IG_genes %>% 
  filter(grepl('^IGK', genes)) 

IGL <- IG_genes %>% 
  filter(grepl('^IGL', genes)) 

IGH <- IG_genes %>% 
  filter(grepl('^IGH'), genes) 


patterns <- c("^IGK", "^IGL", "^IGH")
IG_V <- IG_genes %>% 
  filter(grepl(paste0(patterns, 'V', collapse= "|"), genes))


antiTNF_up_sig_gene_correlation <- read.table("Output_files/Methylation/antiTNF/antiTNF_up_significant_DEG_meth_site_correlation_5000bp_1000rep.txt",
                                              sep = '\t')
antiTNF_down_sig_gene_correlation <- read.table("Output_files/Methylation/antiTNF/antiTNF_down_significant_DEG_meth_site_correlation_5000bp_1000rep.txt", 
                                                sep = '\t')

IG_methyalted <- antiTNF_up_sig_gene_correlation %>% 
  filter(grepl('IG', Gene_name))
unique(IG_methyalted$Gene_name)

TVGs <- read.table("Cleaned_tables/topVar_2000_genes.txt", sep = "\t")

up_corr_tvg <- intersect(antiTNF_up_sig_gene_correlation$Gene, rownames(TVGs))
up_corr_tvg <- antiTNF_up_sig_gene_correlation[antiTNF_up_sig_gene_correlation$Gene %in% up_corr_tvg,]
unique(up_corr_tvg$Gene_name)
