# EZE cohort: Therapy Signatures
  # Number of patients with inactive disease by comparisons and diagnosis
  # Supplementary Figure 1 
  # Author: Graziella Credidio

rm(list = ls())

folder <- "Output_files/Patients_distribution"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading packages ----
library(RColorBrewer)
library(tidyverse)

# Loading files ----
coldata <- read.table("Cleaned_tables/EZECohort_coldata_clean_ord.txt")

# inactive disease patients
aza_noSyst_coldata <- read.table("Cleaned_tables/models/aza/EZECohort_coldata_aza.txt", sep = "\t")
pred_noSyst_coldata <- read.table("Cleaned_tables/models/pred/EZECohort_coldata_pred.txt", sep = "\t")
antiTNF_noBio_coldata <- read.table("Cleaned_tables/models/antiTNF/EZECohort_coldata_antiTNF.txt", sep = "\t")

# all patients
all_aza_noSyst_coldata <- read.table("Cleaned_tables/models/aza/EZECohort_coldata_aza_allPatients.txt", sep = "\t")
all_pred_noSyst_coldata <- read.table("Cleaned_tables/models/pred/EZECohort_coldata_pred_allPatients.txt", sep = "\t")
all_antiTNF_noBio_coldata <- read.table("Cleaned_tables/models/antiTNF/EZECohort_coldata_antiTNF_allPatients.txt", sep = "\t")

# preprocessing existing coldatas ----
aza_noSyst_preprocessing <- function(coldata){
  aza_coldata <- coldata %>% 
    mutate(Comparison = "aza x no syst",
           Treatment = case_when(
             aza_vs_noSyst == 'aza' ~ "Aza",
             aza_vs_noSyst == "no_syst" ~ "No syst (aza)")) %>% 
    dplyr::select(study_id, diagnosis_class, Treatment, Comparison)
  return(aza_coldata)
}
aza_noSyst_coldata <- aza_noSyst_preprocessing(aza_noSyst_coldata)
all_aza_noSyst_coldata <- aza_noSyst_preprocessing(all_aza_noSyst_coldata)

pred_noSyst_preprocessing <- function(coldata){
  pred_coldata <- coldata %>% 
    mutate(Comparison = "pred x no syst",
           Treatment = case_when(
             pred_vs_noSyst == "pred" ~ "Pred",
             pred_vs_noSyst == "no_syst" ~ "No syst (pred)")) %>% 
    dplyr::select(study_id, diagnosis_class, Treatment, Comparison)
  return(pred_coldata)
}
pred_noSyst_coldata <- pred_noSyst_preprocessing(pred_noSyst_coldata)
all_pred_noSyst_coldata <- pred_noSyst_preprocessing(all_pred_noSyst_coldata)

antiTNF_preprocessing <- function(coldata){
  antiTNF_coldata <- coldata %>% 
    mutate(Comparison = "anti-TNF x no biologics",
           Treatment = case_when(
             antiTNF_vs_noBiologics == "anti_tnf" ~ "anti-TNF",
             antiTNF_vs_noBiologics == "no_biologics" ~ 'No biologics')) %>% 
    dplyr::select(study_id, diagnosis_class, Treatment, Comparison)
  return(antiTNF_coldata)
}
antiTNF_noBio_coldata <- antiTNF_preprocessing(antiTNF_noBio_coldata)
all_antiTNF_noBio_coldata <- antiTNF_preprocessing(all_antiTNF_noBio_coldata)

# Creating coldata for other models ----
# inactive disease coldata
coldata_inactive <- coldata %>% 
  filter(remission == "R") %>% 
  filter(crp < 5)

# Aza x no aza
aza_noAza_preprocessing <- function(df){
  aza_noAza_coldata <- df %>%
    filter(diagnosis_class == "CD" | diagnosis_class == "UC") %>% # only IBD patients
    mutate(aza_vs_noAza = case_when(
      (Azathioprin == "1") ~ "Aza",
      (Prednisolon == "1" |
         MTX == "1" |
         X6.MP == "1" |
         Sulfasalazin == "1" |
         Chloroquin == "1" |
         Hydroxychloroquin == "1" |
         Leflunomid == "1" |
         Ciclosporin == "1" |
         Mycophenolat.Mofetil == "1" |
         Cyclophosphamid == "1" |
         Tacrolimus == "1" |
         Apremilast == "1" |
         Mepacrine == "1" |
         No_syst == "1") ~ "No aza")) %>% 
    relocate(aza_vs_noAza, .after = No_syst) %>% 
    mutate(Comparison = "aza x no aza",
           Treatment = aza_vs_noAza) %>% 
    dplyr::select(study_id, diagnosis_class, Treatment, Comparison)
  return(aza_noAza_coldata)
}
aza_noAza_coldata <- aza_noAza_preprocessing(coldata_inactive)
all_aza_noAza_coldata <- aza_noAza_preprocessing(coldata)


# Pred x no pred
pred_noPred_preprocessing <- function(df){
  pred_noPred_coldata <- df %>% 
    mutate(pred_vs_noPred = case_when(
      (Prednisolon == "1") ~ "Pred",
      (Azathioprin == "1" |
         MTX == "1" |
         X6.MP == "1" |
         Sulfasalazin == "1" |
         Chloroquin == "1" |
         Hydroxychloroquin == "1" |
         Leflunomid == "1" |
         Ciclosporin == "1" |
         Mycophenolat.Mofetil == "1" |
         Cyclophosphamid == "1" |
         Tacrolimus == "1" |
         Apremilast == "1" |
         Mepacrine == "1" |
         No_syst == "1") ~ "No pred")) %>% 
    relocate(pred_vs_noPred, .after = No_syst) %>% 
    filter(!(diagnosis_class == "Pso")) %>%  # Exclude Pso patients
    mutate(Comparison = "pred x no pred",
           Treatment = pred_vs_noPred) %>% 
    dplyr::select(study_id, diagnosis_class, Treatment, Comparison)
  return(pred_noPred_coldata)
}
pred_noPred_coldata <- pred_noPred_preprocessing(coldata_inactive)
all_pred_noPred_coldata <- pred_noPred_preprocessing(coldata)


# Anti-TNF x non-anti-TNF
antiTNF_nonAntiTNF_preprocessing <- function(df){
  antiTNF_nonAntiTF_coldata <- df %>% 
    filter(biologics_TNF == "anti_tnf" | biologics_TNF == "non_anti_tnf") %>%
    dplyr::rename(antiTNF_vs_nonAntiTNF = biologics_TNF) %>% 
    mutate(Comparison = "anti-TNF x non-anti-TNF",
           Treatment = case_when(
             antiTNF_vs_nonAntiTNF == "anti_tnf" ~ "anti-TNF",
             antiTNF_vs_nonAntiTNF == "non_anti_tnf" ~ "Non-anti-TNF")) %>% 
    dplyr::select(study_id, diagnosis_class, Treatment, Comparison)
  return(antiTNF_nonAntiTF_coldata)
}
antiTNF_nonAntiTNF_coldata <- antiTNF_nonAntiTNF_preprocessing(coldata_inactive)
all_antiTNF_nonAntiTNF_coldata <- antiTNF_nonAntiTNF_preprocessing(coldata)

# Mtx x no Syst
mtx_noSyst_preprocessing <- function(df){
  mtx_noSyst_coldata <- df %>% 
    mutate(mtx_vs_noSyst = case_when(
      (MTX == "1")~ "mtx",
      (No_syst == "1") ~ "no_syst")) %>% 
    relocate(mtx_vs_noSyst, .after = No_syst) %>% 
    mutate(Comparison = "MTX x no syst",
           Treatment = case_when(
             mtx_vs_noSyst == 'mtx' ~ "MTX",
             mtx_vs_noSyst == 'no_syst' ~ "No syst (MTX)")) %>% 
    filter(diagnosis_class == "RA" | diagnosis_class == "PsA") %>% # only RA and PsA patients
    filter(!(is.na(mtx_vs_noSyst))) %>% 
    dplyr::select(study_id, diagnosis_class, Treatment, Comparison)
  return(mtx_noSyst_coldata)
}
mtx_noSyst_coldata <- mtx_noSyst_preprocessing(coldata_inactive)
all_mtx_noSyst_coldata <- mtx_noSyst_preprocessing(coldata)

# Mtx x no mtx
mtx_noMtx_preprocessing <- function(df){
  mtx_noMtx_coldata <- df %>% 
    mutate(mtx_vs_noMtx = case_when(
      (MTX == "1") ~ "mtx",
      (Prednisolon == "1" |
         Azathioprin == "1" |
         X6.MP == "1" |
         Sulfasalazin == "1" |
         Chloroquin == "1" |
         Hydroxychloroquin == "1" |
         Leflunomid == "1" |
         Ciclosporin == "1" |
         Mycophenolat.Mofetil == "1" |
         Cyclophosphamid == "1" |
         Tacrolimus == "1" |
         Apremilast == "1" |
         Mepacrine == "1" |
         No_syst == "1") ~ "no_mtx")) %>% 
    relocate(mtx_vs_noMtx, .after = No_syst) %>% 
    mutate(Comparison = "MTX x no MTX",
           Treatment = case_when(
             mtx_vs_noMtx == "mtx" ~ "MTX",
             mtx_vs_noMtx == "no_mtx" ~ "No MTX")) %>% 
    filter(diagnosis_class == "RA" | diagnosis_class == "PsA") %>% # only RA and PsA patients
    filter(!(is.na(mtx_vs_noMtx))) %>% 
    dplyr::select(study_id, diagnosis_class, Treatment, Comparison)
  return(mtx_noMtx_coldata)
}
mtx_noMtx_coldata <- mtx_noMtx_preprocessing(coldata_inactive)
all_mtx_noMtx_coldata <- mtx_noMtx_preprocessing(coldata)

# Pre processing ----
inactive_coldata <- rbind(antiTNF_noBio_coldata, antiTNF_nonAntiTNF_coldata, 
                     pred_noSyst_coldata, pred_noPred_coldata,
                     aza_noSyst_coldata, aza_noAza_coldata,
                     mtx_noSyst_coldata, mtx_noMtx_coldata)

all_coldata <- rbind(all_antiTNF_noBio_coldata, all_antiTNF_nonAntiTNF_coldata, 
                     all_pred_noSyst_coldata, all_pred_noPred_coldata,
                     all_aza_noSyst_coldata, all_aza_noAza_coldata,
                     all_mtx_noSyst_coldata, all_mtx_noMtx_coldata)

# Plots ----
# Supplementary Figure 1
plot_comparisons_diagnosis <- function(df, title, filePath ){
  breaks <- seq(from = 0, to = 200, by = 25)
  plot <- ggplot(df, aes(Comparison, fill = diagnosis_class)) +
    geom_bar() + theme_bw() +
    geom_text(aes(label = after_stat(count)), stat = "count", hjust = 0.5, position = position_stack(vjust = 0.5), colour = "#808080", size = 5) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.title.x=element_blank(),
      axis.title.y =element_blank(),
      title = element_text(size = 24),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 14)
      ) +
    scale_y_continuous(breaks = breaks) +
    labs(fill= 'Diagnosis',
         title = title
         ) +
    scale_fill_brewer(type = "qual", palette = "Pastel2") +
    scale_x_discrete(limits= c("anti-TNF x no biologics", "anti-TNF x non-anti-TNF", 
                               "aza x no aza", "aza x no syst", 
                               "MTX x no MTX", "MTX x no syst", 
                               "pred x no pred", "pred x no syst"))
  ggsave(plot, file = filePath, height = 8, width = 15, units = "in", dpi = 300)
  plot
}

plot_comparisons_diagnosis(all_coldata, "All patients", "Output_files/Patients_distribution/allPatients_byComparisons_byDignosis.png")
plot_comparisons_diagnosis(inactive_coldata, "Inactive patients", "Output_files/Patients_distribution/inactive_byComparisons_byDignosis.png")


