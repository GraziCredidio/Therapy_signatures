# EZE cohort: Therapy Signatures
  # Number of patients with inactive disease by treatment and diagnosis
  # Figure 1 
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(RColorBrewer)
library(tidyverse)

# Loading files ----
aza_noSyst_coldata <- read.table("Cleaned_tables/models/aza/EZECohort_coldata_aza.txt", sep = "\t")
pred_noSyst_coldata <- read.table("Cleaned_tables/models/pred/EZECohort_coldata_pred.txt", sep = "\t")
antiTNF_noBio_coldata <- read.table("Cleaned_tables/models/antiTNF/EZECohort_coldata_antiTNF.txt", sep = "\t")

coldata <- read.table("Cleaned_tables/EZECohort_coldata_clean_ord.txt")

# coldata with only inactive disease patients----
coldata <- coldata %>% 
  filter(remission == "R") %>% 
  filter(crp < 5)

# preprocessing existing coldatas ----
aza_noSyst_coldata <- aza_noSyst_coldata %>% 
  mutate(Comparison = "aza x no syst",
         Treatment = case_when(
           aza_vs_noSyst == 'aza' ~ "Aza",
           aza_vs_noSyst == "no_syst" ~ "No syst (aza)")) %>% 
  dplyr::select(study_id, diagnosis_class, Treatment, Comparison)

pred_noSyst_coldata <- pred_noSyst_coldata %>% 
  mutate(Comparison = "pred x no syst",
         Treatment = case_when(
           pred_vs_noSyst == "pred" ~ "Pred",
           pred_vs_noSyst == "no_syst" ~ "No syst (pred)")) %>% 
  dplyr::select(study_id, diagnosis_class, Treatment, Comparison)

antiTNF_noBio_coldata <- antiTNF_noBio_coldata %>% 
  mutate(Comparison = "anti-TNF",
         Treatment = case_when(
           antiTNF_vs_noBiologics == "anti_tnf" ~ "anti-TNF",
           antiTNF_vs_noBiologics == "no_biologics" ~ 'No biologics')) %>% 
  dplyr::select(study_id, diagnosis_class, Treatment, Comparison)
  
# Creating coldata for other models ----
# Aza x no aza
aza_noAza_coldata <- coldata %>%
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
  filter(aza_vs_noAza == "No aza") %>% 
  dplyr::select(study_id, diagnosis_class, Treatment, Comparison)

# Pred x no pred
pred_noPred_coldata <- coldata %>% 
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
  filter(pred_vs_noPred == "No pred") %>% 
  dplyr::select(study_id, diagnosis_class, Treatment, Comparison)

# Anti-TNF x non-anti-TNF
antiTNF_nonAntiTF_coldata <- coldata %>% 
  filter(biologics_TNF == "anti_tnf" | biologics_TNF == "non_anti_tnf") %>%
  dplyr::rename(antiTNF_vs_nonAntiTNF = biologics_TNF) %>% 
  mutate(Comparison = "anti-TNF x non-anti-TNF",
         Treatment = case_when(
           antiTNF_vs_nonAntiTNF == "anti_tnf" ~ "anti-TNF",
           antiTNF_vs_nonAntiTNF == "non_anti_tnf" ~ "Non-anti-TNF")) %>% 
  filter(antiTNF_vs_nonAntiTNF == "non_anti_tnf") %>% 
  dplyr::select(study_id, diagnosis_class, Treatment, Comparison)

# Mtx x no Syst
mtx_noSyst_coldata <- coldata %>% 
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

# Mtx x no mtx
mtx_noMtx_coldata <- coldata %>% 
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
  filter(mtx_vs_noMtx == "no_mtx") %>% 
  dplyr::select(study_id, diagnosis_class, Treatment, Comparison)

# Pre processing ----
all_coldata <- rbind(antiTNF_noBio_coldata, antiTNF_nonAntiTF_coldata, 
                     pred_noSyst_coldata, pred_noPred_coldata,
                     aza_noSyst_coldata, aza_noAza_coldata,
                     mtx_noSyst_coldata, mtx_noMtx_coldata)


# Plots ----
# Figure 1
breaks <- seq(from = 0, to = 200, by = 25)
plot_remitters_patients_treatments <- ggplot(all_coldata, aes(Treatment, fill = diagnosis_class)) +
  geom_bar() + theme_bw() +
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = 0.5, position = position_stack(vjust = 0.5), colour = "#808080", size = 5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x=element_blank(),
    axis.title.y =element_blank(),
    title = element_text(size = 24),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  ) +
  scale_y_continuous(breaks = breaks) +
  labs(fill= 'Diagnosis',
       title = "Inactive disease patients") +
  scale_fill_brewer(type = "qual", palette = "Pastel2") +
  scale_x_discrete(limits= c("anti-TNF", "No biologics", "Non-anti-TNF", "Aza",
                            "No aza", "No syst (aza)", "MTX", "No MTX", "No syst (MTX)", 
                            "Pred", "No pred", "No syst (pred)"))
plot_remitters_patients_treatments 
ggsave(plot_remitters_patients_treatments,
       file = "Output_files/Patients_distribution/inactive_byTreatment_byDiagnosis.png",
       height = 8, width = 15, units = "in", dpi = 300)
