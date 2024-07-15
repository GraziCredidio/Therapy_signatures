# EZE cohort
# Bar plot with number of patients in each comparison by diagnosis (remitters + crp)

graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

library(RColorBrewer)
library(tidyverse)

# Loading files ----
#inactive disease patients
setwd("C:/Documents/Masters thesis/EZE_cohort/Output_files/Maaslin2/Tables/Remission_crp") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort/Output_files/Maaslin2/Tables/Remission_crp") #pc

files_R = list.files(pattern="*.txt")
coldata_R = lapply(files_R, read.delim)

coldata_noBio <- coldata_R[[4]]
coldata_aza_noSyst <- coldata_R[[5]]
coldata_pred_noSyst <- coldata_R[[6]]

# all patients
setwd("C:/Documents/Masters thesis/EZE_cohort/Output_files/Maaslin2/Tables") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort/Output_files/Maaslin2/Tables") #pc

files = list.files(pattern="*.txt")
coldata = lapply(files, read.delim)

coldata_nonTNF <- coldata[[3]]
coldata_noAza <- coldata[[4]]
coldata_noMtx <- coldata[[6]]
coldata_mtx_noSyst <- coldata[[7]]
coldata_noPred <- coldata[[8]]

# filtering to inactive + crp
setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

coldata_nonTNF <- coldata_nonTNF %>% 
  filter(remission == "R" & crp < 5)

coldata_noAza <- coldata_noAza %>% 
  filter(remission == "R" & crp < 5)

coldata_noMtx <- coldata_noMtx %>% 
  filter(remission == "R" & crp < 5)

coldata_mtx_noSyst <- coldata_mtx_noSyst %>% 
  filter(remission == "R" & crp < 5)

coldata_noPred <- coldata_noPred %>% 
  filter(remission == "R" & crp < 5)


# Pre processing ----
all_diagnosis <- data.frame()

diag_noBio <- coldata_noBio %>% 
  mutate(Comparison = "anti-TNF x no biologics",
         Treatment = antiTNF_vs_noBiologics,
         biologics_TNF = antiTNF_vs_noBiologics) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)

diag_nonTNF <- coldata_nonTNF %>% 
  mutate(Comparison = "anti-TNF x non-anti-TNF",
         Treatment = antiTNF_vs_nonAntiTNF,
         biologics_TNF = antiTNF_vs_nonAntiTNF ) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF) 


diag_noPred <- coldata_noPred %>% 
  mutate(Comparison = "pred x no pred",
         Treatment = pred_vs_noPred) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)

diag_pred_noSyst <- coldata_pred_noSyst %>% 
  mutate(Comparison = "pred x no syst",
         Treatment = pred_vs_noSyst) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)


diag_aza_noAza <- coldata_noAza %>% 
  mutate(Comparison = "aza x no aza",
         Treatment = aza_vs_noAza) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)

diag_aza_noSyst <- coldata_aza_noSyst %>% 
  mutate(Comparison = "aza x no syst",
         Treatment = aza_vs_noSyst) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)


diag_mtx_noMtx <- coldata_noMtx %>% 
  mutate(Comparison = "MTX x no MTX",
         Treatment = mtx_vs_noMtx) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)

diag_mtx_noSyst <- coldata_mtx_noSyst %>% 
  mutate(Comparison = "mtx x no syst",
         Treatment = mtx_vs_noSyst) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)

all_diagnosis <- rbind(diag_noBio, diag_nonTNF, diag_noPred, diag_pred_noSyst,
                       diag_aza_noAza,diag_aza_noSyst,diag_mtx_noMtx, diag_mtx_noSyst) 

all_diagnosis$Treatment <- as.factor(all_diagnosis$Treatment)


all_diagnosis<- all_diagnosis %>% 
  mutate(Treatment_binary = case_when(
    (Treatment == "anti_tnf" | Treatment == "pred" | Treatment == "aza" | Treatment == "mtx") ~ "Yes",
    (Treatment == "no_biologics" | Treatment == "non_anti_tnf" | Treatment == "no_mtx" |
       Treatment == "no_syst" | Treatment == "no_aza" | Treatment == "no_pred") ~ "No"
  ))


# Plot by diagnosis ----
breaks <- seq(from = 0, to = 600, by = 50)

plot1 <- ggplot(all_diagnosis, aes(Comparison, fill = diagnosis_class)) +
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
       title = "Patients included in comparisons") +
  scale_fill_brewer(type = "qual", palette = "Pastel2")
plot1 

ggsave(plot1, file = "Output_files/Patient_distribution/R_crp/byComparison_byDiagnosis.png", height = 8, width = 15, units = "in", dpi = 300)



##### Splitting by no syst/nobio and treatment

graphics.off()
rm(list = ls())

setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

library(RColorBrewer)
library(tidyverse)

# Loading files ----
#inactive disease patients
setwd("C:/Documents/Masters thesis/EZE_cohort/Output_files/Maaslin2/Tables/Remission_crp") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort/Output_files/Maaslin2/Tables/Remission_crp") #pc

files_R = list.files(pattern="*.txt")
coldata_R = lapply(files_R, read.delim)

coldata_noBio <- coldata_R[[4]]
coldata_aza_noSyst <- coldata_R[[5]]
coldata_pred_noSyst <- coldata_R[[6]]

# all patients
setwd("C:/Documents/Masters thesis/EZE_cohort/Output_files/Maaslin2/Tables") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort/Output_files/Maaslin2/Tables") #pc

files = list.files(pattern="*.txt")
coldata = lapply(files, read.delim)

coldata_nonTNF <- coldata[[3]]
coldata_noAza <- coldata[[4]]
coldata_noMtx <- coldata[[6]]
coldata_mtx_noSyst <- coldata[[7]]
coldata_noPred <- coldata[[8]]

# filtering to inactive + crp
setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

coldata_nonTNF <- coldata_nonTNF %>% 
  filter(remission == "R" & crp < 5)

coldata_noAza <- coldata_noAza %>% 
  filter(remission == "R" & crp < 5)

coldata_noMtx <- coldata_noMtx %>% 
  filter(remission == "R" & crp < 5)

coldata_mtx_noSyst <- coldata_mtx_noSyst %>% 
  filter(remission == "R" & crp < 5)

coldata_noPred <- coldata_noPred %>% 
  filter(remission == "R" & crp < 5)


# Pre processing ----
all_diagnosis <- data.frame()

diag_noBio <- coldata_noBio %>% 
  mutate(Comparison = "anti-TNF x no biologics",
         Treatment = antiTNF_vs_noBiologics,
         biologics_TNF = antiTNF_vs_noBiologics) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)

diag_nonTNF <- coldata_nonTNF %>% 
  mutate(Comparison = "anti-TNF x non-anti-TNF",
         Treatment = antiTNF_vs_nonAntiTNF,
         biologics_TNF = antiTNF_vs_nonAntiTNF ) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF) 


diag_noPred <- coldata_noPred %>% 
  mutate(Comparison = "pred x no pred",
         Treatment = pred_vs_noPred) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)

diag_pred_noSyst <- coldata_pred_noSyst %>% 
  mutate(Comparison = "pred x no syst",
         Treatment = pred_vs_noSyst) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)


diag_aza_noAza <- coldata_noAza %>% 
  mutate(Comparison = "aza x no aza",
         Treatment = aza_vs_noAza) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)

diag_aza_noSyst <- coldata_aza_noSyst %>% 
  mutate(Comparison = "aza x no syst",
         Treatment = aza_vs_noSyst) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)


diag_mtx_noMtx <- coldata_noMtx %>% 
  mutate(Comparison = "MTX x no MTX",
         Treatment = mtx_vs_noMtx) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)

diag_mtx_noSyst <- coldata_mtx_noSyst %>% 
  mutate(Comparison = "mtx x no syst",
         Treatment = mtx_vs_noSyst) %>% 
  select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)

all_diagnosis <- rbind(diag_noBio, diag_nonTNF, diag_noPred, diag_pred_noSyst,
                       diag_aza_noAza,diag_aza_noSyst,diag_mtx_noMtx, diag_mtx_noSyst) 

all_diagnosis$Treatment <- as.factor(all_diagnosis$Treatment)

# all_diagnosis<- all_diagnosis %>% 
#   mutate(Treatment_binary = case_when(
#     (Treatment == "anti_tnf" | Treatment == "pred" | Treatment == "aza" | Treatment == "mtx") ~ "Yes",
#     (Treatment == "no_biologics" | Treatment == "non_anti_tnf" | Treatment == "no_mtx" |
#        Treatment == "no_syst" | Treatment == "no_aza" | Treatment == "no_pred") ~ "No"
#   ))



diagnosis_treatment <- all_diagnosis %>% 
  filter(Treatment == "anti_tnf" | 
           Treatment == "pred" |
           Treatment == "aza" |
           Treatment == "mtx"
         )


diagnosis_no_treatment <- all_diagnosis %>% 
  filter(Treatment == "no_biologics" |
           Treatment == "non_anti_tnf" |
           Treatment == "no_syst" |
           Treatment == "no_pred" |
           Treatment == "no_aza" |
           Treatment == "no_mtx"
  )


# Plot by diagnosis ----
breaks <- seq(from = 0, to = 200, by = 25)

plot1 <- ggplot(diagnosis_treatment, aes(Comparison, fill = diagnosis_class)) +
  geom_bar() + theme_bw() +
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = 0.5, position = position_stack(vjust = 0.5), colour = "#808080", size = 5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.title.x=element_blank(),
    axis.title.y =element_blank(),
    title = element_text(size = 24),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 14)
  ) +
  scale_y_continuous(breaks = breaks, limits = c(0,150)) +
  labs(fill= 'Diagnosis',
       title = "Patients using treatment") +
  scale_fill_brewer(type = "qual", palette = "Pastel2")
plot1 

ggsave(plot1, file = "Output_files/Patient_distribution/R_crp/byComparison_byDiagnosis_treatment.pdf", height = 8, width = 15, units = "in", dpi = 300)


plot2 <- ggplot(diagnosis_no_treatment, aes(Comparison, fill = diagnosis_class)) +
  geom_bar() + theme_bw() +
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = 0.5, position = position_stack(vjust = 0.5), colour = "#808080", size = 5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.title.x=element_blank(),
    axis.title.y =element_blank(),
    title = element_text(size = 24),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 14)
  ) +
  scale_y_continuous(breaks = breaks, limits = c(0,150)) +
  labs(fill= 'Diagnosis',
       title = "Patients included not using treatments") +
  scale_fill_brewer(type = "qual", palette = "Pastel2")
plot2 

ggsave(plot2, file = "Output_files/Patient_distribution/R_crp/byComparison_byDiagnosis_notreatment.pdf", height = 8, width = 15, units = "in", dpi = 300)
