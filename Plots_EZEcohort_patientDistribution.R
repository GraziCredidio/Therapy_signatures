# EZE cohort
# Description of comparisons

graphics.off()
rm(list = ls())



library(RColorBrewer)
library(tidyverse)

# Loading files ----
coldata_all <- read.csv("Cleaned_tables/EZECohort_coldata_maaslin_complete_16.03.txt", sep="\t")

coldata_noBio_allpatients <- read.table("Output_files/Maaslin2/Tables/Without_pred_dose/coldata_maaslin2_antiTNF_vs_noBiologics_28.02.txt", sep = "\t")

setwd("C:/Documents/Masters thesis/EZE_cohort/Output_files/Maaslin2/Tables") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort/Output_files/Maaslin2/Tables") #pc

files = list.files(pattern="*.txt")
coldata = lapply(files, read.delim)

coldata_noBio <- coldata[[2]]
coldata_nonTNF <- coldata[[3]]

coldata_noAza <- coldata[[4]]
coldata_aza_noSyst <- coldata[[5]]

coldata_noMtx <- coldata[[6]]
coldata_mtx_noSyst <- coldata[[7]]

coldata_noPred <- coldata[[8]]
coldata_pred_noSyst <- coldata[[9]]

# coldata exploration ----
coldata %>% 
  group_by(bmi_class) %>% 
  summarise(no_rows = length(bmi_class))

table(coldata$age_group2)

coldata %>% 
  group_by(biologics_TNF) %>% 
  summarise(no_rows = length(biologics_TNF)) # 24 non_anti_tnf

table(coldata$biologics)

table(coldata$biologics_TNF)

table(coldata$diagnosis_class) 

table(coldata$remission)

table(coldata$treatment_intensification)

table(coldata$Prednisolon)

table(coldata$MTX)

table(coldata$Azathioprin)

table(coldata$No_syst)

nrow(coldata[coldata$Prednisolon == 1 & coldata$Azathioprin == 1, ]) #16 use pred and aza 
nrow(coldata[coldata$Prednisolon == 1 & coldata$MTX == 1, ]) #57 use pred and mtx

table(coldata_pred_noSyst$remission)
table(coldata_aza_noSyst$remission)
table(coldata_mtx_noSyst$remission)

coldata_nonTNF %>% 
  filter(remission == "R" & antiTNF_vs_nonAntiTNF == "non_anti_tnf") %>% 
  nrow()

coldata_noBio %>% 
  filter(remission == "R" & antiTNF_vs_noBiologics == "no_biologics") %>% 
  nrow()

table(coldata$prednisolone_dose)

coldata_pred_noSyst %>% 
  group_by(diagnosis_class) %>% 
  filter(pred_vs_noSyst == "pred") %>% 
  summarise(nrow = length(diagnosis_class))

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

# Plots
setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
setwd("~/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC

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

ggsave(plot1, file = "Output_files/Patient_distribution/byComparison_byDiagnosis.png", height = 8, width = 15, units = "in", dpi = 300)

# Plot by treatment ----
plot2 <- ggplot(all_diagnosis, aes(Comparison, fill = Treatment_binary)) +
  geom_bar() + theme_bw() +
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = 0.5, vjust = 1.5, position ="stack", colour = "white", size = 5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    title = element_text(size = 24),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  ) +
  scale_y_continuous(breaks = breaks) +
  labs(fill= 'Treatment',
       title = "Patients included in comparisons by treatment") +
  scale_fill_brewer(type = "qual", palette = "Pastel1")
plot2

ggsave(plot2, file = "Output_files/Patient_distribution/byComparison_byTreatment.png", height = 8, width = 15, units = "in", dpi = 300)

# plot by age group ----
plot3 <- ggplot(all_diagnosis, aes(Comparison, fill = age_group2)) +
  geom_bar() + theme_bw() +
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = 0.5, position = position_stack(vjust = 0.5), colour = "#808080", size = 5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    title = element_text(size = 24),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  ) +
  scale_y_continuous(breaks = breaks) +
  labs(fill= 'Age group',
       title = "Patients included in comparisons by age group") +
  scale_fill_brewer(type = "qual", palette = "Pastel1")
plot3

ggsave(plot3, file = "Output_files/Patient_distribution/byComparison_byAge.png", height = 8, width = 15, units = "in", dpi = 300)

# Plot comparisons x biologics ----
cols = c("anti_tnf" = "#a6cee3", "non_anti_tnf" = "#1f78b4",
         "no_biologics" = "#b2df8a")

plot4 <- ggplot(all_diagnosis, aes(Comparison, fill = biologics_TNF)) +
  geom_bar() + theme_bw() +
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = 0.5, position = position_stack(vjust = 0.5), colour = "white", size = 5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    title = element_text(size = 24),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 14)
  ) +
  scale_y_continuous(breaks = breaks) +
  labs(fill= 'Use of biologics',
       title = "Patients included in comparisons by biologics treatment") +
  scale_fill_manual(values = cols)
plot4

ggsave(plot4, file = "Output_files/Patient_distribution/byComparison_byBiologics.png", height = 8, width = 15, units = "in", dpi = 300)


# plot of diagnosis x biologics ----
breaks <- seq(from = 0, to = 200, by = 30)

plot5 <- ggplot(coldata_all, aes(diagnosis_class, fill = biologics_TNF)) +
  geom_bar() + theme_bw() +
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = 0.5, position = position_stack(vjust = 0.5), colour = "white", size = 5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    title = element_text(size = 24),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 14)
  ) +
  scale_y_continuous(breaks = breaks) +
  labs(fill= 'Use of biologics',
       title = "Patient distribution by diagnosis and biologics treatment") +
  scale_fill_manual(values = cols)
plot5

ggsave(plot5, file = "Output_files/Patient_distribution/byDiagnosis_byBiologics.png", height = 8, width = 15, units = "in", dpi = 300)

# Plots using systemic therapies ----
## Exploring the number of patients that use aza + pred, mtx + pred:

## preprocessing
coldata_syst <- coldata_all %>% 
  mutate(MTX = recode(MTX,
                      "0" = "No",
                      "1" = "Yes"),
         Prednisolon = recode(Prednisolon,
                              "0" = "No",
                              "1" = "Yes"),
         Azathioprin = recode(Azathioprin,
                              "0" = "No",
                              "1" = "Yes")) %>% 
  mutate(across(c(diagnosis_class, MTX, Prednisolon, Azathioprin), as.factor))

plot6 <- function(syst_therapy, title, path){
  plot <- ggplot(coldata_syst, aes(diagnosis_class, fill = syst_therapy)) +
  geom_bar() + theme_bw() +
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = 0.5, position = position_stack(vjust = 0.5), colour = "#808080", size = 5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    title = element_text(size = 30),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 20)
  ) +
  scale_y_continuous(breaks = breaks) +
  labs(fill= title,
       title = paste("Diagnosis by", title, sep = " ")) +
  scale_fill_brewer(palette = "Pastel2")
  ggsave(plot, file = path, height = 8, width = 15, units = "in", dpi = 300)
  
plot
}

plot6(coldata_syst$MTX, "MTX", "Output_files/Patient_distribution/byDiagnosis_byMTX.png")
plot6(coldata_syst$Prednisolon, "Prednisolon", "Output_files/Patient_distribution/byDiagnosis_byPred.png")
plot6(coldata_syst$Azathioprin, "Azathioprin", "Output_files/Patient_distribution/byDiagnosis_byAza.png")



# Plot anti_tnf x no biologics by use of biologics (with all patients) ----
# diag_noBio_all_patients <- coldata_noBio_allpatients %>% 
#   mutate(Comparison = "anti-TNF x no biologics",
#          Treatment = antiTNF_vs_noBiologics,
#          biologics_TNF = antiTNF_vs_noBiologics) %>% 
#   select(study_id, diagnosis_class, Treatment, Comparison, age_group2, biologics_TNF)
# 

breaks <- seq(from = 0, to = 200, by = 30)

plot7 <- ggplot(diag_noBio, aes(diagnosis_class, fill = biologics_TNF)) +
  geom_bar() + theme_bw() +
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = 0.5, vjust = 2, position ="stack", colour = "white", size = 5) +
  theme(
    axis.text.x = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    title = element_text(size = 24),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 14)
  ) +
  scale_y_continuous(breaks = breaks) +
  labs(fill= 'Use of anti-TNF',
       title = "Patients included in comparison anti-TNF x no Biologics") +
  scale_fill_manual(values = cols)
plot7

ggsave(plot7, file = "Output_files/Patient_distribution/byantiTNF_byDiagnosis_afterFiltering.png", height = 8, width = 15, units = "in", dpi = 300)





### Medication use distribution accross comparisons ----
## How many patients using Aza also use other systemic therapies? (Mtx and pred)
# Exclude no syst patients
coldata_aza <- coldata_aza_noSyst %>% 
  filter(!(No_syst == "1"))


coldata_aza %>% 
  dplyr::select(MTX, Prednisolon) %>% 
  filter(MTX == 1 | Prednisolon == 1) %>% 
  tibble() #out of 51 patients using aza, 14 also use pred


## How many patients using Pred also use other systemic therapies? (Mtx and aza)
# Exclude no syst patients
coldata_pred <- coldata_pred_noSyst %>% 
  filter(!(No_syst == "1"))

coldata_pred %>% 
  dplyr::select(MTX, Azathioprin) %>% 
  filter(MTX == 1 | Azathioprin == 1) %>% 
  tibble() %>% 
  print(n = 80)#out of 192 patients using pred, 73 use either mtx (57) or aza (16) too


## How many patients using MTX also use other systemic therapies? (pred and aza)
# Exclude no syst patients
coldata_mtx <- coldata_mtx_noSyst %>% 
  filter(!(No_syst == "1"))

coldata_mtx %>% 
  dplyr::select(Prednisolon, Azathioprin) %>% 
  filter(Prednisolon == 1 | Azathioprin == 1) %>% 
  tibble() %>% 
  print(n = 80) #out of 95 patients using MTX, 49 also use pred 