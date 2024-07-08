# EZE cohort
# Metadata table for Maaslin2: 
# Comparisons : Anti-tnf x no biologics; Anti-tnf x non-anti-tnf; 
#              systemic therapy x no systemic therapy; systemic therapy x no syst + other syst

graphics.off()
rm(list = ls())

#setwd("C:\\Documents/Masters thesis/EZE_cohort") #laptop
#setwd("D:\\Documentos/Workspace/Masters-Thesis/EZE/EZE_cohort") #PC


# Loading packages ----
library(tidyverse)

# Loading files ----
coldata <- read.table("Cleaned_tables/EZECohort_includedPatients_complete_labelled_09.02.txt", sep = "\t")
names(coldata)
# selecting columns
data_selected <- coldata %>% 
  dplyr::select("study_id", "CD":"Other", "diagnosis_freetext_cleaned", "age", "sex",
         "bmi", "Azathioprin":"No_syst", "remission", 
         "prednisolone_status","prednisolone_dose", "prednisolone_duration", "prednisolone_start", 
         "crp", "Infliximab" : "None_above", "current_biologics",
         "smoking", "comorb_hypertension", "comorb_diabetes", "treatment_intensification") %>% 
  dplyr::rename(diagnosis_class = "diagnosis_freetext_cleaned")

# Recoding ----
## Sex
data_selected <- data_selected %>% 
  mutate(sex = recode(sex,
                      "1" = "Male",
                      "2" = "Female"))

## BMI
data_selected <- data_selected %>% 
  mutate(bmi_class = as.factor(case_when(
    bmi < 18 ~ "Malnutrition",
    bmi >= 18 & bmi < 25 ~ "Healthy",
    bmi >= 25 & bmi < 30 ~ "Overweight",
    bmi >= 30 & bmi < 35 ~ "Obesity_I",
    bmi >= 35 & bmi < 40 ~ "Obesity_II",
    bmi >= 40 ~ "Obesity_III"
  ))) %>%
  relocate(bmi_class, .after = bmi) 

data_selected %>% 
  group_by(bmi_class) %>% 
  summarise(no_rows = length(bmi_class))

## Age groups
data_selected <- data_selected %>% 
  mutate(age_group = cut(age, breaks=c(0,25,65,Inf), include.lowest = FALSE, labels = c("Young", "Adult", "Senior"))) %>% 
  relocate(age_group, .after = age)

data_selected <- data_selected %>% 
  mutate(age_group2 = cut(age, breaks=c(0,20,35,50,65,Inf), include.lowest = FALSE, labels = c("0_20", "21_35", "36_50", "51_65", "65_inf"))) %>% 
  relocate(age_group2, .after = age_group)

table(data_selected$age_group)
table(data_selected$age_group2)

## Prednisolon status
data_selected <- data_selected %>% 
  mutate(prednisolone_status = recode(prednisolone_status,
                      "1" = "Never",
                      "2" = "Former",
                      "3" = "Current"))

table(data_selected$prednisolone_status)

## Biologics (anti-tnf, non anti-tnf, no biologics)
# anti_tnf <- c("Infliximab", "Golimumab", "Etanercept", "Certolizumab_pegol", "Adalimumab")
# non_anti_tnf <- c("Ustekinumab", "Vedolizumab", "Belimumab", "Rituximab", 
#                   "Secukinumab", "Tocilizumab", "Abatacept", "Anakinra")

data_selected <- data_selected %>% 
  mutate(biologics_TNF = case_when(
    (current_biologics == "Infliximab" |
       current_biologics == "Golimumab" |
       current_biologics == "Etanercept" |
       current_biologics == "Certolizumab_pegol" |
       current_biologics == "Adalimumab") ~ "anti_tnf",
    
    (current_biologics == "Ustekinumab"|
       current_biologics == "Vedolizumab" |
       current_biologics == "Belimumab" |
       current_biologics == "Rituximab" |
       current_biologics == "Secukinumab" |
       current_biologics == "Tocilizumab" |
       current_biologics == "Abatacept" |
       current_biologics == "Anakinra") ~"non_anti_tnf",
    
    (current_biologics == "None" | current_biologics == "New_therapy_after_sampling") ~ "no_biologics"
  )) %>%
  relocate(biologics_TNF, .after = current_biologics) 

data_selected %>% 
  group_by(biologics_TNF) %>% 
  summarise(no_rows = length(biologics_TNF)) # only 33 non_anti_tnf

## Biologics (biologics, no biologics)
data_selected <- data_selected %>% 
  mutate(biologics = ifelse(current_biologics == "None" | 
                              current_biologics == "New_therapy_after_sampling", "no_biologics", "biologics")) %>%
  relocate(biologics, .after = current_biologics) 

table(data_selected$biologics)

## Diagnosis
# data_selected$diagnosis_class[data_selected$diagnosis_class == "Pso, CD"] <- "CD, Pso"
# data_selected$diagnosis_class[data_selected$diagnosis_class == "RA, CD"] <- "CD, RA"
# data_selected$diagnosis_class[data_selected$diagnosis_class == "Pso, RA"] <- "RA, Pso"
# data_selected$diagnosis_class[data_selected$diagnosis_class == "SLE, RA"] <- "RA, SLE"
# 
# # renaming other diagnosis as other
# data_selected$diagnosis_class[data_selected$diagnosis_class == "Palindromic rheumatism" |
#                                 data_selected$diagnosis_class == "Still's disease" |
#                                 data_selected$diagnosis_class == "Still's disease, Pso" |
#                                 data_selected$diagnosis_class == "TRAPS" |
#                                 data_selected$diagnosis_class == "SpA" |
#                                 data_selected$diagnosis_class == "SpA, Pso" |
#                                 data_selected$diagnosis_class == "GPA" |
#                                 data_selected$diagnosis_class == "GCA" |
#                                 data_selected$diagnosis_class == "Celiac disease" |
#                                 data_selected$diagnosis_class == "ANCA-associated vasculitis"|
#                                 data_selected$diagnosis_class == "RA, SpA"] <- "Other"
# 
# data_selected$diagnosis_class <- str_replace_all(data_selected$diagnosis_class, ", ", "_")
# 
# data_selected %>% 
#   group_by(diagnosis_class) %>% 
#   summarise(no_rows = length(diagnosis_class)) %>% 
#   print(n=23)

# excluding patients with more than 1 diagnosis
data_selected <- data_selected %>% 
  filter(diagnosis_class == "CD" | diagnosis_class == "PsA" |
           diagnosis_class == "Pso" |diagnosis_class == "RA" |
           diagnosis_class == "SLE" | diagnosis_class == "UC" |
           diagnosis_class == "Arthrosis")
table(data_selected$diagnosis_class)

## Remission
data_selected <- data_selected %>% 
  mutate(remission = recode(remission,
                            "0" = "NR",
                            "1" = "R"))
table(data_selected$remission)

# Excluding PCA outliers, NA crp and BMI ----
pca_outliers <- c("EZE253", "EZE128", "EZE023", "EZE565", "EZE272", "EZE151", "EZE357")

data_selected <- data_selected %>% 
  filter(!(study_id %in% pca_outliers))

data_selected <- data_selected %>% 
  filter(!(is.na(crp))) %>% 
  filter(!(is.na(bmi_class))) 

# CRP transformation ----
data_selected <- data_selected %>% 
  mutate(crp_log = log(crp + 0.001)) %>% 
  relocate(crp_log, .after = crp)
  

# Metada maaslin ----
coldata_maaslin <- data_selected %>% 
  dplyr::select("study_id", "diagnosis_class", "age_group", "age_group2", "age", "sex", 
         "remission", "bmi_class", "bmi", "Azathioprin":"crp_log", "current_biologics", 
         "biologics_TNF", "biologics", "smoking", "comorb_hypertension", "comorb_diabetes",
         "prednisolone_dose", "prednisolone_duration", "prednisolone_start", "treatment_intensification")

#apply(coldata_maaslin[,c(1:6,8:29)], 2, function(x) any(is.na(x))) #checking if there are any NAs in columns used in modeling

write.table(coldata_maaslin, "Cleaned_tables/EZECohort_coldata_maaslin_complete_16.03.txt", sep="\t", row.names = FALSE) 
write.csv(coldata_maaslin, "Cleaned_tables/EZECohort_coldata_maaslin_complete_16.03.csv", row.names = FALSE) 


# Checking correlation of age and crp ----
#install.packages("ggpubr")
library("ggpubr")

coldata <- read.csv("Cleaned_tables/EZECohort_coldata_maaslin_complete_16.03.csv", sep="\t")

## all patients
ggscatter(coldata, x = "age", y = "crp_log", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "age", ylab = "log crp")

shapiro.test(coldata$age) #not normal
shapiro.test(coldata$crp) # not normal

cor.test(coldata$age,coldata$crp, method = "spearman", exact = FALSE) #not correlated



## IBD
coldata_ibd <- coldata %>% 
  filter(diagnosis_class == "UC" | diagnosis_class == "CD") %>% 
  filter(!(age == "91")) #age outlier?
max(coldata_ibd$age)

ggscatter(coldata_ibd, x = "age", y = "crp_log", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "age", ylab = "log crp") #negative association?

shapiro.test(coldata_ibd$age) #not normal
shapiro.test(coldata_ibd$crp) # not normal

cor.test(coldata_ibd$age,coldata_ibd$crp, method = "spearman", exact = FALSE) #p-value = 0.003594

## RA and PsA
coldata_ra <- coldata %>% 
  filter(diagnosis_class == "RA" | diagnosis_class == "PsA")

ggscatter(coldata_ra, x = "age", y = "crp_log", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "age", ylab = "log crp") #weakly correlated

shapiro.test(coldata_ra$age) #not normal
shapiro.test(coldata_ra$crp) # not normal

cor.test(coldata_ra$age,coldata_ra$crp, method = "spearman", exact = FALSE) #p-value = 0.007429


## SLE
coldata_sle <- coldata %>% 
  filter(diagnosis_class == "SLE")

ggscatter(coldata_sle, x = "age", y = "crp_log", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "age", ylab = "log crp") #not correlated

shapiro.test(coldata_sle$age) # normal
shapiro.test(coldata_sle$crp) # not normal

cor.test(coldata_sle$age,coldata_sle$crp, method = "pearson", exact = FALSE) #p-value = 0.4531


## Pso
coldata_pso <- coldata %>% 
  filter(diagnosis_class == "Pso")

ggscatter(coldata_pso, x = "age", y = "crp_log", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "age", ylab = "log crp") #not correlated

shapiro.test(coldata_pso$age) # normal
shapiro.test(coldata_pso$crp) # not normal

cor.test(coldata_sle$age,coldata_sle$crp, method = "pearson", exact = FALSE) #p-value = 0.4531
