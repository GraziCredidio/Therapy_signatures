# EZE cohort
  # Data wrangling - organizing and recategorizing covariates

graphics.off()
rm(list = ls())

setwd("C:\\Documents\\Masters thesis\\EZE_cohort") #laptop
setwd("D:\\Documentos\\Workspace\\Masters-Thesis\\EZE\\EZE_cohort") #PC

# Loading packages ----
library(tidyverse)
library(stringr)
library(rebus)

# loading data ----
coldata <- read.csv("Cleaned_tables\\EZECohort_includedPatients_complete_labelled_09.02.txt", sep = "\t")

# Selecting columns ----
data_selected <- coldata %>% 
  select("study_id","diagnosis_freetext_cleaned", "age", "sex", "remission", 
         "bmi", "smoking", "Azathioprin", "MTX", "Prednisolon", "Leflunomid", "Mycophenolat.Mofetil", "No_syst", 
         "crp", "current_biologics","comorb_hypertension","comorb_diabetes") %>% 
  dplyr::rename(diagnosis_class = "diagnosis_freetext_cleaned")

# Recoding, classifying and replacing NAs----
## Sex
data_selected <- data_selected %>% 
  mutate(sex = recode(sex,
                      "1" = "Male",
                      "2" = "Female"))

## Smoking
data_selected <- data_selected %>% 
  mutate(smoking = recode(smoking,
                          "0" = "Never",
                          "1" = "Ex",
                          "2" = "Current"))
#which(is.na(data_selected$smoking), arr.ind=TRUE) # 12 NA patients: EZE006, 007, 013, 069, 133, 147, 184, 186, 207, 499, 550
data_selected$smoking <- str_replace_na(data_selected$smoking, "NA")

## Remission
data_selected <- data_selected %>% 
  mutate(remission = recode(remission,
                            "0" = "Non_remitter",
                            "1" = "Remitter"))

#which(is.na(data_selected$remission), arr.ind=TRUE) # 5 NA patients: EZE208, 234, 253, 402, 414
data_selected$remission <- str_replace_na(data_selected$remission, "NA")

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

#which(is.na(data_selected$bmi), arr.ind=TRUE) # 7 NA patients: EZE38, 40, 81, 416, 427, 522, 635
data_selected$bmi_class <- str_replace_na(data_selected$bmi_class, "NA") 

data_selected %>% 
  group_by(bmi_class) %>% 
  summarise(no_rows = length(bmi_class))

## Age groups
data_selected <- data_selected %>% 
  mutate(age_group = cut(age, breaks=c(0,25,65,Inf), include.lowest = FALSE, labels = c("Young", "Adult", "Senior"))) %>% 
  relocate(age_group, .after = age) 
   
data_selected %>%
  group_by(age_group) %>%
  summarise(no_rows = length(age_group))

## Biologics (anti-tnf x non anti-tnf)
# anti_tnf <- c("Infliximab", "Golimumab", "Etanercept", "Certolizumab_pegol", "Adalimumab")
# non_anti_tnf <- c("Ustekinumab", "Vedolizumab", "Belimumab", "Rituximab", 
#                   "Secukinumab", "Tocilizumab", "Abatacept", "Anakinra")

data_selected <- data_selected %>% 
  mutate(biologics = ifelse(current_biologics == "None" | current_biologics == "New_therapy_after_sampling", "no_biologics", "biologics")) %>%
  relocate(biologics, .after = current_biologics) 

table(data_selected$biologics)

## Diagnosis 
# reordering names combinations
data_selected$diagnosis_class[data_selected$diagnosis_class == "Pso, CD"] <- "CD, Pso"
data_selected$diagnosis_class[data_selected$diagnosis_class == "RA, CD"] <- "CD, RA"
data_selected$diagnosis_class[data_selected$diagnosis_class == "Pso, RA"] <- "RA, Pso"
data_selected$diagnosis_class[data_selected$diagnosis_class == "SLE, RA"] <- "RA, SLE"

# renaming other diagnosis as other
data_selected$diagnosis_class[data_selected$diagnosis_class == "Palindromic rheumatism" |
                                data_selected$diagnosis_class == "Still's disease" |
                                data_selected$diagnosis_class == "Still's disease, Pso" |
                                data_selected$diagnosis_class == "TRAPS" |
                                data_selected$diagnosis_class == "SpA" |
                                data_selected$diagnosis_class == "SpA, Pso" |
                                data_selected$diagnosis_class == "GPA" |
                                data_selected$diagnosis_class == "GCA" |
                                data_selected$diagnosis_class == "Celiac disease" |
                                data_selected$diagnosis_class == "ANCA-associated vasculitis"|
                                data_selected$diagnosis_class == "RA, SpA"] <- "Other"

data_selected$diagnosis_class <- str_replace_all(data_selected$diagnosis_class, ", ", "_")

data_selected %>% 
  group_by(diagnosis_class) %>% 
  summarise(no_rows = length(diagnosis_class)) %>% 
  print(n =13)

## Comorbidities
#replacing NA value by zero (patient 1)
data_selected$comorb_hypertension[1] <- 0
data_selected$comorb_diabetes[1] <- 0


# Creating coldata ----
coldata_EZE <- data_selected %>% 
  select("study_id", "diagnosis_class", "age_group", "sex", "remission", "bmi_class",
         "smoking", "Azathioprin":"crp", "biologics":"comorb_diabetes")

write.csv(coldata_EZE, "Cleaned_tables/EZECohort_coldata_09.02.csv", row.names = FALSE) #08.02: new age group; 09.02: new coldata coding
