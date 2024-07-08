# EZE cohort: Therapy Signatures
  # Data cleaning, filtering out patients, renaming and recoding columns
  # Author: Graziella Credidio


# Loading packages ----
library(tidyverse)
library(stringr)
library(rebus)

# Loading data ----
data <- read.csv("./Raw_tables/EZECohort_coldata_raw.csv")

# Remove patients that will not be included in analysis (non omics, NA in biologics and combination of therapies) ----
excluded_patients <- c("EZE323","EZE421")
            
data_filtered <- data %>%
  filter(inclusion_omics == 1) %>% 
  filter(!((study_id %in% excluded_patients)))

# Recode columns ----
data_filtered <- data_filtered %>%
  dplyr::rename("CD" = "diagnosis___1",
         "UC" = "diagnosis___2",
         "RA" = "diagnosis___4",
         "PsA" = "diagnosis___5",
         "Pso" = "diagnosis___6",
         "SLE" = "diagnosis___9",
         "Arthrosis" = "diagnosis___10",
         "Other" = "diagnosis___99",
         "Azathioprin" = "current_system_therapy___1",
         "MTX" = "current_system_therapy___2",
         "6-MP" = "current_system_therapy___3",
         "Sulfasalazin" = "current_system_therapy___4",
         "Chloroquin" = "current_system_therapy___5",
         "Hydroxychloroquin" = "current_system_therapy___6",
         "Prednisolon" = "current_system_therapy___7",
         "Leflunomid" = "current_system_therapy___8",
         "Ciclosporin" = "current_system_therapy___9",
         "Mycophenolat-Mofetil" = "current_system_therapy___10",
         "Cyclophosphamid" = "current_system_therapy___11",
         "Tacrolimus" = "current_system_therapy___12",
         "Apremilast" = "current_system_therapy___13",
         "Mepacrine" = "current_system_therapy___14",
         "No_syst" = "current_system_therapy___99")

data_filtered <- data_filtered %>% 
  mutate(Infliximab = as.integer(current_biologics == 1),
         Adalimumab = as.integer(current_biologics == 2),
         Certolizumab_pegol = as.integer(current_biologics == 3),
         Golimumab = as.integer(current_biologics == 4),
         Etanercept = as.integer(current_biologics == 5),
         Vedolizumab = as.integer(current_biologics == 6),
         Ustekinumab = as.integer(current_biologics == 7),
         Risankizumab = as.integer(current_biologics == 8),
         Guselkumab = as.integer(current_biologics == 9),
         Ixekizumab = as.integer(current_biologics == 10),
         Anakinra = as.integer(current_biologics == 12),
         Canakinumab = as.integer(current_biologics == 13),
         Abatacept = as.integer(current_biologics == 14),
         Tocilizumab = as.integer(current_biologics == 15),
         Tofacitinib = as.integer(current_biologics == 16),
         Upadacitinib = as.integer(current_biologics == 17),
         Baricitinib = as.integer(current_biologics == 18),
         Filgotinib = as.integer(current_biologics == 19),
         Rituximab = as.integer(current_biologics == 20),
         Olamkizept = as.integer(current_biologics == 21),
         Denosumab = as.integer(current_biologics == 22),
         Brodalumab = as.integer(current_biologics == 23),
         Secukinumab = as.integer(current_biologics == 24),
         Canakinumab = as.integer(current_biologics == 25),
         Belimumab = as.integer(current_biologics == 26),
         New_therapy_after_sampling = as.integer(current_biologics == 98),
         None_above = as.integer(current_biologics == 99))

data_filtered <- data_filtered %>% 
  mutate(current_biologics = recode(current_biologics,
                                    "1" = "Infliximab",
                                    "2" = "Adalimumab",
                                    "3" = "Certolizumab_pegol",
                                    "4" = "Golimumab",
                                    "5" = "Etanercept",
                                    "6" = "Vedolizumab",
                                    "7" = "Ustekinumab" ,
                                    "8" = "Risankizumab",
                                    "9" = "Guselkumab",
                                    "10" = "Ixekizumab",
                                    "12" = "Anakinra",
                                    "13" = "Canakinumab" ,
                                    "14" = "Abatacept",
                                    "15" = "Tocilizumab",
                                    "16" = "Tofacitinib",
                                    "17" = "Upadacitinib",
                                    "18" = "Baricitinib",
                                    "19" = "Filgotinib",
                                    "20" = "Rituximab",
                                    "21" = "Olamkizept",
                                    "22" = "Denosumab",
                                    "23" = "Brodalumab",
                                    "24" = "Secukinumab",
                                    "25" = "Canakinumab",
                                    "26" = "Belimumab",
                                    "98" = "New_therapy_after_sampling",
                                    "99" = "None"))

data_filtered <- data_filtered %>% 
  mutate(sex = recode(sex,
                      "1" = "Male",
                      "2" = "Female"))

data_filtered <- data_filtered %>% 
  mutate(smoking = recode(smoking,
                          "0" = "Never",
                          "1" = "Ex",
                          "2" = "Current"))

data_filtered <- data_filtered %>% 
  mutate(bmi_class = as.factor(case_when(
    bmi < 18 ~ "Malnutrition",
    bmi >= 18 & bmi < 25 ~ "Healthy",
    bmi >= 25 & bmi < 30 ~ "Overweight",
    bmi >= 30 & bmi < 35 ~ "Obesity_I",
    bmi >= 35 & bmi < 40 ~ "Obesity_II",
    bmi >= 40 ~ "Obesity_III"
  ))) %>%
  relocate(bmi_class, .after = bmi) 

data_filtered <- data_filtered %>% 
  mutate(biologics = ifelse(current_biologics == "None" | current_biologics == "New_therapy_after_sampling", "no_biologics", "biologics")) %>%
  relocate(biologics, .after = current_biologics)

# Create age groups ----
data_filtered <- data_filtered %>% 
  mutate(age_group = cut(age, breaks=c(0,25,65,Inf), include.lowest = FALSE, labels = c("Young", "Adult", "Senior"))) %>% 
  relocate(age_group, .after = age) 

# Cleaning diagnosis columns ----
# Rename diagnosis column
data_filtered <- data_filtered %>% 
  dplyr::rename(diagnosis_class = "diagnosis_freetext_cleaned")

# Reorder diagnosis combinations to match
data_filtered$diagnosis_class[data_filtered$diagnosis_class == "Pso, CD"] <- "CD, Pso"
data_filtered$diagnosis_class[data_filtered$diagnosis_class == "RA, CD"] <- "CD, RA"
data_filtered$diagnosis_class[data_filtered$diagnosis_class == "Pso, RA"] <- "RA, Pso"
data_filtered$diagnosis_class[data_filtered$diagnosis_class == "SLE, RA"] <- "RA, SLE"

# Rename diagnosis as "other"
data_filtered$diagnosis_class[data_filtered$diagnosis_class == "Palindromic rheumatism" |
                                data_filtered$diagnosis_class == "Still's disease" |
                                data_filtered$diagnosis_class == "Still's disease, Pso" |
                                data_filtered$diagnosis_class == "TRAPS" |
                                data_filtered$diagnosis_class == "SpA" |
                                data_filtered$diagnosis_class == "SpA, Pso" |
                                data_filtered$diagnosis_class == "GPA" |
                                data_filtered$diagnosis_class == "GCA" |
                                data_filtered$diagnosis_class == "Celiac disease" |
                                data_filtered$diagnosis_class == "ANCA-associated vasculitis"|
                                data_filtered$diagnosis_class == "RA, SpA"] <- "Other"
# cleaning special characters
data_filtered$diagnosis_class <- str_replace_all(data_filtered$diagnosis_class, ", ", "_")

# Replacing NAs ----
data_filtered$smoking <- str_replace_na(data_filtered$smoking, "NA")
data_filtered$remission <- str_replace_na(data_filtered$remission, "NA")
data_filtered$bmi_class <- str_replace_na(data_filtered$bmi_class, "NA") 


# Selecting data to be saved----
coldata_EZE <- data_filtered %>% 
  select("study_id", "diagnosis_class", "age_group", "sex", "remission", "bmi_class",
         "smoking", "Azathioprin":"crp", "biologics":"comorb_diabetes")

# Writing and loading processed data tables ----
write.table(coldata_EZE, "Cleaned_tables/EZECohort_coldata_clean.txt", sep="\t", row.names = TRUE)