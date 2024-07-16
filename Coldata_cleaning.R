# EZE cohort: Therapy Signatures
  # Data cleaning, filtering out patients, renaming and recoding columns
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(tidyverse)
library(stringr)

folder <- "Cleaned_tables/"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# Loading data ----
data <- read.csv("Raw_tables/EZECohort_coldata_raw.csv")

# Remove patients that will not be included in analysis (non omics, outliers) ----
excluded_patients <- c("EZE323","EZE421", "EZE253", "EZE128", "EZE023", "EZE565", "EZE272", "EZE151", "EZE357")

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


data_filtered <- data_filtered %>% 
  mutate(biologics = ifelse(current_biologics == "None" | current_biologics == "New_therapy_after_sampling", "no_biologics", "biologics")) %>%
  relocate(biologics, .after = current_biologics)

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
  mutate(prednisolone_status = recode(prednisolone_status,
                                      "1" = "Never",
                                      "2" = "Former",
                                      "3" = "Current"))

data_filtered <- data_filtered %>% 
  mutate(remission = recode(remission,
                            "0" = "NR",
                            "1" = "R"))

# Create age groups ----
data_filtered <- data_filtered %>% 
  mutate(age_group = cut(age, breaks=c(0,20,35,50,65,Inf), include.lowest = FALSE, labels = c("0_20", "21_35", "36_50", "51_65", "65_inf"))) %>% 
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

# Clean special characters ----
data_filtered$diagnosis_class <- str_replace_all(data_filtered$diagnosis_class, ", ", "_")

# NAs ----
# Replace empty cells with NA
data_filtered$smoking <- str_replace_na(data_filtered$smoking, "NA")
data_filtered$remission <- str_replace_na(data_filtered$remission, "NA")
#data_filtered$bmi_class <- str_replace_na(data_filtered$bmi_class, "NA") 

# Exclude patients without crp and bmi data
data_filtered <- data_filtered %>% 
  filter(!(is.na(crp))) %>% 
  filter(!(is.na(bmi_class)))

# Exclude patients with more than one diagnosis
data_filtered <- data_filtered %>% 
  filter(diagnosis_class == "CD" | diagnosis_class == "PsA" |
           diagnosis_class == "Pso" |diagnosis_class == "RA" |
           diagnosis_class == "SLE" | diagnosis_class == "UC" |
           diagnosis_class == "Arthrosis")

# Log transformation of CRP levels ----
data_filtered <- data_filtered %>% 
  mutate(crp_log = log(crp + 0.001)) %>% 
  relocate(crp_log, .after = crp)

# Selecting data that will be used for models ----
coldata_maaslin <- data_filtered %>% 
  dplyr::select("study_id", "diagnosis_class", "age_group", "sex", "remission", "bmi_class", 
                "Azathioprin", "Prednisolon", "No_syst", "crp_log", "current_biologics", 
                "biologics_TNF", "biologics")

# Saving cleaned data tables ----
write.table(data_filtered, "Cleaned_tables/EZECohort_coldata_clean.txt", sep="\t", row.names = TRUE)
write.table(coldata_maaslin, "Cleaned_tables/EZECohort_coldata_models.txt", sep="\t", row.names = FALSE) 