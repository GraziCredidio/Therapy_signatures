# EZE cohort
  # Filtering out patients, renaming columns

graphics.off()
rm(list = ls())

setwd("C:\\Documents\\Masters thesis\\EZE_cohort") #laptop
setwd("D:\\Documentos\\Workspace\\Masters-Thesis\\EZE\\EZE_cohort") #PC

# Loading packages ----
library(tidyverse)
# loading data ----
data <- read.csv("Raw_tables\\EZECohort-GraziellaEZE_DATA_2023-01-23_1313.csv")

# removing patients (non omics, NA in biologics and comination of therapies not included) ----
excluded_patients <- c("EZE323","EZE421") # NA biologics
            
data_filtered <- data %>%
  filter(inclusion_omics == 1) %>% 
  filter(!((study_id %in% excluded_patients)))

# relabelling and recoding columns ----
data_filtered_labelled <- data_filtered %>%
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

data_filtered_labelled <- data_filtered_labelled %>% 
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


data_filtered_labelled <- data_filtered_labelled %>% 
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



data_filtered_labelled %>% 
  group_by(current_biologics) %>% 
  summarise(no_rows = length(current_biologics)) 

# Writing and loading processed data tables ----
write.csv(data_filtered_labelled, "Cleaned_tables/EZECohort_includedPatients_complete_labelled_09.02.csv", row.names = FALSE) #31.01: exclusion of combination therapy patients
write.table(data_filtered_labelled, "Cleaned_tables/EZECohort_includedPatients_complete_labelled_09.02.txt", sep="\t", row.names = TRUE)
