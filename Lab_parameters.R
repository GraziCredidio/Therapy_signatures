# EZE cohort
# Boxplot with lab values - r_crp therapy patients

rm(list = ls())

library(ggpubr)

coldata_R_antiTNF <- read.table("Output_files/Maaslin2/Tables/Remission/coldata_maaslin2_antiTNF_vs_noBio_R_17.03.txt", sep = "\t")
redcap <- read.csv("Raw_tables/EZECohort-GraziellaEZE_DATA_2023-04-12_1328.csv")

# Extracting leukocytes, thrombocytes and erythrocytes from redcap coldata (updated)
redcap_patients <- redcap[redcap$study_id %in% coldata_R_antiTNF$study_id,]
redcap_patients <- redcap_patients %>% 
  dplyr::select(study_id, leucocytes, erythrocytes, thrombocytes)

idx <- match(coldata_R_antiTNF$study_id, redcap_patients$study_id) #matching order of patients
ord_redcap_patients <- redcap_patients[idx,]

coldata_R_antiTNF <- left_join(coldata_R_antiTNF, ord_redcap_patients, by = "study_id")



male <- coldata_R_antiTNF %>% 
  filter(sex == "Male")

mean(male$thrombocytes)
sd(male$thrombocytes)

mean(male$erythrocytes)
sd(male$erythrocytes)

mean(male$leucocytes)
sd(male$leucocytes)

female <- coldata_R_antiTNF %>% 
  filter(sex == "Female") %>% 
  filter(!(is.na(thrombocytes)))

mean(female$thrombocytes)
sd(female$thrombocytes)

mean(female$erythrocytes)
sd(female$erythrocytes)

mean(female$leucocytes)
sd(female$leucocytes)


write.table(male, "Output_files/Lab_parameters/Males_antiTNF.txt", sep = "\t")
write.table(female, "Output_files/Lab_parameters/Females_antiTNF.txt", sep = "\t")



# Loading data of other therapies
male_antiTNF <- read.table("Output_files/Lab_parameters/Males_antiTNF.txt", sep = "\t")
male_antiTNF <- male_antiTNF[,c("study_id", "thrombocytes", "erythrocytes", "leucocytes")]
male_antiTNF$analysis <- "antiTNF"
female_antiTNF <- read.table("Output_files/Lab_parameters/Females_antiTNF.txt", sep = "\t")
female_antiTNF <- female_antiTNF[,c("study_id", "thrombocytes", "erythrocytes", "leucocytes")]
female_antiTNF$analysis <- "antiTNF"

male_aza <- read.table("Output_files/Lab_parameters/Males_aza.txt", sep = "\t")
male_aza <- male_aza[,c("study_id", "thrombocytes", "erythrocytes", "leucocytes")]
male_aza$analysis <- "aza"
female_aza <- read.table("Output_files/Lab_parameters/Females_aza.txt", sep = "\t")
female_aza <- female_aza[,c("study_id", "thrombocytes", "erythrocytes", "leucocytes")]
female_aza$analysis <- "aza"

male_pred <- read.table("Output_files/Lab_parameters/Males_pred.txt", sep = "\t")
male_pred <- male_pred[,c("study_id", "thrombocytes", "erythrocytes", "leucocytes")]
male_pred$analysis <- "pred"
female_pred <- read.table("Output_files/Lab_parameters/Females_pred.txt", sep = "\t")
female_pred <- male_pred[,c("study_id", "thrombocytes", "erythrocytes", "leucocytes")]
female_pred$analysis <- "pred"


#Female
female <- rbind(female_antiTNF, female_aza, female_pred)

female_plot_long <- female %>% 
  pivot_longer(cols = c( leucocytes, erythrocytes),
               names_to = "lab",
               values_to = "counts")


female_plot <- ggplot(female_plot_long, aes(analysis, counts)) + 
  geom_boxplot() + 
  xlab("") +
  ylab("Counts/10^-9 L") +
  labs(title = element_text("Lab parameters - female")) +
  stat_compare_means(label = "p.adj", comparisons = list(c("pred", "aza", "antiTNF")), method = "t.test", show.legend = T) +
  facet_grid("~lab") +
  theme_bw()



# Male
male <- rbind(male_antiTNF, male_aza, male_pred)

male_plot_long <- male %>% 
  pivot_longer(cols = c( leucocytes, erythrocytes),
               names_to = "lab",
               values_to = "counts")


male_plot <- ggplot(male_plot_long, aes(analysis, counts)) + 
  geom_boxplot() + 
  xlab("") +
  ylab("Counts/10^-9 L") +
  labs(title = element_text("Lab parameters - male")) +
  stat_compare_means(label = "p.adj", comparisons = list(c("pred", "aza", "antiTNF")), method = "wilcox.test", show.legend = T) +
  facet_grid("~lab") +
  theme_bw()

