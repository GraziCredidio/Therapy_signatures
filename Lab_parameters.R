# EZE cohort: Therapy Signatures
  # Female and Male erythrocytes and leukocytes counts - therapy patients (pred, aza, anti-TNF)
  # Supplementary Figure 2 (Master's thesis)
  # Author: Graziella Credidio

rm(list = ls())

# Loading packages ----
library(ggpubr)

# Loading files ----
coldata_antiTNF <- read.table("Cleaned_tables/models/antiTNF/EZECohort_coldata_antiTNF.txt", sep = "\t")
coldata_pred <- read.table("Cleaned_tables/models/pred/EZECohort_coldata_pred.txt", sep = "\t")
coldata_aza <- read.table("Cleaned_tables/models/aza/EZECohort_coldata_aza.txt", sep = "\t")

# Filtering data ----
# Males using therapies
male_antiTNF <- coldata_antiTNF %>% 
  filter(antiTNF_vs_noBiologics == 'anti_tnf') %>% 
  filter(sex == 'Male') %>% 
  select(c("study_id", "erythrocytes", "leucocytes"))
male_antiTNF$analysis <- "antiTNF"

male_pred <- coldata_pred %>% 
  filter(pred_vs_noSyst == 'pred') %>% 
  filter(sex == 'Male') %>% 
  select(c("study_id", "erythrocytes", "leucocytes"))
male_pred$analysis <- "pred"

male_aza <- coldata_aza %>% 
  filter(aza_vs_noSyst == 'aza') %>% 
  filter(sex == 'Male') %>% 
  select(c("study_id", "erythrocytes", "leucocytes"))
male_aza$analysis <- "aza"

# Females using therapies
female_antiTNF <- coldata_antiTNF %>% 
  filter(antiTNF_vs_noBiologics == 'anti_tnf') %>% 
  filter(sex == 'Female') %>% 
  select(c("study_id", "erythrocytes", "leucocytes"))
female_antiTNF$analysis <- "antiTNF"

female_pred <- coldata_pred %>% 
  filter(pred_vs_noSyst == 'pred') %>% 
  filter(sex == 'Female') %>% 
  select(c("study_id", "erythrocytes", "leucocytes"))
female_pred$analysis <- "pred"

female_aza <- coldata_aza %>% 
  filter(aza_vs_noSyst == 'aza') %>% 
  filter(sex == 'Female') %>% 
  select(c("study_id", "erythrocytes", "leucocytes"))
female_aza$analysis <- "aza"

# Creating dataframes to be plotted ----
male <- rbind(male_antiTNF, male_aza, male_pred)
female <- rbind(female_antiTNF, female_aza, female_pred)

female_plot_long <- female %>% 
  pivot_longer(cols = c(leucocytes, erythrocytes),
               names_to = "lab",
               values_to = "counts")

male_plot_long <- male %>% 
  pivot_longer(cols = c(leucocytes, erythrocytes),
               names_to = "lab",
               values_to = "counts")
# NAs in columns
sum(is.na(female_plot_long$counts)) #98 NAs
sum(is.na(male_plot_long$counts)) #34 NAs

# Plots ----
female_plot <- ggplot(female_plot_long, aes(analysis, counts)) + 
  geom_boxplot() + 
  xlab("") +
  ylab("Counts/10^-9 L") +
  labs(title = element_text("Lab parameters - Females")) +
  stat_compare_means(label = "p.adj", comparisons = list(c("pred", "aza", "antiTNF")), method = "t.test", show.legend = T) +
  facet_grid("~lab") +
  theme_bw()
female_plot
ggsave(female_plot, file = "Output_files/Lab_params/females_erythro_leuko.png", height = 8, width = 10, units = "in", dpi = 300)

male_plot <- ggplot(male_plot_long, aes(analysis, counts)) + 
  geom_boxplot() + 
  xlab("") +
  ylab("Counts/10^-9 L") +
  labs(title = element_text("Lab parameters - Males")) +
  stat_compare_means(label = "p.adj", comparisons = list(c("pred", "aza", "antiTNF")), method = "t.test", show.legend = T) +
  facet_grid("~lab") +
  theme_bw()
male_plot
ggsave(male_plot, file = "Output_files/Lab_params/males_erythro_leuko.png", height = 8, width = 10, units = "in", dpi = 300)
