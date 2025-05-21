library(tidyverse)
library(lubridate)
library(RMySQL)
library(utils)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

##read both csvs!
raw_data <- read.csv("Aditi 2024 Data Pull 2025_04_17_RAW.csv")
summarized_data <- read.csv("Aditi 2024 Data Pull 2025_04_17_SUMMARIZED.csv")

##raw data
raw_data <- raw_data %>%
  arrange(animalID) %>%
  group_by(animalID, date) %>%
  slice_min(order_by = entry, n = 1, with_ties = FALSE) %>%
  ungroup() #this gets rid of double resights (same animal observed twice on same day)

##create a metadata table combining summarized and raw data
metadata <- summarized_data %>%
  filter(!is.na(BirthDate)) %>% #only known birthdate moms
  filter(Precision <= 7) %>% #precision within a week
  filter(AgeYears <= 24) %>% #removes strange age animals
  left_join(raw_data, by = c("animalID", "season"))

##calculate number of days resighted for each animalID per season
total_resights <- metadata %>%
  group_by(season, animalID) %>%
  filter(as.Date(date) >= as.Date(BirthDate)) %>%
  summarize(total_resights = n()) %>%
  filter(total_resights >= 15)

##calculate number of times mom was seen with 1 pup
count_1_pup <- metadata %>%
  filter(withpup == 1) %>%
  group_by(animalID, season) %>%
  summarize(count_1_pup = n())

#calculate mother-pup association strength (1s out of total resights)
Proportion_MPA <- total_resights %>%
  left_join(count_1_pup, by = c("animalID", "season")) %>%
  mutate(count_1_pup = replace_na(count_1_pup, 0),
         proportion = count_1_pup / total_resights)

#update the metadata with MPA
metadata <- metadata %>%
  left_join(Proportion_MPA, by = c("animalID", "season"))

#violin plot for distribution of MPA across all seasons
ggplot(data = Proportion_MPA, aes(x = factor(season), y = proportion)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.5) + 
  labs(title = "MPA distribution by season", x = "Season", y = "Proportion") +
  theme_few()

#sample size per season
sample_size_per_season <- Proportion_MPA %>%
  group_by(season) %>%
  summarise(sample_size = n_distinct(animalID))

#create a plot for sample size per season
ggplot(data = sample_size_per_season, aes(x = factor(season), y = sample_size)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(title = "Sample size for each season", x = "Season", y = "Sample size") +
  theme_few()

#sample size per Age
sample_size_per_age <- metadata %>%
  group_by(AgeYears) %>%
  summarise(sample_size = n_distinct(animalID))

#create a plot for sample size by age per season
ggplot(data = sample_size_per_age, aes(x = AgeYears, y = sample_size)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(title = "Sample Size for each Age", x = "Age", y = "Sample Size") +
  theme_few()

#Assign each female to a harem in a given season
#First count number of times each animal was seen in each area
harem_assignment <- metadata %>%
group_by(animalID, season, area) %>%
  summarise(count = n(), .groups = "drop")

# For each animal, find the location with the maximum count
most_frequent_harem <- harem_assignment %>%
  group_by(animalID) %>%
  filter(count == max(count)) %>%
  ungroup()

#FOR FUN LM
#Calculate arrival date
arrival_date <- metadata %>%
  group_by(animalID) %>%
  summarize(ArrivalDate = first(as.Date(date))) 

ggplot(data = metadata, aes(x = AgeYears, y = proportion)) +
  geom_point() +
  labs(title = "Proportion MPA by Age", x = "Age", y = "MPA") +
  theme_few()

mod_Age_by_Proportion <- lm(AgeYears ~ proportion, data = metadata)
summary(mod_Age_by_Proportion)
