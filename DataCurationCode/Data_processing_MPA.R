library(tidyverse)
library(lubridate)
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
  filter(as.Date(date) >= as.Date(BirthDate)) %>% #only count resights after Birth Date
  summarize(total_resights = n()) %>%
  filter(total_resights >= 14) #stringency for total resights post-birth

##calculate number of times mom was seen with 1 pup
count_1_pup <- metadata %>% 
  filter(withpup == 1) %>%
  group_by(animalID, season) %>%
  summarize(count_1_pup = n())

##calculate mother-pup association strength (1s out of total resights)
Proportion_MPA <- total_resights %>%
  left_join(count_1_pup, by = c("animalID", "season")) %>%
  mutate(count_1_pup = replace_na(count_1_pup, 0),
         proportion = count_1_pup / total_resights)

##Update the metadata with MPA
metadata <- metadata %>%
  left_join(Proportion_MPA, by = c("animalID", "season"))

##assign each female to a harem in a given season
area_counts <- metadata %>%
  group_by(animalID, season, area) %>%
  summarize(count = n(), .groups = "drop") #count number of times each animal was seen in each area

#for each animal, find the area with the maximum count
harem_assignment <- area_counts %>%
  group_by(animalID, season) %>%
  filter(count == max(count)) %>% #these are assigned areas for each female
  ungroup()

##calculate year born for each female
age_in_first_season <- metadata %>%
  group_by(animalID) %>%
  filter(season == min(season)) %>%
  slice(1) %>%
  ungroup() %>%
  select(animalID, AgeYears, season) #AgeYears in the first season each animal was seen

#Year born = first season seen - Age in that season
year_born <- age_in_first_season %>%
  group_by(animalID) %>%
  summarize(year_born = season - AgeYears)

##Update the metadata with year_born
metadata <- metadata %>%
  left_join(year_born, by = "animalID")

##Total pups in lifetime
total_pups_lifetime <- metadata %>%
  select(animalID, season) %>%
  distinct() %>%
  group_by(animalID) %>%
  summarize(total_pups = n(), .groups = "drop")

intrinsic_variables <- intrinsic_variables %>%
  left_join(total_pups_lifetime, by = "animalID")

##Make table with intrinsic factors
intrinsic_variables <- metadata %>%
  select(animalID, season, AgeYears, BirthDate, year_born, proportion, total_resights) %>%
  distinct() %>% #eliminate the metadata multiple rows; we only need one per animalID per season
  filter(!is.na(proportion), proportion > 0) %>% #eliminate NA and 0s for proportion
  mutate(BirthDate = as.Date(BirthDate)) %>%
  group_by(animalID) %>%  
  mutate(year_born_fct = factor(year_born), #make year_born a factor in separate column
         animalID_fct = factor(animalID), #make animalID a factor in separate column
         season_fct = factor(season),
         age_last_seen = max(AgeYears), #calculate age at last observation
         BirthDate = format(BirthDate, "%m-%d")) %>% #format BirthDate as month-day
  ungroup()

##Calculate arrival date
arrival <- metadata %>%
  group_by(animalID, season) %>%
  summarize(arrival_date = first(as.Date(date))) %>%
  mutate(arrival_md = format(arrival_date, "%m-%d"))

##Visualizing sample size is useful for determining stringency
##sample size per season
sample_size_per_season <- intrinsic_variables %>%
  group_by(season) %>%
  summarize(sample_size = n_distinct(animalID))

##sample size per age
sample_size_per_age <- intrinsic_variables %>%
  group_by(AgeYears) %>%
  summarize(sample_size = n_distinct(animalID))

##create a plot for sample size per season
ggplot(data = sample_size_per_season, aes(x = season, y = sample_size)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(title = "Sample size for each season", x = "Season", y = "Sample size") +
  theme_few() + 
  scale_y_continuous(n.breaks = 10)

##sample size per age
sample_size_per_age <- intrinsic_variables %>%
  group_by(AgeYears) %>%
  summarize(sample_size = n_distinct(animalID))

##create a plot for sample size by age
ggplot(data = sample_size_per_age, aes(x = AgeYears, y = sample_size)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(title = "Sample Size for each Age", x = "Age", y = "Sample Size") +
  theme_few() +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 20)

ggplot(intrinsic_variables, aes(x = BirthDate)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(x = "Birth Date", y = "Number of Animals",
       title = "Distribution of Birth Dates") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Proportion MPA by birthdate
ggplot(data = intrinsic_variables, aes(x = BirthDate, y = proportion)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(x = "Birth Date", y = "Proportion",
       title = "MPA distribution by Birth Date") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##General distribution of proportions
ggplot(data = intrinsic_variables, aes(x = proportion)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "darkblue") +
  labs(title = "Histogram of Proportion MPA",
       x = "Proportion MPA", y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_few()

##Graph age against proportion
ggplot(data = intrinsic_variables, aes(x = AgeYears, y = proportion)) +
  geom_point() +
  labs(title = "Proportion MPA by Age", x = "Age", y = "MPA") +
  theme_few()

##Graph season against proportion
ggplot(data = intrinsic_variables, aes(x = season_fct, y = proportion)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(title = "MPA distribution by season", x = "Season", y = "Proportion")

##Graph year born against proportion
ggplot(data = intrinsic_variables, aes(x = year_born_fct, y = proportion)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.5) + 
  labs(title = "MPA distribution by season", x = "Season", y = "Proportion")

