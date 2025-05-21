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
  filter(as.Date(date) >= as.Date(BirthDate)) %>% #only count resights after Birth Date
  summarize(total_resights = n()) %>%
  filter(total_resights >= 12) #stringency for total resights post-birth

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

##violin plot for distribution of MPA across all seasons
ggplot(data = Proportion_MPA, aes(x = factor(season), y = proportion)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.5) + 
  labs(title = "MPA distribution by season", x = "Season", y = "Proportion") +
  theme_few()

##assign each female to a harem in a given season
area_counts <- metadata %>%
  group_by(animalID, season, area) %>%
  summarise(count = n(), .groups = "drop") #count number of times each animal was seen in each area

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
  summarise(year_born = season - AgeYears)

##Update the metadata with year_born
metadata <- metadata %>%
  left_join(year_born, by = "animalID")

##Make table with intrinsic factors
intrinsic_variables <- metadata %>%
  select(animalID, season, AgeYears, BirthDate, year_born, proportion) %>%
  distinct()

##Make a table with extrinsic factors
#for now, we only have harem assignment
#NEED: weather/tidal, harem density for each harem
extrinsic_variables <- metadata

#Visualizing trends in arrival and pupping dates
#ask Maddie how to visualize this!
#Calculate arrival date
arrival_date <- metadata %>%
  group_by(animalID, season) %>%
  summarize(arrival_date = first(as.Date(date)))

#trend for arrival date across years
ggplot(data = arrival_date, aes(x = factor(season), y = arrival_date)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.5) + 
  labs(title = "Arrival Date by Season", x = "Season", y = "Arrival Date") +
  theme_few()
  
#For fun; modeling different things
##Graph age against proportion
ggplot(data = intrinsic_variables, aes(x = AgeYears, y = proportion)) +
  geom_point() +
  labs(title = "Proportion MPA by Age", x = "Age", y = "MPA") +
  theme_few()

#Linear model- this won't work as proportion is skewed, not normally distributed
mod_Age_by_Proportion <- lm(proportion ~ AgeYears, data = metadata)
summary(mod_Age_by_Proportion)

#The qqplot to prove it in case u don't believe me
qqnorm(metadata$proportion, pch = 1, frame = FALSE)
qqline(metadata$proportion, col = "steelblue", lwd = 2)

##Poisson model for age and proportion - a little better
mod_Poisson_Age_by_Proportion <- glm(proportion ~ AgeYears, data = metadata, family = "poisson")
summary(mod_Poisson_Age_by_Proportion)

##sample size per season
sample_size_per_season <- Proportion_MPA %>%
  group_by(season) %>%
  summarise(sample_size = n_distinct(animalID))

##create a plot for sample size per season
ggplot(data = sample_size_per_season, aes(x = factor(season), y = sample_size)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(title = "Sample size for each season", x = "Season", y = "Sample size") +
  theme_few() + 
  scale_y_continuous(n.breaks = 10)

##sample size per Age
sample_size_per_age <- metadata %>%
  group_by(AgeYears) %>%
  summarise(sample_size = n_distinct(animalID))

##create a plot for sample size by age
ggplot(data = sample_size_per_age, aes(x = AgeYears, y = sample_size)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(title = "Sample Size for each Age", x = "Age", y = "Sample Size") +
  theme_few()






