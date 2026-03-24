##Load all necessary packages

library(tidyverse)
library(lubridate)
library(utils)
library(dplyr)
library(glmmTMB)
library(lme4)
library(DHARMa)
library(ggeffects)
library(performance)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(parallel)
library(binom)
library(gt)
library(flextable)
library(grid)
library(cowplot)
library(viridisLite)

################################ Defining MPA + intrinsic data processing ###############################

##read both csvs!
raw_data <- read.csv("./RawData/Aditi 2024 Data Pull 2025_04_17_RAW.csv") 
# the raw resight data: includes animalID, season, date, area, entry date, entry (person), withpup
summarized_data <- read.csv("./RawData/Aditi 2024 Data Pull 2025_04_17_SUMMARIZED.csv")
# specific per-animal data: includes animalID, age, season, birth date with precision

summary_counts <- raw_data %>%
  mutate(date = as.Date(date)) %>%
  summarise(n_seals = n_distinct(animalID), n_sightings = n())

##ensure that each animalID is only seen once per day (accounting for double resights)
raw_data <- raw_data %>%
  arrange(animalID) %>%
  group_by(animalID, date) %>%
  slice_min(order_by = entry, n = 1, with_ties = FALSE) %>% #gets rid of double resights
  filter(withpup %in% c(0, 1)) %>% #restrict to only values of wtihpup score = 0 or 1
  ungroup() 

##calculate prior pupping experience from raw data
pupping_exp <- raw_data %>%
  group_by(animalID, season) %>%
  summarize(had_pup = as.integer(any(withpup == 1)), #any at least 1 observation of 1 pup counts as pupping experience
            .groups = "drop") %>%
  group_by(animalID) %>%
  arrange(season, .by_group = TRUE) %>%  #so that season is chronological
  mutate(experience_prior = lag(cumsum(had_pup == 1), default = 0)) %>% #calculate prior pup experience
  ungroup()

#join with raw data to include prior experience
raw_data <- raw_data %>%
  left_join(pupping_exp, by = c("animalID", "season"))

##create a metadata table combining summarized and raw data
metadata <- summarized_data %>%
  filter(!is.na(BirthDate)) %>%       # only known birthdate moms
  filter(Precision <= 7) %>%          # precision within a week
  filter(AgeYears <= 24) %>%          # removes strange age animals
  left_join(raw_data, by = c("animalID", "season")) %>%   # <-- missing %>% before
  group_by(AgeYears) %>%
  ungroup()

metadata <- metadata %>%
  group_by(AgeYears) %>%
  mutate(mean_exp_age = mean(experience_prior, na.rm = TRUE)) %>%
  ungroup()

metadata <- metadata %>%
  mutate(residual_experience = experience_prior - mean_exp_age)

### Resight effort and calculating MPA ###

# 1) calculate number of days resighted for each animalID per season (only post-bith date to standardize across birth timing)
total_resights <- metadata %>%
  group_by(season, animalID) %>%
  filter(as.Date(date) >= as.Date(BirthDate)) %>% #only count resights after Birth Date
  summarize(total_resights = n()) %>% #Count total_resights after filtering
  filter(total_resights >= 14) #stringency for total resights post-birth; at least 14 days is a good middle ground 

# 2) calculate number of times mom was seen with 1 pup
count_1_pup <- metadata %>% 
  filter(withpup == 1) %>%
  group_by(animalID, season) %>%
  summarize(count_1_pup = n())

# 3) calculate mom-pup association strength:
#Compute number of with-pup scores of 1 out of entire total resights calculated above 
#Result: proportion of pup association strength for each female each breeding season (variable = proportion)
Proportion_MPA <- total_resights %>%
  left_join(count_1_pup, by = c("animalID", "season")) %>%
  mutate(proportion = count_1_pup / total_resights)

##update the metadata to include MPA
metadata <- metadata %>%
  left_join(Proportion_MPA, by = c("animalID", "season"))

### Assigning each female harem locations ###

# 1) Count sightings for each female in each area 
area_counts <- metadata %>%
  group_by(animalID, season, area) %>%
  summarize(count = n(), .groups = "drop") #count number of times each animal was seen in each area

## 2) for each animal, find the area with the maximum count
harem_assignment <- area_counts %>%
  group_by(animalID, season) %>%
  filter(count == max(count)) %>% #these are assigned areas for each female
  ungroup() %>%
  select(animalID, season, dominant_area = area) #call dominant area to avoid confusion

### Calculate cohort for each female ###

# 1) Calculate year born for each female
age_in_first_season <- metadata %>%
  # Identify each female's first observed season
  group_by(animalID) %>%
  filter(season == min(season)) %>%
  slice(1) %>%
  ungroup() %>%
  select(animalID, AgeYears, season) %>%
  mutate(
    # Ensure season and age are numeric to calculate year born in the next step
    season  = as.numeric(season),
    AgeYears = as.numeric(AgeYears))

# 2) Year born = first season seen - Age in that season
year_born <- age_in_first_season %>%
  mutate(year_born = season - AgeYears) %>%
  select(animalID, year_born) 

## update metadata with year_born
metadata <- metadata %>%
  left_join(year_born, by = "animalID")

metadata <- metadata %>%
  filter(!is.na(proportion))

########################### Wave energy and tide data processing #############################
##Identify seasonal date range: 
##Calculate earliest and latest dates are in the breeding season
##Used for: selecting dates when uploading the wave and tide data
dates <- metadata %>% 
  group_by(season) %>%
  summarize(earliest_date = min(date, na.rm = TRUE), latest_date = max(date, na.rm = TRUE))

# 1) Read in tide data 
#list only CSV files containing "tide" in the filename
tide_files <- list.files(path = "./RawData/", pattern = "tide.*\\.csv$", full.names = TRUE)

#read and combine all CSVs
tide_data_all <- tide_files %>%
  lapply(read_csv) %>%   # or read_csv2() depending on delimiter
  bind_rows()

#clean and create date variables
tide_data_all <- tide_data_all %>%
  mutate(Date = as.Date(as.character(Date)),
         MM   = month(Date),
         YY   = year(Date),
         season = case_when(MM == 12 ~ YY + 1,
                            MM %in% c(1, 2, 3) ~ YY,
                            TRUE ~ NA_real_))
# December observations belong to the *following* breeding season;
# Jan–Mar belong to the current year’s season

#select necessary columns from the tide data 
#Date: 
#Time: 
#Verified: 
#season: 
tide_data_clean <- tide_data_all %>%
  select("Date", "Time (GMT)", "Verified (ft)", "season")

# 2) Read in wave data 
#list only CSV files containing "wave" in the filename
wave_files <- list.files(path = "./RawData/", pattern = "wave.*\\.csv$", full.names = TRUE)

#read and combine only the wave data files
wave_data_all <- wave_files %>%
  map_dfr(~ read.csv(.x, colClasses = "character")) #everything as character makes it easier to join for now

#convert all wave columns to numeric
wave_data_all <- wave_data_all %>%
  mutate(YY = as.numeric(trimws(YY)), #trimws to make sure there are no weird blank spaces and make each column numeric
         MM = as.numeric(trimws(MM)),
         DD = as.numeric(trimws(DD)),
         hh = as.numeric(trimws(hh)),
         WVHT = as.numeric(trimws(WVHT)),
         DPD = as.numeric(trimws(DPD))) %>%
  filter(WVHT != 99, DPD != 99) #99 is used as an NA placeholder in the data, so we can filter these out

#fit wave data to specific breeding seasons
wave_data_all <- wave_data_all %>%
  mutate(season = case_when(MM == 12 ~ YY + 1, #December data assigned to next year's breeding season
                            MM %in% c(1,2,3) ~ YY, #Jan-Mar data stay in current breeding season
                            TRUE ~ NA_real_)) #Other months not in any breeding season, assign NA

##suppress scientific notation
options(scipen = 999)

#select necessary columns from the wave data
#YY: 
#MM: 
#DD: 
#hh:
#WVHT: 
#DPD:
#season: 
wave_data_clean <- wave_data_all %>%
  select("YY", "MM", "DD", "hh", "WVHT", "DPD", "season") %>%
  mutate(wave_power = 0.49 * WVHT^2 * DPD) #calculate wave power from the wave power equation: (kW/m) = 0.49 * H_s^2 * T_p

# 3) Build aligned datetime columns and join the tide and wave data 
#Tide: Combine data and GMT tide into POSIXct 
tide_data_clean <- tide_data_clean %>%
  mutate(tide_datetime = as.POSIXct(paste(Date, `Time (GMT)`), #combine date + time
                                    format = "%Y-%m-%d %H:%M", #match format
                                    tz = "UTC")) #set time zone to UTC

#Wave: Combine data and GMT tide into POSIXct 
wave_data_clean <- wave_data_clean %>%
  mutate(wave_datetime = as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00", #Format string to create a datetime: "YYYY-MM-DD HH:00:00"
                                            YY, #Year
                                            MM, #Month
                                            DD, #Day of month
                                            hh), #Hour (0–23)
                                    format = "%Y-%m-%d %H:%M:%S", #match format
                                    tz = "UTC")) #Set the timezone to UTC

##join the data sets by season
tide_wave <- left_join(wave_data_clean, tide_data_clean, by = c("wave_datetime" = "tide_datetime", "season")) %>%
  filter(!is.na(`Verified (ft)`)) #make sure there are no NAs for tide height

quantile(tide_wave$wave_power,
         probs = c(0.9, 0.95, 0.975, 0.99),
         na.rm = TRUE)

quantile(tide_wave$`Verified (ft)`,
         probs = c(0.9, 0.95, 0.975, 0.99),
         na.rm = TRUE)

# 4) Flag and categorize extreme events
##set extreme wave and tide threshold levels
extreme_wave_threshold <- 110.8388 #based on 0.95
extreme_tide_threshold <- 5.53 #based on 0.95

##flag cases where both wave power and tide were extreme
tide_wave_flagged <- tide_wave %>%
  mutate(extreme_wave = wave_power >= extreme_wave_threshold,
         extreme_tide = `Verified (ft)` >= extreme_tide_threshold,
         extreme_both = extreme_wave & extreme_tide) %>%
  group_by(season) %>%
  summarize(n_extreme_both = sum(extreme_both),
            n_extreme_tide = sum(extreme_tide),
            n_extreme_wave = sum(extreme_wave)) #summarize occurrences of both per season

##plot extreme tide and wave events per season
ggplot(tide_wave_flagged, aes(x = season, y = n_extreme_both)) +
  geom_line(linewidth = 1.2, color = "#1f78b4") +
  geom_point(color = "darkblue") +
  scale_x_continuous(n.breaks = 28) +
  labs(title = "Extreme Tide and Storm Events (1996–2025)",
       x = "Year", y = "Number of Extreme Events") +
  theme_minimal()

############################# Harem density data processing #################################
##run the source code from harem density processing
source("./DataCurationCode/Updated_UAS_harem_density.R")

# 1) seal density csv read
seal.density <- read.csv("./RawData/seal.density.csv")

# 2a) Calculate average density per area-season for modeling, then link to assigned harems
area_density <- seal.density %>%
  rename(dominant_area = Beach) %>% 
  mutate(date = ymd(date), season = year(date)) %>% 
  #transform the date column so that it's in date form with year = season to match the resight data
  group_by(dominant_area, season) %>%
  summarize(avg_density = mean(density)) %>% #calculate mean density per area in each season 
  ungroup() |>
  left_join(harem_assignment, by = c("season", "dominant_area")) %>% 
  filter(!is.na(animalID)) #only keep animals observed in the 2016-2025 dataset to match the drone data

# 2b) Data for plotting (one row per beach x season)
area_density_plot <- area_density %>%
  group_by(dominant_area, season) %>%
  summarize(
    avg_density = first(avg_density),
    n_animals   = n_distinct(animalID), 
    .groups = "drop")

##Quick check: area density per-season at each location
ggplot(area_density_plot, aes(x = dominant_area,
                         y = avg_density,
                         fill = avg_density)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ season, scales = "free_x") +
  scale_fill_gradient(low  = "#DEEBF7",
                      high = "navy") +
  labs(x = "Beach",
       y = "Average Density",
       fill = "Average Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3) Clean up the area density data
# Identify animal-season cases with multiple harems (ties)
harem_duplicates <- area_density %>%
  group_by(season, animalID) %>%
  filter(n() > 1) %>%
  arrange(animalID, season)

# Break ties by selecting the harem with highest average density
area_density <- area_density %>%
  group_by(season, animalID) %>%
  slice_max(avg_density, with_ties = FALSE) %>%
  ungroup()

# Check for remaining area mismatches (for manual cleaning if needed)
table(area_density$dominant_area)
table(seal.density$Beach)

##remove all harem density NAs because there are very few (<10) instances of areas that don't match the map locations
area_density <- area_density %>%
  filter(!is.na(avg_density))

###################### Build final modeling table ###################

##make table with intrinsic variables
intrinsic_variables <- metadata %>%
  ## left_join(harem_assignment, by = c("animalID", "season")) %>% 
  select(animalID, season, AgeYears, BirthDate, year_born, residual_experience, experience_prior, proportion, total_resights, count_1_pup) %>%
  distinct(animalID, season, .keep_all = TRUE) %>%  # get rid of duplicates so only 1 row per animal-season
  filter(!is.na(proportion), proportion > 0) %>%
  mutate(BirthDate = as.Date(BirthDate),
         birth_doy = lubridate::yday(BirthDate)) %>%
        ## group = case_when(grepl("^N", dominant_area) | dominant_area %in% c("BBNN","BBNS","BBN") ~ "NP",
                           # TRUE ~ "SP")) %>% #create a grouping factor for NP/SP location differences
  group_by(animalID) %>%
  mutate(animalID_fct = factor(animalID), #animalID factor
         season_fct = factor(season)) %>% #season factor
       ## optional group_fct = factor(group)) %>% #group factor
  ungroup() %>%
  left_join(area_density %>% select(animalID, season, dominant_area, avg_density), by = c("animalID", "season")) %>%
  left_join(tide_wave_flagged, by = "season")

######################## Modifying age variables needed for the model approach ############################

##Setting senescence threshold
age_senesce <- 9
  
##Setting threshold in data
intrinsic_variables <- intrinsic_variables %>%
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old"))) %>%
  mutate(age10 = (AgeYears - age_senesce) / 10) #scaled numeric version of age centered at senescence threshold

