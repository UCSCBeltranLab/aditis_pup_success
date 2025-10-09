##Load all necessary packages

library(tidyverse)
library(lubridate)
library(utils)
library(dplyr)
library(glmmTMB)
library(splines)
library(DHARMa)
library(ggeffects)
library(performance)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

################################ Creating proportion MPA and intrinsic variables ###############################
##setwd
setwd("/Users/aditimary/Documents/aditis_pup_success/RawData")

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

##update the metadata with MPA
metadata <- metadata %>%
  left_join(Proportion_MPA, by = c("animalID", "season"))

##assign each female to a harem in a given season
area_counts <- metadata %>%
  group_by(animalID, season, area) %>%
  summarize(count = n(), .groups = "drop") #count number of times each animal was seen in each area

##for each animal, find the area with the maximum count
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

##year born = first season seen - Age in that season
year_born <- age_in_first_season %>%
  group_by(animalID) %>%
  summarize(year_born = season - AgeYears)

##update the metadata with year_born
metadata <- metadata %>%
  left_join(year_born, by = "animalID")

##make table with intrinsic variables
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

##calculate arrival date
arrival <- metadata %>%
  group_by(animalID, season) %>%
  summarize(arrival_date = first(as.Date(date))) %>%
  mutate(arrival_md = format(arrival_date, "%m-%d"))

######################### Plotting general distributions #############################

##visualizing sample size is useful for determining stringency

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

##plot for BirthDate
ggplot(intrinsic_variables, aes(x = BirthDate)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(x = "Birth Date", y = "Number of Animals",
       title = "Distribution of Birth Dates") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##proportion MPA by birthdate
ggplot(data = intrinsic_variables, aes(x = BirthDate, y = proportion)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(x = "Birth Date", y = "Proportion",
       title = "MPA distribution by Birth Date") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##general distribution of proportions
ggplot(data = intrinsic_variables, aes(x = proportion)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "darkblue") +
  labs(title = "Histogram of Proportion MPA",
       x = "Proportion MPA", y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_few()

##graph age against proportion
ggplot(data = intrinsic_variables, aes(x = AgeYears, y = proportion)) +
  geom_point() +
  labs(title = "Proportion MPA by Age", x = "Age", y = "MPA") +
  theme_few()

##graph season against proportion
ggplot(data = intrinsic_variables, aes(x = season_fct, y = proportion)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(title = "MPA distribution by season", x = "Season", y = "Proportion")

##graph year born against proportion
ggplot(data = intrinsic_variables, aes(x = year_born_fct, y = proportion)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.5) + 
  labs(title = "MPA distribution by season", x = "Season", y = "Proportion")

######################## Modifying intrinsic variables needed for the model approach ############################

##first, modify intrinsic_variables so it has a column for proportion = 1 (1) OR not 1 (0)
intrinsic_variables <- intrinsic_variables %>%
  mutate(is_one = as.numeric(proportion == 1))

##create a version of intrinsic_variables with proportion < 1
intrinsic_variables_sub <- intrinsic_variables %>%
  filter(proportion < 1)

##transform proportion for this filtered data to make it suitable for gamma distribution 
flipped_prop <- max(intrinsic_variables_sub$proportion) - intrinsic_variables_sub$proportion #subtract all proportions from the maximum
flipped_prop <- flipped_prop - min(flipped_prop) + 0.0000001 #ensure no 0s by adding small constant
intrinsic_variables_sub$flipped_prop <- flipped_prop

##histogram of flipped proportion
ggplot(data = intrinsic_variables_sub, aes(x = flipped_prop)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "darkblue") +
  labs(title = "Histogram of Flipped Proportion MPA <1",
       x = "Flipped Proportion MPA", y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_few()

########################### Wave energy and tide data #############################

##first, see what the earliest and latest dates are in the breeding season
dates <- metadata %>%
  group_by(season) %>%
  summarize(earliest_date = min(date, na.rm = TRUE), latest_date = max(date, na.rm = TRUE))

##list only CSV files containing "tide" in the filename
tide_files <- list.files(pattern = "tide.*\\.csv$")

##read and combine only the tide data files
tide_data_all <- tide_files %>%
  map_dfr(read_csv)

##fit tide data to specific breeding seasons
tide_data_all <- tide_data_all %>%
  mutate(
    Date = as.Date(Date),
    MM = lubridate::month(Date),
    YY = lubridate::year(Date),
    season = case_when(
      MM == 12 ~ YY + 1, #December belongs to next year's season
      MM %in% c(1, 2, 3) ~ YY))  #Jan–Mar stay in current year

##select necessary columns from the tide data
tide_data_clean <- tide_data_all %>%
  select("Date", "Time (GMT)", "Verified (ft)", "season")

##list only CSV files containing "wave" in the filename
wave_files <- list.files(pattern = "wave.*\\.csv$")

##read and combine only the wave data files
wave_data_all <- wave_files %>%
  map_dfr(~ read.csv(.x, colClasses = "character"))

##convert all wave columns to numeric
wave_data_all <- wave_data_all %>%
  mutate(
    YY = as.numeric(trimws(YY)),
    MM = as.numeric(trimws(MM)),
    DD = as.numeric(trimws(DD)),
    hh = as.numeric(trimws(hh)),
    WVHT = as.numeric(trimws(WVHT)),
    DPD = as.numeric(trimws(DPD))) %>%
  filter(WVHT != 99, DPD != 99) #99 is used as an NA placeholder in the data, so we can filter these out

##fit wave data to specific breeding seasons
wave_data_all <- wave_data_all %>%
  mutate(season = case_when(
    MM == 12 ~ YY + 1,   #December data assigned to next year's breeding season
    MM %in% c(1,2,3) ~ YY,  #Jan-Mar data stay in current breeding season
    TRUE ~ NA_real_))      #Other months not in any breeding season, assign NA

##suppress scientific notation
options(scipen = 999)

##select necessary columns from the wave data
wave_data_clean <- wave_data_all %>%
  select("YY", "MM", "DD", "hh", "WVHT", "DPD", "season") %>%
  mutate(wave_power = 0.49 * WVHT^2 * DPD) #calculate wave power from the wave power equation

##need date and time columns to match in tide and wave data
#cleaned tide data
tide_data_clean <- tide_data_clean %>%
  mutate(tide_datetime = as.POSIXct(
    paste(Date, `Time (GMT)`), #Combine the `Date` and `Time (GMT)` columns into one string- "2023-10-02 14:00:00"
    format = "%Y-%m-%d %H:%M:%S", #Tell R the format of the string being converted
    tz = "UTC")) #Set the timezone to UTC

#cleaned wave data
wave_data_clean <- wave_data_clean %>%
  mutate(wave_datetime = as.POSIXct(
    sprintf("%04d-%02d-%02d %02d:00:00",  #Format string to create a datetime: "YYYY-MM-DD HH:00:00"
            YY, #Year (e.g., 2023), zero-padded to 4 digits
            MM, #Month (e.g., 1 becomes 01), zero-padded to 2 digits
            DD, #Day of month, zero-padded to 2 digits
            hh), #Hour (0–23), zero-padded to 2 digits
    format = "%Y-%m-%d %H:%M:%S", #Tell R the format of the string being converted
    tz = "UTC")) #Set the timezone to UTC

##join the data sets by season
tide_wave <- left_join(wave_data_clean, tide_data_clean, by = c("wave_datetime" = "tide_datetime", "season")) %>%
  filter(!is.na(`Verified (ft)`)) #make sure there are no NAs for tide height

##set extreme wave and tide threshold levels
extreme_wave_threshold <- 30 #30 Kw/h
extreme_tide_threshold <- 6 #6 ft  

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
  geom_line(linewidth = 1.2) +
  geom_point() +
  labs(title = "Extreme Tide and Storm Events (1996–2025)",
       x = "Year", y = "Number of Extreme Events") +
  theme_minimal()

############################# Harem density #################################

##seal density csv read
seal.density <- read.csv("seal.density.csv")

##avg area density calculation
area_density <- seal.density %>%
  rename(area = Beach) %>%
  mutate(date = ymd(date), season = year(date)) %>%
  group_by(area, season) %>%
  summarize(avg_density = mean(density), .groups = "drop")

##check what area density looks like per-season at each location
ggplot(area_density, aes(x = area, y = avg_density, fill = area)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ season, scales = "free_x") +
  labs(x = "Location",
       y = "Average Density",
       title = "Harem Density per Location by Season") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

##only need season after 2016 to match the harem density data
metadata_filtered <- metadata %>%
  filter(season >= 2016)

##assign each female to a harem in a given season
area_counts <- metadata_filtered %>%
  group_by(animalID, season, area) %>%
  summarize(count = n(), .groups = "drop") #count number of times each animal was seen in each area

##for each animal, find the area with the maximum count
harem_assignment <- area_counts %>%
  group_by(animalID, season) %>%
  filter(count == max(count)) %>% #these are assigned areas for each female
  ungroup()

##join extrinsic harem data with metadata to include harem assignment and area density
extrinsic_harem <- metadata_filtered %>%
  select(animalID, season, proportion, total_resights, AgeYears) %>%
  distinct() %>%
  filter(!is.na(proportion), proportion > 0) %>%
  left_join(harem_assignment, by = c("animalID", "season")) %>%
  left_join(area_density, by = c("area", "season")) %>%
  mutate(is_one = as.numeric(proportion == 1))

##check whether there are ties for assigned harems
harem_duplicates <- extrinsic_harem %>%
  group_by(season, animalID) %>%
  filter(n() > 1) %>%
  arrange(animalID, season) #there are duplicates

##break ties by the harem with higher density
extrinsic_harem <- extrinsic_harem %>%
  group_by(season, animalID) %>%
  slice_max(avg_density, with_ties = FALSE) %>%
  ungroup()

##check locations in the data that don't match with the ones we have density for from the Ano map
unique(extrinsic_harem$area)
unique(seal.density$Beach) #there are location mismatches

##remove all harem density NAs for now
extrinsic_harem <- extrinsic_harem %>%
  filter(!is.na(avg_density))

##create extrinsic variables with tide/wave data and harem density
extrinsic_variables <- extrinsic_harem %>%
  left_join(tide_wave_flagged, by = "season") %>%
  mutate(season_fct = factor(season)) %>%
  mutate(area_fct = factor(area),
         animalID_fct = factor(animalID))

########################### Modifying extrinsic variables needed for the model approach ################################

##first, modify intrinsic_variables so it has a column for proportion = 1 (1) OR not 1 (0)
extrinsic_variables <- extrinsic_variables %>%
  mutate(is_one = as.numeric(proportion == 1))

##create a version of extrinsic_variables with proportion < 1
extrinsic_variables_sub <- extrinsic_variables %>%
  filter(proportion < 1)

##transform proportion for this filtered data to make it suitable for gamma distribution 
flipped_prop <- max(extrinsic_variables_sub$proportion) - extrinsic_variables_sub$proportion #subtract all proportions from the maximum
flipped_prop <- flipped_prop - min(flipped_prop) + 0.0000001 #ensure no 0s by adding small constant
extrinsic_variables_sub$flipped_prop <- flipped_prop

##histogram of flipped proportion
ggplot(data = extrinsic_variables_sub, aes(x = flipped_prop)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "darkblue") +
  labs(title = "Histogram of Flipped Proportion MPA <1",
       x = "Flipped Proportion MPA", y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_few()

