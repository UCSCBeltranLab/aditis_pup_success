library(tidyverse)
library(lubridate)
library(utils)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

####################### wave energy + tides #############################

#First, see what the earliest and latest dates are in the breeding season
dates <- metadata %>%
  group_by(season) %>%
  summarize(earliest_date = min(date, na.rm = TRUE), latest_date = max(date, na.rm = TRUE))

#List only CSV files containing "tide" in the filename
tide_files <- list.files(pattern = "tide.*\\.csv$")

#Read and combine only the tide data files
tide_data_all <- tide_files %>%
  map_dfr(read_csv)

#Fit tide data to specific breeding seasons
tide_data_all <- tide_data_all %>%
  mutate(
    Date = as.Date(Date),
    MM = lubridate::month(Date),
    YY = lubridate::year(Date),
    season = case_when(
      MM == 12 ~ YY + 1, #December belongs to next year's season
      MM %in% c(1, 2, 3) ~ YY))  #Jan–Mar stay in current year

#Select necessary columns from the tide data
tide_data_clean <- tide_data_all %>%
  select("Date", "Time (GMT)", "Verified (ft)", "season")

#List only CSV files containing "wave" in the filename
wave_files <- list.files(pattern = "wave.*\\.csv$")

#Read and combine only the wave data files
wave_data_all <- wave_files %>%
  map_dfr(~ read.csv(.x, colClasses = "character"))

#Convert all wave columns to numeric
wave_data_all <- wave_data_all %>%
  mutate(
    YY = as.numeric(trimws(YY)),
    MM = as.numeric(trimws(MM)),
    DD = as.numeric(trimws(DD)),
    hh = as.numeric(trimws(hh)),
    WVHT = as.numeric(trimws(WVHT)),
    DPD = as.numeric(trimws(DPD))) %>%
  filter(WVHT != 99, DPD != 99) #99 is used as an NA placeholder in the data, so we can filter these out

#Fit wave data to specific breeding seasons
wave_data_all <- wave_data_all %>%
  mutate(season = case_when(
    MM == 12 ~ YY + 1,   #December data assigned to next year's breeding season
    MM %in% c(1,2,3) ~ YY,  #Jan-Mar data stay in current breeding season
    TRUE ~ NA_real_))      #Other months not in any breeding season, assign NA

#suppress scientific notation
options(scipen = 999)

#Select necessary columns from the wave data
wave_data_clean <- wave_data_all %>%
  select("YY", "MM", "DD", "hh", "WVHT", "DPD", "season") %>%
  mutate(wave_power = 0.49 * WVHT^2 * DPD) #calculate wave power from the wave power equation

##Need date and time columns to match in tide and wave data
tide_data_clean <- tide_data_clean %>%
  mutate(tide_datetime = as.POSIXct(
  paste(Date, `Time (GMT)`), #Combine the `Date` and `Time (GMT)` columns into one string- "2023-10-02 14:00:00"
  format = "%Y-%m-%d %H:%M:%S", #Tell R the format of the string being converted
  tz = "UTC")) #Set the timezone to UTC

wave_data_clean <- wave_data_clean %>%
  mutate(wave_datetime = as.POSIXct(
    sprintf("%04d-%02d-%02d %02d:00:00",  #Format string to create a datetime: "YYYY-MM-DD HH:00:00"
            YY, #Year (e.g., 2023), zero-padded to 4 digits
            MM, #Month (e.g., 1 becomes 01), zero-padded to 2 digits
            DD, #Day of month, zero-padded to 2 digits
            hh), #Hour (0–23), zero-padded to 2 digits
    format = "%Y-%m-%d %H:%M:%S", #Tell R the format of the string being converted
    tz = "UTC")) #Set the timezone to UTC

#Join the data sets by season
tide_wave <- left_join(wave_data_clean, tide_data_clean, by = c("wave_datetime" = "tide_datetime", "season")) %>%
  filter(!is.na(`Verified (ft)`)) #make sure there are no NAs for tide height

#Set extreme wave and tide threshold levels
extreme_wave_threshold <- 30 #30 Kw/h
extreme_tide_threshold <- 6 #6 ft  

#Flag cases where both wave power and tide were extreme
tide_wave_flagged <- tide_wave %>%
  mutate(extreme_wave = wave_power >= extreme_wave_threshold,
    extreme_tide = `Verified (ft)` >= extreme_tide_threshold,
    extreme_both = extreme_wave & extreme_tide) %>%
  group_by(season) %>%
  summarize(n_extreme_both = sum(extreme_both),
            n_extreme_tide = sum(extreme_tide),
            n_extreme_wave = sum(extreme_wave)) #summarize occurrences of both per season

#Plot extreme tide and wave events per season
ggplot(tide_wave_flagged, aes(x = season, y = n_extreme_both)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  labs(title = "Extreme Tide and Storm Events (1996–2025)",
       x = "Year", y = "Number of Extreme Events") +
  theme_minimal()

############################# harem density #################################

#Seal density csv read
seal.density <- read.csv("seal.density.csv")

#Avg area density calculation
area_density <- seal.density %>%
  rename(area = Beach) %>%
  mutate(date = ymd(date), season = year(date)) %>%
  group_by(area, season) %>%
  summarize(avg_density = mean(density), .groups = "drop")

#check what area density looks like per-season at each location
ggplot(area_density, aes(x = area, y = avg_density, fill = area)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ season, scales = "free_x") +
  labs(x = "Location",
       y = "Average Density",
       title = "Harem Density per Location by Season") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

#Only need season after 2016 to match the harem density data
metadata_filtered <- metadata %>%
  filter(season >= 2016)

#Assign each female to a harem in a given season
area_counts <- metadata_filtered %>%
  group_by(animalID, season, area) %>%
  summarize(count = n(), .groups = "drop") #count number of times each animal was seen in each area

#For each animal, find the area with the maximum count
harem_assignment <- area_counts %>%
  group_by(animalID, season) %>%
  filter(count == max(count)) %>% #these are assigned areas for each female
  ungroup()

#Join extrinsic harem data with metadata to include harem assignment and area density
extrinsic_harem <- metadata_filtered %>%
  select(animalID, season, proportion, total_resights, AgeYears) %>%
  distinct() %>%
  filter(!is.na(proportion), proportion > 0) %>%
  left_join(harem_assignment, by = c("animalID", "season")) %>%
  left_join(area_density, by = c("area", "season")) %>%
  mutate(is_one = as.numeric(proportion == 1))

#Check whether there are ties for assigned harems
harem_duplicates <- extrinsic_harem %>%
  group_by(season, animalID) %>%
  filter(n() > 1) %>%
  arrange(animalID, season) #there are duplicates

#Break ties by the harem with higher density
extrinsic_harem <- extrinsic_harem %>%
  group_by(season, animalID) %>%
  slice_max(avg_density, with_ties = FALSE) %>%
  ungroup()

#Check locations in the data that don't match with the ones we have density for from the Ano map
unique(extrinsic_harem$area)
unique(seal.density$Beach) #there are location mismatches

#Remove NAs for now (come back to this later)
extrinsic_harem <- extrinsic_harem %>%
  filter(!is.na(avg_density))

#create extrinsic variables with tide/wave data and harem density
extrinsic_variables <- extrinsic_harem %>%
  left_join(tide_wave_flagged, by = "season") %>%
  mutate(season_fct = factor(season)) %>%
  mutate(area_fct = factor(area))

#create a version of extrinsic_variables with proportion < 1
extrinsic_variables_sub <- extrinsic_variables %>%
  filter(proportion < 1)

##transform proportion for this filtered data to make it suitable for gamma distribution 
flipped_prop <- max(extrinsic_variables_sub$proportion) - extrinsic_variables_sub$proportion #subtract all proportions from the maximum
flipped_prop <- flipped_prop - min(flipped_prop) + 0.0000001 #ensure no 0s by adding small constant
extrinsic_variables_sub$flipped_prop <- flipped_prop









##Data in case we want a separate weather model
#Join extrinsic weather data with metadata for separate modeling
extrinsic_weather <- metadata %>%
  select(animalID, season, proportion, total_resights, AgeYears) %>%
  distinct() %>%
  filter(!is.na(proportion), proportion > 0) %>%
  left_join(tide_wave_flagged, by = "season") %>%
  mutate(is_one = as.numeric(proportion == 1),
         season_fct = factor(season))

#create a version of extrinsic_variables with proportion < 1
extrinsic_weather_sub <- extrinsic_weather %>%
  filter(proportion < 1)

##transform proportion for this filtered data to make it suitable for gamma distribution 
flipped_prop <- max(extrinsic_weather_sub$proportion) - extrinsic_weather_sub$proportion #subtract all proportions from the maximum
flipped_prop <- flipped_prop - min(flipped_prop) + 0.0000001 #ensure no 0s by adding small constant
extrinsic_weather_sub$flipped_prop <- flipped_prop





