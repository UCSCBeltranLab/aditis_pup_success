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
library(rptR)

################################ Calculating MPA and intrinsic data processing ###############################

##take away the set working directory

##read both csvs!
raw_data <- read.csv("./RawData/Aditi 2024 Data Pull 2025_04_17_RAW.csv")
summarized_data <- read.csv("./RawData/Aditi 2024 Data Pull 2025_04_17_SUMMARIZED.csv")

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

##filter for only values of withpup == 1 and 0
metadata <- metadata %>%
  filter(withpup %in% c(0, 1))

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

##calculate mom-pup association strength (1s out of total resights)
Proportion_MPA <- total_resights %>%
  left_join(count_1_pup, by = c("animalID", "season")) %>%
  mutate(proportion = count_1_pup / total_resights)

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
  select(animalID, AgeYears, season)

##ensure numeric columns before computing
age_in_first_season <- age_in_first_season %>%
  mutate(season = as.numeric(season),
    AgeYears = as.numeric(AgeYears))

## year born = first season seen - Age in that season
year_born <- age_in_first_season %>%
  mutate(year_born = season - AgeYears) %>%
  select(animalID, year_born)

## update metadata with year_born
metadata <- metadata %>%
  left_join(year_born, by = "animalID")

##make table with intrinsic variables
intrinsic_variables <- metadata %>%
  select(animalID, season, AgeYears, BirthDate, year_born, proportion, total_resights, count_1_pup) %>%
  distinct() %>% #eliminate the metadata multiple rows; we only need one per animalID per season
  filter(!is.na(proportion), proportion > 0) %>% #eliminate NA and 0s for proportion
  mutate(BirthDate = as.Date(BirthDate)) %>%
  group_by(animalID) %>%  
  mutate(year_born_num = as.numeric(year_born), #numeric year_born allows it to be included in the model
         animalID_fct = factor(animalID), #make animalID a factor in separate column
         season_fct = factor(season),
         age_last_seen = max(AgeYears), #calculate age at last observation
         BirthDate = format(BirthDate, "%m-%d")) %>% #format BirthDate as month-day
  ungroup()

######################## Modifying biotic variables needed for the model approach ############################

##Setting senescence threshold at 11 (Allison paper!)
age_senesce <- 11

##Setting threshold in data
intrinsic_variables <- intrinsic_variables %>%
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old"))) %>%
  mutate(age10 = (AgeYears - age_senesce) / 10) #scaled numeric version of Age centered at senescence threshold

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

########################### Wave energy and tide data processing #############################

##first, see what the earliest and latest dates are in the breeding season
##this helps with selecting dates when uploading the wave and tide data
dates <- metadata %>%
  group_by(season) %>%
  summarize(earliest_date = min(date, na.rm = TRUE), latest_date = max(date, na.rm = TRUE))

##list only CSV files containing "tide" in the filename
tide_files <- list.files(path = "./RawData/", pattern = "tide.*\\.csv$", full.names = TRUE)

##read and combine only those tide data files
tide_data_all <- tide_files %>%
  map_dfr(~ read.csv(.x, check.names = FALSE))

##fit tide data to specific breeding seasons
tide_data_all <- tide_data_all %>%
  mutate(Date = as.Date(as.character(Date)),
         MM = lubridate::month(Date),
         YY = lubridate::year(Date),
         season = case_when(MM == 12 ~ YY + 1,
                            MM %in% c(1, 2, 3) ~ YY))

##select necessary columns from the tide data
tide_data_clean <- tide_data_all %>%
  select("Date", "Time (GMT)", "Verified (ft)", "season")

##list only CSV files containing "wave" in the filename
wave_files <- list.files(path = "./RawData/", pattern = "wave.*\\.csv$", full.names = TRUE)

##read and combine only the wave data files
wave_data_all <- wave_files %>%
  map_dfr(~ read.csv(.x, colClasses = "character")) #everything as character makes it easier to join for now

##convert all wave columns to numeric
wave_data_all <- wave_data_all %>%
  mutate(YY = as.numeric(trimws(YY)),
         MM = as.numeric(trimws(MM)),
         DD = as.numeric(trimws(DD)),
         hh = as.numeric(trimws(hh)),
         WVHT = as.numeric(trimws(WVHT)),
         DPD = as.numeric(trimws(DPD))) %>%
  filter(WVHT != 99, DPD != 99) #99 is used as an NA placeholder in the data, so we can filter these out

##fit wave data to specific breeding seasons
wave_data_all <- wave_data_all %>%
  mutate(season = case_when(MM == 12 ~ YY + 1, #December data assigned to next year's breeding season
                            MM %in% c(1,2,3) ~ YY, #Jan-Mar data stay in current breeding season
                            TRUE ~ NA_real_)) #Other months not in any breeding season, assign NA

##suppress scientific notation
options(scipen = 999)

##select necessary columns from the wave data
wave_data_clean <- wave_data_all %>%
  select("YY", "MM", "DD", "hh", "WVHT", "DPD", "season") %>%
  mutate(wave_power = 0.49 * WVHT^2 * DPD) #calculate wave power from the wave power equation

##need date and time columns to match in tide and wave data
#cleaned tide data
tide_data_clean <- tide_data_clean %>%
  mutate(tide_datetime = as.POSIXct(paste(Date, `Time (GMT)`), #combine date + time
                                    format = "%Y-%m-%d %H:%M", #match format
                                    tz = "UTC")) #set time zone to UTC

#cleaned wave data
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

##set extreme wave and tide threshold levels
extreme_wave_threshold <- 30 #30 kW/m
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
  geom_line(linewidth = 1.2, color = "#1f78b4") +
  geom_point(color = "darkblue") +
  labs(title = "Extreme Tide and Storm Events (1996–2025)",
       x = "Year", y = "Number of Extreme Events") +
  theme_minimal()

##First we need the wave and tide datasets together with the proportion data!
abiotic_variables <- metadata %>%
  select(animalID, season, proportion, total_resights, AgeYears) %>%
  distinct() %>%
  filter(!is.na(proportion), proportion > 0) %>%
  left_join(harem_assignment, by = c("animalID", "season")) %>%
  left_join(tide_wave_flagged, by = "season") %>%
  mutate(season_fct = as.factor(season),
         animalID_fct = as.factor(animalID),
         area_fct = as.factor(area))

############################# Harem density data processing #################################

##seal density csv read
seal.density <- read.csv("./RawData/seal.density.csv")

##avg area density calculation
area_density <- seal.density %>%
  rename(area = Beach) %>%
  mutate(date = ymd(date), season = year(date)) %>%
  group_by(area, season) %>%
  summarize(avg_density = mean(density), .groups = "drop") %>%
  left_join(harem_assignment, by = c("season", "area")) %>% 
  filter(!is.na(animalID)) #only keep animals observed in the 2016-2025 dataset

##check what area density looks like per-season at each location
ggplot(area_density, aes(x = area, y = avg_density, fill = area)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ season, scales = "free_x") +
  labs(x = "Location",
       y = "Average Density",
       title = "Harem Density per Location by Season") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##check whether there are ties for assigned harems
harem_duplicates <- area_density %>%
  group_by(season, animalID) %>%
  filter(n() > 1) %>%
  arrange(animalID, season) #there are duplicates

##break ties by the harem with higher density
area_density <- area_density %>%
  group_by(season, animalID) %>%
  slice_max(avg_density, with_ties = FALSE) %>%
  ungroup()

##check locations in the data that don't match with the ones we have density for from the Ano map
unique(area_density$area)
unique(seal.density$Beach) #there are location mismatches

##remove all harem density NAs for now
area_density <- area_density %>%
  filter(!is.na(avg_density))

########################### Modifying abiotic variables needed for the model approach ################################

##check location names in the dataset
unique(abiotic_variables$area)

##create groups for South and North point
abiotic_variables <- abiotic_variables %>%
  mutate(group = case_when(
    grepl("^N", area_fct) | area_fct %in% c("BBNN", "BBNS", "BBN") ~ "NP",
    TRUE ~ "SP"))

##first, modify abiotic_variables so it has a column for proportion = 1 (1) OR not 1 (0)
abiotic_variables <- abiotic_variables %>%
  mutate(is_one = as.numeric(proportion == 1))

##create a version of abiotic_variables with proportion < 1
abiotic_variables_sub <- abiotic_variables %>%
  filter(proportion < 1)

##transform proportion for this filtered data to make it suitable for gamma distribution 
flipped_prop <- max(abiotic_variables_sub$proportion) - abiotic_variables_sub$proportion #subtract all proportions from the maximum
flipped_prop <- flipped_prop - min(flipped_prop) + 0.0000001 #ensure no 0s by adding small constant
abiotic_variables_sub$flipped_prop <- flipped_prop


######################### Plotting distributions and predictors #############################

##visualizing sample size is useful for determining stringency

##sample size per season
sample_size_per_season <- intrinsic_variables %>%
  group_by(season) %>%
  summarize(sample_size = n_distinct(animalID))

##create a plot for sample size per season
ggplot(data = sample_size_per_season, aes(x = season, y = sample_size)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(x = "Season", 
       y = "Sample size") +
  theme_few() + 
  scale_y_continuous(n.breaks = 10)

##sample size per age
sample_size_per_age <- intrinsic_variables %>%
  group_by(AgeYears) %>%
  summarize(sample_size = n_distinct(animalID))

##create a plot for sample size by age
ggplot(data = sample_size_per_age, aes(x = AgeYears, y = sample_size)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(x = "Age", 
       y = "Sample Size") +
  theme_few() +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 20)

##histogram of proportions
ggplot(data = intrinsic_variables, aes(x = proportion)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "darkblue") +
  labs(x = "Proportion MPA", 
       y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_few()

##histogram of flipped proportion
ggplot(data = intrinsic_variables_sub, aes(x = flipped_prop)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "darkblue") +
  labs(title = "Histogram of Flipped Proportion MPA <1",
       x = "Flipped Proportion MPA", y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_few()

##plot age against proportion
ggplot(data = intrinsic_variables, aes(x = AgeYears, y = proportion)) +
  geom_point() +
  labs(x = "Age", 
       y = "Proportion") +
  theme_few()

##plot year born against proportion
ggplot(data = intrinsic_variables, aes(x = year_born_num, y = proportion)) +
  geom_point() +
  labs(x = "Year Born", 
       y = "Proportion")

##plot season against proportion
ggplot(data = intrinsic_variables, aes(x = season_fct, y = proportion)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(x = "Season", 
       y = "Proportion")

##variation in proportion by yr born
intrinsic_summary <- intrinsic_variables %>%
  group_by(year_born) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE),
            conf_int = sd(proportion, na.rm = TRUE) / sqrt(n()))

##plot yr born against mean proportion with CIs
ggplot(intrinsic_summary, aes(x = year_born, y = mean_prop)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_prop - conf_int, ymax = mean_prop + conf_int), width = 0.2) +
  labs(x = "Cohort (Year Born)",
       y = "Mean Proportion Mom–Pup Association") +
  theme_minimal()
