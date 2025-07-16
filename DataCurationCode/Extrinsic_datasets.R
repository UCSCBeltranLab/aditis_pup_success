##Create files for wave energy and tides
library(tidyverse)
library(lubridate)
library(utils)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

#First, see what the earliest and latest dates are in the breeding season
dates <- metadata %>%
  group_by(season) %>%
  summarise(
    earliest_date = min(date, na.rm = TRUE),
    latest_date = max(date, na.rm = TRUE)
  )


#List only CSV files containing "tide" in the filename
tide_files <- list.files(pattern = "tide.*\\.csv$")

#Read and combine only the tide data files
tide_data_all <- tide_files %>%
  map_dfr(read_csv)

tide_data_all <- tide_data_all %>%
  mutate(
    Date = as.Date(Date),  # ensure Date is in Date format
    MM = lubridate::month(Date),
    YY = lubridate::year(Date),
    season = case_when(
      MM == 12 ~ YY + 1, # December belongs to next year's season
      MM %in% c(1, 2, 3) ~ YY))  # Jan–Mar stay in current year


tide_data_clean <- tide_data_all %>%
  select("Date", "Time (GMT)", "Verified (ft)", "season")

#List only CSV files containing "wave" in the filename
wave_files <- list.files(pattern = "wave.*\\.csv$")

#Read and combine only the wave data files
wave_data_all <- wave_files %>%
  map_dfr(~ read.csv(.x, colClasses = "character"))

wave_data_all <- wave_data_all %>%
  mutate(
    YY = as.numeric(trimws(YY)),
    MM = as.numeric(trimws(MM)),
    DD = as.numeric(trimws(DD)),
    hh = as.numeric(trimws(hh)),
    WVHT = as.numeric(trimws(WVHT)),
    DPD = as.numeric(trimws(DPD))) %>%
  filter(WVHT != 99, DPD != 99)


wave_data_all <- wave_data_all %>%
  mutate(
    season = case_when(
      MM == 12 ~ YY + 1,   # December data assigned to next year's season
      MM %in% c(1,2,3) ~ YY,  # Jan-Mar data stay in current year
      TRUE ~ NA_real_      # Other months not in any season, assign NA
      ))

options(scipen = 999)

wave_data_clean <- wave_data_all %>%
  select("YY", "MM", "DD", "hh", "WVHT", "DPD", "season") %>%
  mutate(wave_power = 0.49 * WVHT^2 * DPD)


tide_data_clean <- tide_data_clean %>%
  mutate(tide_datetime = as.POSIXct(paste(Date, `Time (GMT)`),
                               format = "%Y-%m-%d %H:%M:%S",
                               tz = "UTC"))

wave_data_clean <- wave_data_clean %>%
mutate(
  wave_datetime = as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00", YY, MM, DD, hh),
                           format = "%Y-%m-%d %H:%M:%S",
                           tz = "UTC"))

tide_wave <- left_join(wave_data_clean, tide_data_clean, by = c("wave_datetime" = "tide_datetime", "season")) %>%
filter(!is.na(`Verified (ft)`))
