##run the source code from Data Processing
source("./DataCurationCode/Data_processing_MPA.R")

##run the source code from GLMMs
source("./AnalysisCode/Final_GLMMs.R")

########## Sample sizes + proportion distribution #########

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

########### Raw plots: predictors vs. mother-offspring association ##################
##plot age (1996-2025) against proportion
ggplot(data = intrinsic_variables, aes(x = AgeYears, y = proportion)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Age", 
       y = "MOA") +
  theme_minimal()

##violin plot age vs. proportion
ggplot(data = intrinsic_variables %>%
         mutate(age_bin = factor(AgeYears)),
       aes(x = age_bin, y = proportion)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
  labs(x = "Age", 
       y = "MOA") +
  theme_minimal()

##subset (2016-2023) age against proportion
ggplot(data = intrinsic_2016_2023, aes(x = AgeYears, y = proportion)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Age",
       y = "MOA") +
  theme_minimal()

##plot season against proportion
ggplot(intrinsic_variables, aes(x = season_fct, y = proportion)) +
  geom_violin(fill = "lightblue", alpha = 0.6, scale = "width") +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
  labs(x = "Season", 
       y = "MOA") +
  theme_minimal()

season_means <- intrinsic_variables %>%
  group_by(season_fct) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE),
            .groups = "drop")

##plot avg_density against proportion
ggplot(intrinsic_variables,
       aes(x = avg_density, y = proportion)) +
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.6) +
  labs(x = "Average location density",
       y = "MOA",
       color = "Age") +
  theme_minimal()

##plot n_extreme_both against proportion
ggplot(intrinsic_variables,
       aes(x = n_extreme_both, y = proportion)) +
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.6) +
  labs(x = "Number of extreme wave and tide events",
       y = "MOA",
       color = "Age") +
  theme_minimal()

##change in extreme wave and tide events across years
ggplot(tide_wave_flagged, aes(x = season, y = n_extreme_both)) +
  geom_line(linewidth = 1.2, color = "#1f78b4") +
  geom_point(color = "darkblue") +
  labs(x = "Year", 
       y = "Number of Extreme Wave and Tide Events") +
  theme_minimal()

##per-year extreme wave and tide events bar plot
ggplot(tide_wave_flagged, aes(x = season, y = n_extreme_both)) +
  geom_bar(stat = "identity", fill = "lightpink") +
  labs(x = "Year", 
       y = "Number of extreme wave and tide events") +
  theme_few() +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 23)

##proportion vs. extreme wave and tide
ggplot(intrinsic_variables,
       aes(x = n_extreme_both, y = proportion)) +
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.6) +
  labs(x = "Number of extreme wave and tide events",
       y = "MOA",
       color = "Age") +
  theme_minimal()

##proportion for n_extreme both
ggplot(intrinsic_variables, aes(x = n_extreme_both, y = proportion)) +
  stat_summary(fun.data = mean_cl_boot,
               geom = "pointrange",
               size = 0.8) +
  stat_summary(fun = mean,
               geom = "line",
               aes(group = 1),
               linewidth = 0.8) +
  stat_summary(fun = mean,
               geom = "point",
               size = 2) +
  labs(x = "Number of extreme wave and tide events", 
       y = "MOA") +
  theme_minimal()

#create a summary for geom_text
season_summary <- intrinsic_variables %>%
  group_by(season, season_fct) %>%
  summarize(mean_prop = mean(proportion, na.rm = TRUE),
            n_extreme_both = first(n_extreme_both),
            .groups = "drop")

#plot season by proportion with labels above for n_extreme_both for each year
ggplot(intrinsic_variables, aes(x = season_fct, y = proportion)) +
  stat_summary(fun.data = mean_cl_boot,
               geom = "pointrange",
               size = 0.8) +
  stat_summary(fun = mean,
               geom = "line",
               aes(group = 1),
               linewidth = 0.8) +
  stat_summary(fun = mean,
               geom = "point",
               size = 2) +
  geom_text(data = season_summary,
            aes(x = season_fct, y = mean_prop, label = n_extreme_both),
            color = "black",
            nudge_y = 0.03,
            size = 3) +
  labs(x = "Season",
       y = "MOA") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


############## Raw data plots for Marm ############

#age colored by season (raw)
ggplot(intrinsic_variables, aes(x = AgeYears, y = proportion, color = season_fct)) +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  labs(x = "Age",
       y = "MOA",
       color = "Season") +
  theme_minimal()

#density colored by season (raw)
ggplot(intrinsic_variables, aes(x = avg_density, y = proportion, color = season_fct)) +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  labs(x = "Density",
       y = "MOA",
       color = "Season") +
  theme_minimal()

#extreme events w/ season colors
ggplot(intrinsic_variables, 
       aes(x = n_extreme_both, y = proportion, color = season_fct)) +
  geom_point(alpha = 0.35) +
  geom_smooth(aes(group = 1),
              method = "lm", 
              se = TRUE, 
              color = "black", 
              linewidth = 1) +
  labs(x = "Extreme events",
       y = "MOA",
       color = "Season") +
  theme_minimal()

#age faceted by season (raw)
ggplot(intrinsic_variables, aes(x = AgeYears, y = proportion)) +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  facet_wrap(~ season_fct) +
  labs(x = "Age",
       y = "MOA") +
  theme_minimal()

#density faceted by season (raw)
ggplot(intrinsic_variables, aes(x = avg_density, y = proportion)) +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  facet_wrap(~ season_fct) +
  labs(x = "Density",
       y = "MOA") +
  theme_minimal()

#support for threshold using raw data (1996-2025)
ggplot(intrinsic_variables, aes(x = AgeYears, y = proportion, color = age_cat)) +
  geom_point(alpha = 0.45) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_vline(xintercept = 9, linetype = "dashed") +
  labs(x = "Age",
       y = "MOA",
       color = "Age class") +
  theme_minimal()

#support for linear using raw data (2016-2023)
ggplot(intrinsic_2016_2023, aes(x = AgeYears, y = proportion, color = age_cat)) +
  geom_point(alpha = 0.45) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_vline(xintercept = 9, linetype = "dashed") +
  labs(x = "Age",
       y = "MOA",
       color = "Age class") +
  theme_minimal()

############## raw data support for interactions ###############

#interaction with age and density using age_cat
ggplot(intrinsic_variables, aes(x = avg_density, y = proportion)) +
  geom_point(alpha = 0.45) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  facet_wrap(~ age_cat) +
  labs(x = "Average location density", y = "MOA") +
  theme_minimal()

#interaction with age and extreme events using age_cat
ggplot(intrinsic_variables, aes(x = n_extreme_both, y = proportion)) +
  geom_point(alpha = 0.45) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  facet_wrap(~ age_cat) +
  labs(x = "Number of extreme wave and tide events", y = "MOA") +
  theme_minimal()

#with linear age as a color for each line
ggplot(intrinsic_variables, 
       aes(x = avg_density, y = proportion, color = factor(AgeYears))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Average location density",
       y = "MOA",
       color = "Age (years)") +
  theme_minimal()

#with linear age as a color for each line
ggplot(intrinsic_variables, 
       aes(x = n_extreme_both, y = proportion, color = factor(AgeYears))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Number of extreme wave and tide events",
       y = "MOA",
       color = "Age (years)") +
  theme_minimal()

#interaction age and density with linear age
ggplot(intrinsic_variables, aes(x = avg_density, y = proportion)) +
  geom_point(alpha = 0.45) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  facet_wrap(~ AgeYears) +
  labs(x = "Average location density", y = "MOA") +
  theme_minimal()

#interaction with age and extreme with linear age
ggplot(intrinsic_variables, aes(x = n_extreme_both, y = proportion)) +
  geom_point(alpha = 0.45) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  facet_wrap(~ AgeYears) +
  labs(x = "Number of extreme wave and tide events", y = "MOA") +
  theme_minimal()








