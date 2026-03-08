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

################## Colors for season effects (age and density) ######################

### density ###

season_pal <- c("#F4A3A3","#F6C177", "#F9E2AF", "#A6D189","#81C8BE", "#8CAAEE","#9CC7E4", "#F5BDE6")
names(season_pal) <- levels(intrinsic_2016_2023$season_fct)

plot_density_by_season <- ggplot(intrinsic_2016_2023,
                                 aes(avg_density, proportion, color = season_fct)) +
  geom_point(
    alpha = 0.45, size = 2,
    position = position_jitter(height = 0.008, width = 0)
  ) +
  geom_ribbon(
    data = pred_density_2016_2023,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "navy", alpha = 0.18, inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_density_2016_2023,
    aes(x = x, y = predicted),
    color = "navy", linewidth = 1.2, inherit.aes = FALSE
  ) +
  scale_color_manual(values = season_pal, name = "Season") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_few() +
  labs(
    x = "Within-colony location density",
    y = "Mother-offspring association"
  )

plot_density_by_season

# age with season colored
plot_age_season_col_2016_2023 <- ggplot() +
  geom_line(
    data = season_lines_2016_2023,
    aes(AgeYears, pred, group = season_fct, color = season_fct),
    alpha = 0.45
  ) +
  geom_ribbon(
    data = ci_2016_2023,
    aes(AgeYears, ymin = lo, ymax = hi),
    fill = "#83C05A", alpha = 0.28
  ) +
  geom_line(
    data = ci_2016_2023,
    aes(AgeYears, pred),
    color = "#83C05A", linewidth = 1.2
  ) +
  geom_pointrange(
    data = observed_data_2016_2023,
    aes(AgeYears, avg_prop, ymin = lwr, ymax = upr),
    color = "black"
  ) +
  geom_text(
    data = observed_data_2016_2023,
    aes(AgeYears, 1.01, label = n_age),
    vjust = -0.5
  ) +
  scale_color_manual(values = season_pal, name = "Season") +
  scale_x_continuous(breaks = ages) +
  coord_cartesian(ylim = c(0.75, 1.02), clip = "off") +
  theme_few() +
  labs(x = "Age (Years)", y = "Mother-offspring association")

plot_age_season_col_2016_2023

#age with season colored for piecewise


############### Extreme event model plot with color jitter by season ############

# 1) predicted effect of extreme events by season
pred_extreme_2016_2023 <- ggpredict(mod_binom_2016_2023,
                                    terms = c("n_extreme_both [all]", "season_fct"))

# 2) force season order in the raw data
intrinsic_2016_2023 <- intrinsic_2016_2023 %>%
  mutate(season_fct = factor(season_fct,
                             levels = sort(unique(as.numeric(as.character(season_fct))))))

# 3) force the same order in ggpredict output (its season column is `group`)
pred_extreme_2016_2023$group <- factor(pred_extreme_2016_2023$group, levels = seasons)

# 5) plot: points + ribbons + lines, with ordered legend
ggplot(pred_extreme_2016_2023, aes(x = x, y = predicted, group = group, color = group, fill = group)) +
  geom_jitter(data = intrinsic_2016_2023,
              aes(x = n_extreme_both, y = proportion, color = season_fct),
              alpha = 0.6, size = 2, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.12, color = NA, fill = "#A30B37") +
  geom_line(linewidth = 1.2, color = "#A30B37") +
  scale_color_manual(values = pal, breaks = seasons, name = "Year") +
  scale_fill_manual(values = pal, breaks = seasons, name = "Year") +
  labs(x = "Number of extreme wave and tide events",
       y = "Mother-offspring association") +
  theme_few()
