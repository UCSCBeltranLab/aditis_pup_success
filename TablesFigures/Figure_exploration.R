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

################## Plots for Marm: exploring random effect of season ######################

intrinsic_variables <- intrinsic_variables %>%
  mutate(season_fct = factor(season, levels = sort(unique(season))))

#age faceted by season
ggplot(intrinsic_variables,
       aes(x = AgeYears, y = proportion, color = season_fct)) +
  geom_point(alpha = 0.6) +
  geom_smooth(se = FALSE) +
  facet_wrap(~ season_fct) +
  labs(x = "Age (years)", 
       y = "MOA", 
       color = "Season") +
  theme_minimal()

#age colored by season
ggplot(intrinsic_variables,
       aes(x = AgeYears, y = proportion, color = season_fct)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Maternal age (years)", 
       y = "MOA") +
  theme_minimal()

#avg density colored by season
ggplot(intrinsic_variables,
       aes(x = avg_density, y = proportion, color = season_fct)) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Average harem density",
       y = "MOA",
       color = "Season") +
  theme_minimal()

#number of extreme wave and tide events colored by season
ggplot(intrinsic_variables,
       aes(x = n_extreme_both, y = proportion, color = season)) +  
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_viridis_c(option = "mako") +
  labs(x = "Number of extreme wave + tide events", 
       y = "MOA", 
       color = "Season") +
  theme_minimal()

############### Individual level consistency figures (individual ID effects) ###############

# 1) Build a dataframe pre-merge for modelling to calculate individual-level consistency
season_level <- Proportion_MPA %>% #has animalID, season, count_1_pup, total_resights
  mutate(proportion = count_1_pup / total_resights) %>% #Includes ALL data
  left_join(metadata %>% select(animalID, season, AgeYears),
            by = c("animalID","season")) %>%
  distinct() %>% #reduces to 1 animalID/season combo per observation (no duplicates)
  mutate(animalID_fct = factor(animalID),
         season_fct = factor(season))

#helper for weighted variance calculations
wvar <- function(x, w) {
  w <- w / sum(w) # (1) Normalize weights so they sum to 1
  mu <- sum(w * x) # (2) Compute the weighted mean of x
  sum(w * (x - mu)^2) # (3) Compute weighted average squared deviation from mean
}

# 2) Calculate cosnsitency score for MPA for each individual
##Summarize per individual:
#  - n_seasons: number of seasons observed
#  - n_obs_total: total number of resights
#  - mean_assoc: weighted mean proportion of association (0–1)
#  - sd_assoc: weighted SD across seasons, indicating behavioral variability
#  - frac_zero / frac_one: fraction of seasons where proportion = 0 or 1
#  - consistency: rescaled measure (1 − SD/0.5), bounded between 0–1 where 1 = perfectly consistent, 0 = highly variable
individual_consistency <- season_level %>%
  group_by(animalID_fct) %>%
  mutate(n_seasons   = n(),
         n_obs_total = sum(total_resights),
         mean_assoc  = weighted.mean(proportion, w = total_resights, na.rm = TRUE),     
         sd_assoc    = sqrt(wvar(proportion, w = total_resights)),     
         frac_zero   = mean(proportion == 0),
         frac_one    = mean(proportion == 1),
         .groups = "drop") %>%
  mutate(
    #Consistency: standardized so that SD = 0 (no variation) → 1, and SD = 0.5 (max) → 0
    consistency = pmax(pmin(1 - (sd_assoc / 0.5), 1), 0)   # 0–1, higher = more consistent
  ) %>%
  filter(n_seasons >= 3) #require at least 3 breeding seasons per individual

### Foundational consistency figures (allows us to look at specific animalIDs) ###

##make a data table to create a star label for animals with consistent prop == 1
star_data <- individual_consistency %>%
  group_by(animalID_fct) %>%
  summarize(all_one = all(sd_assoc == 0),
            max_y = max(proportion)) %>%
  filter(all_one) %>%
  mutate(y_position = max_y + 0.05)

## boxplot with color gradient for standard deviation + stars for "highly consistent" MOA indiciduals
ggplot(individual_consistency, aes(x = animalID_fct, y = proportion, fill = sd_assoc)) +
  geom_boxplot(width = 0.7, outlier.size = 0.6, outlier.color = "#04BBB2") +
  geom_jitter(size = 1, alpha = 0.6, color = "#04BBB2") +
  scale_fill_gradient(low = "#D8FFF7", high = "#073481", name = "SD of Proportion") +
  geom_text(data = star_data, aes(x = animalID_fct, y = y_position, label = "*"), color = "#C699E1", size = 7, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 1.03)) +
  scale_y_continuous(n.breaks = 10) +
  theme_few() +
  labs(x = "Individual ID", 
       y = "Mother-offspring association") +
  theme(axis.text.x = element_text(angle = 90, size = 6))

############ Age-corrected consistency figure ###########

## 1) Compute age-specific means and deviations (season level)
#   For each age (AgeYears), compute the mean maternal association across all seals,
#   then calculate each observation's deviation from that age-specific mean.
season_level_age_centered <- individual_consistency %>%
  group_by(AgeYears) %>%
  mutate(
    mean_prop_age = mean(proportion, na.rm = TRUE),          # population mean at this age
    dev_from_age  = proportion - mean_prop_age               # deviation from age-specific mean
  ) %>%
  ungroup()

# Visualize within-age variation in maternal association w/ boxplots for each age
ggplot(season_level_age_centered,
       aes(x = factor(AgeYears), y = dev_from_age)) +
  geom_boxplot(fill = "#A1D3B2", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +           # 0 = age-specific mean
  labs(
    x = "Age (years)",
    y = "Deviation from age-specific mean",
    title = "Within-age variation in maternal association"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))


## 1) Mean age per individual (weighted by resights)
#   Compute a weighted mean age per individual, weighting by total resights,
#   so individuals with more resight effort contribute more to their mean age.
age_by_ind <- season_level_age_centered %>%
  group_by(animalID_fct) %>%
  summarise(
    mean_age = weighted.mean(AgeYears, w = total_resights, na.rm = TRUE),
    .groups = "drop"
  )

## 2) Classify direction and consistency of deviation
#   dev_eps: threshold for how far from 0 a mean deviation must be to count as
#   "meaningfully above" or "meaningfully below" age peers.
dev_eps      <- 0.02   # deviation > 0.02 or < -0.02 is considered meaningful

# cons_thresh: threshold for high consistency in deviation across years.
# (e.g., consistency_dev >= 0.8)
cons_thresh  <- 0.7   

ind_consistency_age <- season_level_age_centered %>%
  group_by(animalID_fct) %>%
  summarize(
    n_seasons   = n_distinct(season),
    n_obs_total = sum(total_resights, na.rm = TRUE),
    mean_assoc  = mean(proportion, na.rm = TRUE),
    mean_dev    = mean(dev_from_age, na.rm = TRUE),
    sd_dev      = sd(dev_from_age, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Consistency metric (0..1): high when sd_dev is small relative to signal
    consistency_dev = 1 - (sd_dev / (sd_dev + abs(mean_dev) + 1e-6))
  )

# Using the individual-level object ind_consistency_age, which already contains:
# - mean_dev: mean of dev_from_age across seasons
# - consistency_dev: consistency in that deviation across seasons
# - mean_assoc: lifetime mean maternal association

ind_consistency_age <- ind_consistency_age %>%
  mutate(
    # Direction of mean deviation relative to age peers
    dev_dir = case_when(
      mean_dev >  dev_eps  ~ "Above peers",
      mean_dev < -dev_eps  ~ "Below peers",
      TRUE                 ~ "Near age mean"
    ),
    
    # Consistency of that deviation across years
    dev_consistency_class = case_when(
      dev_dir == "Above peers"    & consistency_dev >= cons_thresh ~ "Consistently above peers",
      dev_dir == "Below peers"    & consistency_dev >= cons_thresh ~ "Consistently below peers",
      dev_dir == "Near age mean"  & consistency_dev >= cons_thresh ~ "Consistently near mean",
      TRUE ~ "Inconsistent"  # fluctuates around peers, no stable position
    ),
    
    # Combine direction + consistency into a single age-normalized strategy label
    strategy_cross = case_when(
      dev_consistency_class == "Consistently above peers" ~ "Consistently above peers",
      dev_dir == "Above peers"    & dev_consistency_class == "Inconsistent" ~ "Inconsistent above peers",
      dev_dir == "Below peers"    & dev_consistency_class == "Inconsistent" ~ "Inconsistent below peers",
      dev_dir == "Near age mean"  & dev_consistency_class == "Inconsistent" ~ "Inconsistent near mean",
      TRUE ~ "Other"  # catch-all; should ideally be rare or empty
    ),
    
    # Order strategy levels for plotting (top: specialists, bottom: below peers)
    strategy_cross = factor(
      strategy_cross,
      levels = c(
        "Consistently above peers",
        "Inconsistent above peers",
        "Inconsistent near mean",
        "Inconsistent below peers"
      )
    )
  )


## 3) Summary statistics for strategy classes
# Overall summary: how many individuals in each consistency class, and their mean maternal association
summary_overall <- ind_consistency_age %>%
  group_by(dev_consistency_class) %>%
  summarise(
    n         = n(),                                 # number of individuals
    prop      = n / sum(n),                          # proportion of sample
    mean_assoc = mean(mean_assoc, na.rm = TRUE),     # average maternal association
    .groups   = "drop"
  ); summary_overall

# Direction-only summary: how many individuals are on average above, below, etc. and their mean MPA
summary_dir <- ind_consistency_age %>%
  group_by(dev_dir) %>%
  summarise(
    n         = n(),
    prop      = n / sum(n),
    mean_assoc = mean(mean_assoc, na.rm = TRUE),
    .groups   = "drop"
  ); summary_dir

# Cross-tab: direction x consistency class: how many individuals are in different strategies
summary_cross <- ind_consistency_age %>%
  group_by(dev_dir, dev_consistency_class) %>%
  summarise(
    n               = n(),
    prop_within_dir = n / sum(n),                    # proportion within each dev_dir
    mean_assoc      = mean(mean_assoc, na.rm = TRUE),
    .groups         = "drop"
  ); summary_cross


## 4) Join individual-level consistency back to season rows
individual_consistency_age_rows <- season_level_age_centered %>%
  left_join(
    ind_consistency_age %>%
      select(animalID_fct, consistency_dev),
    by = "animalID_fct"
  )

names(individual_consistency_age_rows)


## 5) Boxplot of within-individual variation colored by consistency
#   For each individual, show distribution of maternal association (proportion) across seasons
ggplot(individual_consistency_age_rows,
       aes(x = animalID_fct, y = proportion, fill = consistency_dev)) +
  geom_boxplot(width = 0.7, outlier.size = 1.5, outlier.color = "black") +
  geom_jitter(size = 1, alpha = 0.6, color = "#04BBB2") +
  scale_fill_gradient(
    low  = "#D8FFF7",
    high = "#073481",
    name = "Age-normalized\nconsistency"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(n.breaks = 10) +
  theme_few() +
  labs(
    x = "Animal ID",
    y = "Proportion maternal association"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = "right"
  )


## 6) Strategy scatterplot in "strategy space"
#   Create a compact dataframe for plotting strategy classes
ind_consistency_age_strat <- ind_consistency_age %>%
  select(animalID_fct, mean_assoc, consistency_dev, n_obs_total, strategy_cross)

# Plot individuals in "strategy space":
#   - x-axis: mean maternal association (0–1)
#   - y-axis: age-normalized consistency
#   - vertical dashed lines: thresholds for mean association
#   - horizontal dashed line: threshold for high consistency
ggplot(ind_consistency_age_strat,
       aes(x = mean_assoc, y = consistency_dev, color = strategy_cross)) +
  geom_point(aes(size = n_obs_total, alpha = 0.8)) +
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = cons_thresh, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c(
    "Consistently above peers"  = "#1f78b4",
    "Inconsistent above peers"  = "#6F73D2",
    "Inconsistent near mean"    = "#66C2A5",
    "Inconsistent below peers"  = "#B74F6F"
  )) +
  scale_size_continuous(name = "Total observations", range = c(2, 10)) +
  labs(
    x     = "Mean maternal association (0–1)",
    y     = "Consistency in deviation from age mean",
    color = "Age-normalized\nstrategy class",
    title = "Age-normalized maternal strategies",
    subtitle = "Deviation from age peers (x-axis) and consistency in that deviation (y-axis)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.box      = "vertical",
    legend.title    = element_text(face = "bold")
  )

## check who the most consistent individuals are!
ind_consistency_age$animalID_fct[ind_consistency_age$strategy_cross == "Consistently above peers"]

ind_consistency_age_perfect <- ind_consistency_age %>%
  filter(strategy_cross == "Consistently above peers")
