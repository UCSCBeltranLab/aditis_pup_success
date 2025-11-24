
##run the source code from Data Processing
source("./DataCurationCode/Data_processing_MPA.R")

######################## Piecewise binomial (1996 - 2025) ##################################

### Piecewise model for is_one - full dataset ###
mod_binom_1996_2025 <- glmmTMB(is_one ~ age10 : age_cat + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                               family = binomial(link = "logit"),
                               data = intrinsic_variables)
summary(mod_binom_1996_2025) #model summary
exp(fixef(mod_binom_1996_2025)$cond) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_binom_1996_2025, plot = TRUE) #plot all residuals 
testDispersion(simulateResiduals(mod_binom_1996_2025)) #check for overdispersion
check_collinearity(mod_binom_1996_2025) #check predictor VIFs

### Visualizing and plotting the effects ###

# 1) A table for raw observed proportion to overlay on plot
observed_data <- intrinsic_variables %>%
  group_by(AgeYears) %>%
  summarize(n = n(),
            n_one = sum(is_one), 
            perc_one = n_one / n, #actual proportion of is_one = 1 per age
            lwr = binom.test(n_one, n)$conf.int[1], #CIs
            upr = binom.test(n_one, n)$conf.int[2])

# 2) Create ANNUAL (rand effect of season only) predictions
pred_season <- ggpredict(mod_binom_1996_2025, 
                         terms = c("age10 [all]", "age_cat", "season_fct"), 
                         type = "random") %>%
  as_tibble() %>%
  mutate(AgeYears = x * 10 + age_senesce,
         age_cat = factor(group, levels = c("Young", "Old")),
         season_fct = factor(facet)) %>%
  filter((age_cat == "Young" & AgeYears < age_senesce) | (age_cat == "Old" & AgeYears >= age_senesce))

# 3) Create a prediction grid to visualize the Age effect
pred_grid <- expand.grid(AgeYears = 3:22,
                         season_fct = levels(intrinsic_variables$season_fct),
                         animalID_fct = NA) %>% #placeholder since we won't predict per individual
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old")),
         age10 = (AgeYears - age_senesce) / 10,
         total_resights = mean(intrinsic_variables$total_resights, na.rm = TRUE),
         age_last_seen = mean(intrinsic_variables$age_last_seen, na.rm = TRUE)) #use means for predictors we want to hold constant

##identify number of observations per year
seal_years <- intrinsic_variables %>%
  count(season_fct)

# 4) Make predictions for age effect
pred_grid <- pred_grid %>%
  mutate(pred_logit = predict(mod_binom_1996_2025, newdata = pred_grid, re.form = NA), #linear predictor (logit scale) using fixed effects only
         se_logit = predict(mod_binom_1996_2025, newdata = pred_grid, re.form = NA, se.fit = TRUE)$se.fit) %>% #standard error of prediction on logit scale
  mutate(conf_lo = pred_logit - 1.96 * se_logit, 
         conf_hi = pred_logit + 1.96 * se_logit) %>% #Compute 95% confidence intervals (still on logit scale)
  left_join(seal_years, by = "season_fct") %>% #Add sample size weights for each season
  group_by(AgeYears, age_cat) %>%
  #Compute weighted summary statistics on probability scale
  summarise(predicted = weighted.mean(plogis(pred_logit), n), #avg predicted probability
            conf_lo = weighted.mean(plogis(conf_lo), n), #CI lower
            conf_hi = weighted.mean(plogis(conf_hi), n), #CI upper
            .groups = "drop")

# 5) Plot predictions with senescent threshold
ggplot(pred_grid, aes(x = AgeYears, y = predicted, color = age_cat, fill = age_cat)) +
  
  #Predictions for individual seasons
  geom_line(aes(group = interaction(age_cat, season_fct)),
            pred_season,
            alpha = 0.3, color = "#7EAAC1") +
  
  #Main ribbon for model w/ CI
  geom_ribbon(data = pred_grid,
              aes(x = AgeYears, ymin = conf_lo, ymax = conf_hi, fill = age_cat),
              alpha = 0.3, color = NA, inherit.aes = FALSE) +
  
  #Thick colored prediction line
  geom_line(data = pred_grid,
            aes(x = AgeYears, y = predicted, color = age_cat),
            linewidth = 1.2, inherit.aes = FALSE) +
  
  #Raw observed points
  geom_pointrange(data = observed_data,
                  aes(x = AgeYears, y = perc_one, ymin = lwr, ymax = upr),
                  inherit.aes = FALSE, color = "black") +
  
  #Overlay sample size for each age
  geom_text(data = observed_data,
            aes(x = AgeYears, y = 1.03, label = n),
            inherit.aes = FALSE, size = 3) +
  
  #Formatting, colors
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_y_continuous(n.breaks = 5) + 
  scale_x_continuous(n.breaks = 20) +
  scale_color_manual(values = c("Young" = "#66C2A5", "Old" = "#D5A5C9")) +
  scale_fill_manual(values = c("Young" = "#66C2A5", "Old" = "#D5A5C9")) +
  labs(x = "Age (Years)", 
       y = "Probability of 1", 
       color = "Age class", 
       fill = "Age class") +
  theme_few()

############################ Piecewise gamma model (1996-2025) ####################################

### Piecewise model for intermediate proportions - full dataset ###
mod_gamma_1996_2025 <- glmmTMB(flipped_prop ~ age10 : age_cat + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                               family = Gamma(link = "log"),
                               data = intrinsic_variables_sub)
summary(mod_gamma_1996_2025) #model summary
exp(fixef(mod_gamma_1996_2025)$cond) #gamma (log link) = log-scale coef → multiplicative effects
simulateResiduals(mod_gamma_1996_2025, plot = TRUE) #plot all residuals 
testDispersion(simulateResiduals(mod_gamma_1996_2025)) #check for overdispersion
check_collinearity(mod_gamma_1996_2025) #check predictor VIFs

##################################### Piecewise binomial (2016 - 2025) ########################################

# 1) Subset the data for only seasons 2016 - 2025
intrinsic_2016_2025 <- intrinsic_variables %>%
  filter(season >= 2016)

##Drop levels from the model with all years!
intrinsic_2016_2025 <- intrinsic_2016_2025 %>% 
  droplevels()

# 2) Add harem density data for 2016-2025
##Result: intrinsic data for 2016_2025 with average density per season/area
intrinsic_2016_2025 <- intrinsic_2016_2025 %>%
  left_join(area_density, by = c("animalID", "season"))

##Quick check: Histogram of filtered data's proportions
ggplot(data = intrinsic_2016_2025, aes(x = proportion)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "darkblue") +
  labs(x = "Proportion MPA", 
       y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_few()
  
### Piecewise model for is_one 2016-2025 ###
mod_binom_2016_2025 <- glmmTMB(is_one ~ age10 : age_cat + avg_density + age_last_seen + total_resights + (1 | season_fct) + (1 | animalID_fct),
                               family = binomial(link = "logit"),
                               data = intrinsic_2016_2025)
summary(mod_binom_2016_2025) #model summary
exp(fixef(mod_binom_2016_2025)$cond) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_binom_2016_2025, plot = TRUE) #plot all residuals 
testDispersion(simulateResiduals(mod_binom_2016_2025)) #check for overdispersion
check_collinearity(mod_binom_2016_2025) #check predictor VIFs

# 1) A table for raw observed proportion to overlay on plot
observed_data_2016_2025 <- intrinsic_2016_2025 %>%
  group_by(AgeYears) %>%
  summarize(n = n(),
            n_one = sum(is_one), 
            perc_one = n_one / n, #actual proportion of is_one = 1 per age
            lwr = binom.test(n_one, n)$conf.int[1], #CIs
            upr = binom.test(n_one, n)$conf.int[2])

# 2) Create ANNUAL (rand effect of season only) predictions
pred_season_2016_2025 <- ggpredict(mod_binom_2016_2025, 
                         terms = c("age10 [all]", "age_cat", "season_fct"), 
                         type = "random") %>%
  as_tibble() %>%
  mutate(AgeYears = x * 10 + age_senesce,
         age_cat = factor(group, levels = c("Young", "Old")),
         season_fct = factor(facet)) %>%
  filter((age_cat == "Young" & AgeYears < age_senesce) | (age_cat == "Old" & AgeYears >= age_senesce))

# 3) Create a prediction grid to visualize the Age effect
pred_grid_2016_2025 <- expand.grid(AgeYears = 3:22,
                         season_fct = levels(intrinsic_2016_2025$season_fct),
                         animalID_fct = NA) %>% #placeholder since we won't predict per individual
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old")),
         age10 = (AgeYears - age_senesce) / 10,
         total_resights = mean(intrinsic_2016_2025$total_resights, na.rm = TRUE),
         age_last_seen = mean(intrinsic_2016_2025$age_last_seen, na.rm = TRUE),
         avg_density = mean(intrinsic_2016_2025$avg_density, na.rm = TRUE)) #use means for predictors we want to hold constant

##identify number of observations per year
seal_years_2016_2025 <- intrinsic_2016_2025 %>%
  count(season_fct)

# 4) Make predictions for age effect
pred_grid_2016_2025 <- pred_grid_2016_2025 %>%
  mutate(pred_logit = predict(mod_binom_2016_2025, newdata = pred_grid_2016_2025, re.form = NA), #linear predictor (logit scale) using fixed effects only
         se_logit = predict(mod_binom_2016_2025, newdata = pred_grid_2016_2025, re.form = NA, se.fit = TRUE)$se.fit) %>% #standard error of prediction on logit scale
  mutate(conf_lo = pred_logit - 1.96 * se_logit, 
         conf_hi = pred_logit + 1.96 * se_logit) %>% #Compute 95% confidence intervals (still on logit scale)
  left_join(seal_years_2016_2025, by = "season_fct") %>% #Add sample size weights for each season
  group_by(AgeYears, age_cat) %>%
  #Compute weighted summary statistics on probability scale
  summarise(predicted = weighted.mean(plogis(pred_logit), n), #avg predicted probability
            conf_lo = weighted.mean(plogis(conf_lo), n), #CI lower
            conf_hi = weighted.mean(plogis(conf_hi), n), #CI upper
            .groups = "drop")

# 5) Plot predictions with senescent threshold
ggplot(pred_grid_2016_2025, aes(x = AgeYears, y = predicted, color = age_cat, fill = age_cat)) +
  
  #Predictions for individual seasons
  geom_line(aes(group = interaction(age_cat, season_fct)),
            pred_season_2016_2025,
            alpha = 0.3, color = "#7EAAC1") +
  
  #Main ribbon for model w/ CI
  geom_ribbon(data = pred_grid_2016_2025,
              aes(x = AgeYears, ymin = conf_lo, ymax = conf_hi, fill = age_cat),
              alpha = 0.3, color = NA, inherit.aes = FALSE) +
  
  #Thick colored prediction line
  geom_line(data = pred_grid_2016_2025,
            aes(x = AgeYears, y = predicted, color = age_cat),
            linewidth = 1.2, inherit.aes = FALSE) +
  
  #Raw observed points
  geom_pointrange(data = observed_data_2016_2025,
                  aes(x = AgeYears, y = perc_one, ymin = lwr, ymax = upr),
                  inherit.aes = FALSE, color = "black") +
  
  #Overlay sample size for each age
  geom_text(data = observed_data_2016_2025,
            aes(x = AgeYears, y = 1.03, label = n),
            inherit.aes = FALSE, size = 3) +
  
  #Formatting, colors
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_y_continuous(n.breaks = 5) + 
  scale_x_continuous(n.breaks = 20) +
  scale_color_manual(values = c("Young" = "#66C2A5", "Old" = "#D5A5C9")) +
  scale_fill_manual(values = c("Young" = "#66C2A5", "Old" = "#D5A5C9")) +
  labs(x = "Age (Years)", 
       y = "Probability of 1", 
       color = "Age class", 
       fill = "Age class") +
  theme_few()

############################# Piecewise gamma model (2016-2025) #####################################

##Subset the intermediate proportion data for years 2016-2025
intrinsic_sub_2016_2025 <- intrinsic_variables_sub %>%
  filter(season >= 2016) %>%
  left_join(area_density, by = c("animalID", "season")) #join with area density data 

### Piecewise model for intermediate proportions 2016-2025 ###
mod_gamma_2016_2025 <- glmmTMB(flipped_prop ~ age_cat : age10 + avg_density + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                               family = Gamma(link = "log"),
                               data = intrinsic_sub_2016_2025)
summary(mod_gamma_2016_2025) #model summary
exp(fixef(mod_gamma_2016_2025)$cond) #gamma (log link) = log-scale coef → multiplicative effects
simulateResiduals(mod_gamma_2016_2025, plot = TRUE) #plot all residuals 
testDispersion(simulateResiduals(mod_gamma_2016_2025)) #check for overdispersion
check_collinearity(mod_gamma_2016_2025) #check predictor VIFs

# 1) predicted effect of density on flipped proportion
pred_density <- ggpredict(mod_gamma_2016_2025, terms = "avg_density [all]")

# 2) plot flipped proportion and average density
ggplot(data = pred_density, aes(x = x, y = predicted)) +
  
  #Raw points
  geom_point(data = intrinsic_sub_2016_2025,
             aes(x = avg_density, y = flipped_prop),
             color = "darkblue",
             alpha = 0.4,
             size = 2,
             inherit.aes = FALSE) +
  
  #Confidence ribbon
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "skyblue",
              alpha = 0.3) +
  
  #Prediction line (inherits global aes)
  geom_line(color = "darkblue",
            linewidth = 1.2) +
  
  #Formatting
  scale_y_reverse() +
  labs(x = "Average Density",
       y = "Flipped Proportion") +
  theme_minimal(base_size = 14)

################# Abiotic binomial and gamma models (1996-2023) ##################

### Piecewise binomial model with per-year extreme wave/tide events ###
mod_abiotic_bin <- glmmTMB(is_one ~ n_extreme_both + group_fct + total_resights + (1 | season_fct) + (1 | animalID_fct),
                           data = abiotic_variables,
                           family = binomial(link = "logit"))
summary(mod_abiotic_bin) #model summary
exp(fixef(mod_abiotic_bin)$cond) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_abiotic_bin, plot = TRUE) #plot all residuals 
testDispersion(simulateResiduals(mod_abiotic_bin)) #check for overdispersion
check_collinearity(mod_abiotic_bin) #check predictor VIFs
 
### Plot effects of group (NP or SP) on proportion MPA ###

# 1) predicted effect of group on is_one
effect_group <- ggpredict(mod_abiotic_bin, terms = "group_fct")

# 2) plot predicted probability of 1
ggplot(effect_group, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  geom_point(size = 3) +
  labs(x = "Group",
       y = "Predicted probability of 1") +
  scale_y_continuous(n.breaks = 8) +
  theme_minimal(base_size = 14)

### Piecewise gamma model with per-year extreme wave/tide events ###
mod_abiotic_gam <- glmmTMB(flipped_prop ~ n_extreme_both + group_fct + total_resights + (1 | season_fct) + (1 | animalID_fct),
                           data = abiotic_variables_sub,
                           family = Gamma(link = "log"))
summary(mod_abiotic_gam) #model summary
exp(fixef(mod_abiotic_gam)$cond) #gamma (log link) = log-scale coef → multiplicative effects
simulateResiduals(mod_abiotic_gam, plot = TRUE) #plot all residuals 
testDispersion(simulateResiduals(mod_abiotic_gam)) #check for overdispersion
check_collinearity(mod_abiotic_gam) #check predictor VIFs

##################### Consistency models and figures ###################

##Subset by age to avoid need to account for age/specific piecewise split 
#For all data
young_data <- intrinsic_variables %>% filter(age_cat == "Young") 
old_data   <- intrinsic_variables %>% filter(age_cat == "Old")

#For 0 - 1 (gamma model) data
young_data_sub <- intrinsic_variables_sub %>% filter(age_cat == "Young")
old_data_sub   <- intrinsic_variables_sub %>% filter(age_cat == "Old")

### Population-level repeatability calculations ###

#Repeatability for ALL individuals from all years, for is_one binomial 
rpt_all <- rptBinary( #Use Binary version of rpt
  is_one ~ AgeYears + age_last_seen + total_resights +
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = intrinsic_variables,
  link     = "logit",
  nboot    = 1000,
  npermut  = 0
); summary(rpt_all) #summarize results

#Repeatability for young individuals from all years, for is_one binomial 
rpt_young <- rptBinary( #Use Binary version of rpt
  is_one ~ age10 + age_last_seen + total_resights + 
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = young_data,
  link     = "logit",
  nboot    = 1000,
  npermut  = 0
); summary(rpt_young) #summarize results

#Repeatability for old individuals from all years, for is_one binomial 
rpt_old <- rptBinary( #Use Binary version of rpt
  is_one ~ age10 + age_last_seen + total_resights + 
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = old_data,
  link     = "logit",
  nboot    = 1000,
  npermut  = 0
); summary(rpt_old) #summarize results

#Repeatability for ALL individuals, all years, for between 0 and 1 
rpt_sub_all <- rpt( #Use regular rpt that works for normally distributed data
  flipped_prop ~ age10 + age_last_seen + total_resights +
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = intrinsic_variables_sub,
  datatype = "Gamma",
  nboot    = 1000,
  npermut  = 0
); summary(rpt_sub_all) #summarize results

#Repeatability for young individuals, all years, for between 0 and 1 
rpt_young_sub <- rpt( #Use regular rpt that works for normally distributed data
  flipped_prop ~ age10 + age_last_seen + total_resights +
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = young_data_sub,
  datatype = "Gamma",
  nboot    = 1000,
  npermut  = 0
); summary(rpt_young_sub) #summarize results

#Repeatability for old individuals, all years, for between 0 and 1 
rpt_old_sub <- rpt( #Use regular rpt that works for normally distributed data
  flipped_prop ~ age10 + age_last_seen + total_resights +
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = old_data_sub,
  datatype = "Gamma",
  nboot    = 1000,
  npermut  = 0
); summary(rpt_old_sub) #summarize results

### Individual level consistency ###

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
  filter(n_seasons >= 2) #require at least two repeated breeding seasons per individual

# 3) Plot consistency across entire dataset
ggplot(individual_consistency,
       aes(x = mean_assoc, y = consistency)) +
  geom_point(aes(size = n_obs_total), alpha = 0.75) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey60") +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_size_continuous(name = "Total observations") + #This is TOTAL observations across all breeding seasons (not seasons)
  labs(
    x = "Mean maternal association (0–1)",
    y = "Consistency (1 − SD/0.5)",
    title = "Individual-level consistency"
  ) +
  theme_minimal(base_size = 14)

# 4) Classify individuals into behavioral strategies 
##Strategy classification:
# Maternal specialist: consistently associated (mean > 0.8, consistency > 0.8)
# Non-maternal specialist: consistently unassociated (mean < 0.2, consistency > 0.8)
# Plastic: fluctuates in association (consistency < 0.8)
# Intermediate specialist: edge cases with moderate consistency near thresholds

individual_consistency <- individual_consistency %>%
  mutate(strategy = case_when(
    mean_assoc <= 0.2 & consistency >= 0.8 ~ "Non-maternal specialist",
    mean_assoc >= 0.8 & consistency >= 0.8 ~ "Maternal specialist",
    consistency < 0.8 ~ "Plastic",
    TRUE ~ "Intermediate specialist"))

# 5) Plot these strategies with total observations (not season!) in consistency / MPA space 
ggplot(individual_consistency,
       aes(x = mean_assoc, y = consistency, color = strategy)) +
  geom_point(aes(size = n_obs_total, alpha = 0.8)) + #size of points weighted by number of observations
  geom_vline(xintercept = c(0.2, 0.8), linetype = "dashed", color = "grey60") + #MPA cutoff for behavioral strategy
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey60") + #consistency cutoff or behavioral strategy
  scale_color_manual(values = c(
    "Maternal specialist" = "#1f78b4",
    "Non-maternal specialist" = "#B74F6F",
    "Plastic" = "#6F73D2",
    "Intermediate specialist" = "#66C2A5")) +
  scale_size_continuous(name = "Total observations", range = c(2,10)) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x = "Mean maternal association (0–1)",
    y = "Consistency (0–1)",
    color = "Behavioral strategy",
    title = "Individual consistency and maternal behavior strategy") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"))

### Foundational consistency figures (allows us to look at specific animalIDs) ###

##make a data table to create a star label for animals with consistent prop == 1
star_data <- individual_consistency %>%
  group_by(animalID_fct) %>%
  summarize(all_one = all(sd_assoc == 0),
            max_y = max(proportion)) %>%
  filter(all_one) %>%
  mutate(y_position = max_y + 0.05)

## boxplot with color gradient for proportion and standard deviation 
ggplot(individual_consistency, aes(x = animalID_fct, y = proportion, fill = sd_assoc)) +
  geom_boxplot(width = 0.7, outlier.size = 1.5, outlier.color = "black") +
  geom_jitter(size = 1, alpha = 0.6, color = "#04BBB2") +
  scale_fill_gradient(low = "#D8FFF7", high = "#073481", name = "SD of Proportion") +
  geom_text(data = star_data, aes(x = animalID_fct, y = y_position, label = "*"), color = "#C699E1", size = 7, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 1.03)) +
  scale_y_continuous(n.breaks = 10) +
  theme_few() +
  labs(x = "Animal ID", 
       y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, size = 6))

# 1) one row per animal with its consistency
animal_order <- individual_consistency %>%
  distinct(animalID_fct, consistency) %>%
  arrange(desc(consistency)) %>%
  pull(animalID_fct)

# 2) order animalID to consistency score
individual_consistency <- individual_consistency %>%
  mutate(animalID_fct = factor(animalID_fct, levels = animal_order))

# 3) plot same as above but ordered from most -> least consistent
ggplot(individual_consistency, aes(x = animalID_fct, y = proportion, fill = consistency)) +
  geom_boxplot(width = 0.7, outlier.size = 1.5, outlier.color = "black") +
  geom_jitter(size = 1, alpha = 0.6, color = "#04BBB2") +
  scale_fill_gradient(low = "#D8FFF7", high = "#073481",
                      name = "Individual consistency score") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(n.breaks = 10) +
  theme_few() +
  labs(x = "Animal ID (ordered by consistency)", 
       y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, size = 6))
