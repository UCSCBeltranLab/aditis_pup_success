
##run the source code from Data Processing
source("./DataCurationCode/Data_processing_MPA.R")

######################## Piecewise binomial (1996 - 2025) ##################################

##Piecewise model with animalID, season, year_born
mod_binom_1996_2025 <- glmmTMB(is_one ~ age10 : age_cat + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                               family = binomial(link = "logit"),
                               data = intrinsic_variables)
summary(mod_binom_1996_2025)
exp(fixef(mod_binom_1996_2025)$cond) 
DHARMa::simulateResiduals(mod_binom_1996_2025, plot = TRUE) 
DHARMa::testDispersion(DHARMa::simulateResiduals(mod_binom_1996_2025))
check_collinearity(mod_binom_1996_2025)

##Visualizing and plotting the effects!

##A table for raw observed proportion to overlay on plot
observed_data <- intrinsic_variables %>%
  group_by(AgeYears) %>%
  summarize(n = n(),
            n_one = sum(is_one),
            perc_one = n_one / n,
            lwr = binom.test(n_one, n)$conf.int[1],
            upr = binom.test(n_one, n)$conf.int[2])

##Create ANNUAL (i.e. rand effect of season only) predictions
pred_season <- ggpredict(mod_binom_1996_2025, 
  terms = c("age10 [all]", "age_cat", "season_fct"),
  type = "random") %>% 
  as_tibble() %>% 
  mutate(AgeYears = x * 10 + age_senesce,
         age_cat = factor(group, levels = c("Young", "Old")),
         season_fct = factor(facet)) %>%
  filter((age_cat == "Young" & AgeYears < age_senesce) | (age_cat == "Old" & AgeYears >= age_senesce))

##Need a prediction grid to visualize the Age effect
pred_grid <- expand.grid(AgeYears = 3:22,
  season_fct = levels(intrinsic_variables$season_fct),
  animalID_fct = NA) %>%  # placeholder since we won't predict per individual
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                     levels = c("Young", "Old")),
    age10 = (AgeYears - age_senesce) / 10,
    total_resights = mean(intrinsic_variables$total_resights, na.rm = TRUE),
    age_last_seen = mean(intrinsic_variables$age_last_seen, na.rm = TRUE),
    year_born_num = mean(intrinsic_variables$year_born_num, na.rm = TRUE))

seal_years <- intrinsic_variables %>%
  count(season_fct)

##Make predictions
pred_grid <- pred_grid %>%
  mutate(pred_logit = predict(mod_binom_1996_2025, newdata = pred_grid, re.form = NA),
         se_logit = predict(mod_binom_1996_2025, newdata = pred_grid, re.form = NA, se.fit = TRUE)$se.fit,
         pred_fixed = predict(mod_binom_1996_2025, newdata = pred_grid, re.form = NA)) %>%
  mutate(conf_lo = pred_logit - 1.96 * se_logit, 
         conf_hi = pred_logit + 1.96 * se_logit) %>%
  left_join(seal_years, by = "season_fct") %>%
  group_by(AgeYears, age_cat) %>%
  summarise(predicted = weighted.mean(plogis(pred_logit), n),
            predicted_pop = weighted.mean(plogis(pred_fixed), n),
            conf_lo = weighted.mean(plogis(conf_lo), n),
            conf_hi = weighted.mean(plogis(conf_hi), n),
            .groups = "drop")

##Plot predictions with senescent threshold
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

mod_gamma_1996_2025 <- glmmTMB(flipped_prop ~ age10 : age_cat + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                               family = Gamma(link = "log"),
                               data = intrinsic_variables_sub)
summary(mod_gamma_1996_2025)
exp(fixef(mod_gamma_1996_2025)$cond) 
DHARMa::simulateResiduals(mod_gamma_1996_2025, plot = TRUE)
DHARMa::testDispersion(DHARMa::simulateResiduals(mod_gamma_1996_2025))
check_collinearity(mod_gamma_1996_2025)

##################################### Piecewise binomial (2016 - 2025) ########################################

##Subset the data for only seasons 2016 - 2025
intrinsic_2016_2025 <- intrinsic_variables %>%
  filter(season >= 2016)

##Drop levels from the model with all years!
intrinsic_2016_2025 <- intrinsic_2016_2025 %>% droplevels()

##Add harem density data for 2016-2025
intrinsic_2016_2025 <- intrinsic_2016_2025 %>%
  left_join(area_density, by = c("animalID", "season"))

##Histogram of filtered data's proportions
ggplot(data = intrinsic_2016_2025, aes(x = proportion)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "darkblue") +
  labs(x = "Proportion MPA", 
       y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_few()
  
##Piecewise model for 2016-2025
mod_binom_2016_2025 <- glmmTMB(is_one ~ age10 : age_cat + avg_density + age_last_seen + total_resights + (1 | season_fct) + (1 | animalID_fct),
                               family = binomial(link = "logit"),
                               data = intrinsic_2016_2025)
summary(mod_binom_2016_2025)
exp(fixef(mod_binom_2016_2025)$cond)
DHARMa::simulateResiduals(mod_binom_2016_2025, plot = TRUE)
DHARMa::testDispersion(DHARMa::simulateResiduals(mod_binom_2016_2025))
check_collinearity(mod_binom_2016_2025)

##A table for raw observed proportion to overlay on plot
observed_data_2016_2025 <- intrinsic_2016_2025 %>%
  group_by(AgeYears) %>%
  summarize(n = n(),
            n_one = sum(is_one),
            perc_one = n_one / n,
            lwr = binom.test(n_one, n)$conf.int[1],
            upr = binom.test(n_one, n)$conf.int[2])

##Create ANNUAL (i.e. rand effect of season only) predictions
pred_season_2016_2025 <- ggpredict(
  mod_binom_2016_2025, 
  terms = c("age10 [all]", "age_cat", "season_fct"),
  type = "random") %>% 
  as_tibble() %>% 
  mutate(AgeYears = x * 10 + age_senesce,
         age_cat = factor(group, levels = c("Young", "Old")),
         season_fct = factor(facet),
         animalID_fct = factor(facet)) %>%
  filter((age_cat == "Young" & AgeYears < age_senesce) | (age_cat == "Old" & AgeYears >= age_senesce))

##Need a prediction grid to visualize the age effect
pred_grid_2016_2025 <- expand.grid(
  AgeYears = 3:22,
  season_fct = levels(intrinsic_2016_2025$season_fct),
  animalID_fct = NA) %>%
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old")),
         age10 = (AgeYears - age_senesce) / 10,
         total_resights = mean(intrinsic_2016_2025$total_resights, na.rm = TRUE),
         age_last_seen = mean(intrinsic_2016_2025$age_last_seen, na.rm = TRUE),
         avg_density = mean(intrinsic_2016_2025$avg_density, na.rm = TRUE))

seal_years_2016_2025 <- intrinsic_2016_2025 %>%
  count(season_fct)

##Make predictions
pred_grid_2016_2025 <- pred_grid_2016_2025 %>%
  mutate(pred_logit = predict(mod_binom_2016_2025, newdata = pred_grid_2016_2025, re.form = NULL),
         se_logit = predict(mod_binom_2016_2025, newdata = pred_grid_2016_2025, re.form = NA, se.fit = TRUE)$se.fit,
         pred_fixed = predict(mod_binom_2016_2025, newdata = pred_grid_2016_2025, re.form = NA)) %>%
  mutate(conf_lo = pred_logit - 1.96 * se_logit, 
         conf_hi = pred_logit + 1.96 * se_logit) %>%
  left_join(seal_years_2016_2025, by = "season_fct") %>%
  group_by(AgeYears, age_cat) %>%
  summarise(predicted = weighted.mean(plogis(pred_logit), n),
            predicted_pop = weighted.mean(plogis(pred_fixed), n),
            conf_lo = weighted.mean(plogis(conf_lo), n),
            conf_hi = weighted.mean(plogis(conf_hi), n),
            .groups = "drop")

##Plot predictions with senescent threshold
ggplot(pred_grid_2016_2025, aes(x = AgeYears, y = predicted, color = age_cat, fill = age_cat)) +
  
  #Predictions for individual seasons
  geom_line(aes(group = interaction(age_cat, season_fct)),
            pred_season_2016_2025, linewidth = 0.5,
            alpha = 0.9, color = "black") +
  
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

##Subset the data for years 2016-2025
intrinsic_sub_2016_2025 <- intrinsic_variables_sub %>%
  filter(season >= 2016) %>%
  left_join(area_density, by = c("animalID", "season"))

##Model with animalID and season
mod_gamma_2016_2025 <- glmmTMB(flipped_prop ~ age_cat : age10 + avg_density + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                               family = Gamma(link = "log"),
                               data = intrinsic_sub_2016_2025)
summary(mod_gamma_2016_2025)
exp(fixef(mod_gamma_2016_2025)$cond) 
DHARMa::simulateResiduals(mod_gamma_2016_2025, plot = TRUE)
DHARMa::testDispersion(DHARMa::simulateResiduals(mod_gamma_2016_2025))
check_collinearity(mod_gamma_2016_2025)

#must fit a simple beta model just to visualize density effect on proportion (unflipped)
library(betareg)

mod_beta <- betareg(proportion ~ avg_density,
                    data = intrinsic_sub_2016_2025)
summary(mod_beta)

pred_density <- ggpredict(mod_beta, terms = "avg_density [all]")

ggplot(data = pred_density, aes(x = x, y = predicted)) +
  
  #Raw points
  geom_point(data = intrinsic_sub_2016_2025,
             aes(x = avg_density, y = proportion),
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
  labs(x = "Average Density",
       y = "Proportion") +
  theme_minimal(base_size = 14)

##################### Consistency models and figures ###################

#Subset by age to avoid need to account for age/specific piecewise split 
#For all data
young_data <- intrinsic_variables %>% filter(age_cat == "Young") 
old_data   <- intrinsic_variables %>% filter(age_cat == "Old")

#For 0 - 1 (gamma model) data
young_data_sub <- intrinsic_variables_sub %>% filter(age_cat == "Young")
old_data_sub   <- intrinsic_variables_sub %>% filter(age_cat == "Old")

##Population-level repeatability calculations 

#Repeatability for ALL individuals from all years, for is_one binomial 
rpt_all <- rptBinary( #Use Binary version of rpt
  is_one ~ age10 + age_last_seen + total_resights +
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

##Individual level consistency 

#Build a dataframe pre-merge for modelling to calculate individual-level consistency
season_level <- Proportion_MPA %>% # has animalID, season, count_1_pup, total_resights
  mutate(proportion = count_1_pup / total_resights) %>% #Includes ALL data
  left_join(metadata %>% select(animalID, season, AgeYears),
            by = c("animalID","season")) %>%
  distinct() %>%
  mutate(animalID_fct = factor(animalID),
    season_fct   = factor(season))

#helper for weighted variance calculations
wvar <- function(x, w) {
  w <- w / sum(w) # (1) Normalize weights so they sum to 1
  mu <- sum(w * x) # (2) Compute the weighted mean of x
  sum(w * (x - mu)^2) # (3) Compute weighted average squared deviation from mean
}

#Calculate cosnsitency score for MPA for each individual
# Steps: Summarize per individual
#  - n_seasons: number of seasons observed
#  - n_obs_total: total number of resights
#  - mean_assoc: weighted mean proportion of association (0–1)
#  - sd_assoc: weighted SD across seasons, indicating behavioral variability
#  - frac_zero / frac_one: fraction of seasons where proportion = 0 or 1
#  - consistency: rescaled measure (1 − SD/0.5), bounded between 0–1
#                 where 1 = perfectly consistent, 0 = highly variable
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
  filter(n_seasons >= 3) #require at least two repeated breeding seasons per individual

#Plot consistency across entire dataset
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

#Classify individuals into behavioral strategies 
# Strategy classification logic:
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

#Plot these strategies with total observations (not season!) in consistency / MPA space 
ggplot(individual_consistency,
       aes(x = mean_assoc, y = consistency, color = strategy)) +
  geom_point(aes(size = n_obs_total, alpha = 0.8)) +
  geom_vline(xintercept = c(0.2, 0.8), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey60") +
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

##Additional consistency figures

##make a data table to create a star label for animals with consistent prop == 1
star_data <- individual_consistency %>%
  group_by(animalID_fct) %>%
  summarize(all_one = all(sd_assoc == 0),
            max_y = max(proportion)) %>%
  filter(all_one) %>%
  mutate(y_position = max_y + 0.05)

##boxplot with 2 color gradient for proportion standard deviation 
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

##one row per animal with its consistency
animal_order <- individual_consistency %>%
  distinct(animalID_fct, consistency) %>%
  arrange(desc(consistency)) %>%
  pull(animalID_fct)

##order animalID to consistency score
individual_consistency <- individual_consistency %>%
  mutate(animalID_fct = factor(animalID_fct, levels = animal_order))

##plot
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

################# Abiotic binomial and gamma models (1996-2023) ##################

##binomial model with extreme events per-year and group (SP or NP)
mod_abiotic_bin <- glmmTMB(is_one ~ n_extreme_both + group + total_resights + (1 | season_fct) + (1 | animalID_fct),
                           data = abiotic_variables,
                           family = binomial(link = "logit"))
summary(mod_abiotic_bin)
exp(fixef(mod_abiotic_bin)$cond) 
DHARMa::simulateResiduals(mod_abiotic_bin, plot = TRUE)
DHARMa::testDispersion(DHARMa::simulateResiduals(mod_abiotic_bin))
check_collinearity(mod_abiotic_bin)

##gamma model with extreme events per-year and group (SP or NP)
mod_abiotic_gam <- glmmTMB(flipped_prop ~ n_extreme_both + group + total_resights + (1 | season_fct) + (1 | animalID_fct),
                           data = abiotic_variables_sub,
                           family = Gamma(link = "log"))
summary(mod_abiotic_gam)
exp(fixef(mod_abiotic_gam)$cond) 
DHARMa::simulateResiduals(mod_abiotic_gam, plot = TRUE)
DHARMa::testDispersion(DHARMa::simulateResiduals(mod_abiotic_gam))
check_collinearity(mod_abiotic_gam)

