
source("./DataCurationCode/Data_processing_MPA.R") ##do the source code from Data Processing

######################## Piecewise binomial 1996 - 2025 ##################################

##Setting senescence threshold at 11 (Allison paper!)
age_senesce <- 11

##Setting threshold in data
intrinsic_variables <- intrinsic_variables %>%
mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
       levels = c("Young", "Old"))) %>%
mutate(age10 = (AgeYears - age_senesce) / 10) #scaled numeric version of Age centered at senescence threshold

#Make year_born continuous so we can scale and add as a fixed effect 
intrinsic_variables <- intrinsic_variables %>%
  mutate(year_born_cont = as.numeric(year_born_fct))

##trying multiple versions of random effects to compare AIC
#NOTE FROM MADDIE: I would stick with the third model (year_born as a numeric fixed effect), which seems to allow for year born & animalID to be controlled for 

##Piecewise model with only season- lowest AIC
mod_piecewise_season <- glmmTMB(is_one ~ age10 : age_cat + age_last_seen + total_resights + scale(year_born_cont) + (1 | season_fct),
                                     family = binomial(link = "logit"),
                                     data = intrinsic_variables)
summary(mod_piecewise_season)
exp(fixef(mod_piecewise_season)$cond) 
DHARMa::simulateResiduals(mod_piecewise_season, plot = TRUE)
#residuals look good here
check_collinearity(mod_piecewise_season)

##Piecewise model with year born and season
#mod_piecewise_yrborn_season <- glmmTMB(is_one ~ age10 : age_cat + age_last_seen + total_resights + (1 | year_born_fct) + (1 | season_fct), 
#                                   family = binomial(link = "logit"), 
#                                   data = intrinsic_variables)
#summary(mod_piecewise_yrborn_season)
#exp(fixef(mod_piecewise_yrborn_season)$cond) 
#DHARMa::simulateResiduals(mod_piecewise_yrborn_season, plot = TRUE)
#residuals look a bit worse than above model
#check_collinearity(mod_piecewise_yrborn_season)

##Piecewise model with animalID, season, year_born
mod_piecewise_all <- glmmTMB(is_one ~ age10 : age_cat + age_last_seen + total_resights + scale(year_born_cont) + (1 | animalID_fct) + (1 | season_fct),
                             family = binomial(link = "logit"),
                             data = intrinsic_variables)
summary(mod_piecewise_all)
exp(fixef(mod_piecewise_all)$cond) 
DHARMa::simulateResiduals(mod_piecewise_all, plot = TRUE) #Slight deviation in residuals (check dispersion)
DHARMa::testDispersion(DHARMa::simulateResiduals(mod_piecewise_all)) #No effect of over/under dispersion, which means we can ignore slight deviance in model quantiles!
check_collinearity(mod_piecewise_all) # <-- I would keep this one and comment out/remove the others!

##Visualizing and plotting the effects! -- Note from Maddie: I adjusted this code to be the same _all model as I suggested to keep above and it seems to have the same result 

##A table for raw observed proportion to overlay on plot
observed_data <- intrinsic_variables %>%
  group_by(AgeYears) %>%
  summarize(n = n(),
            n_one = sum(is_one),
            perc_one = n_one / n,
            lwr = binom.test(n_one, n)$conf.int[1],
            upr = binom.test(n_one, n)$conf.int[2])

##Create ANNUAL (i.e. rand effect of season only) predictions
pred_season <- ggpredict(
  mod_piecewise_all, 
  terms = c("age10 [all]", "age_cat", "season_fct"),
  type = "random") %>% 
  as_tibble() %>% 
  mutate(AgeYears = x * 10 + age_senesce,
         age_cat = factor(group, levels = c("Young", "Old")),
         season_fct = factor(facet)) %>%
  filter((age_cat == "Young" & AgeYears < age_senesce) | (age_cat == "Old" & AgeYears >= age_senesce))

##Need a prediction grid to visualize the Age effect
pred_grid <- expand.grid(
  AgeYears = 3:22,
  season_fct = levels(intrinsic_variables$season_fct),
  animalID_fct = NA  # placeholder since we won't predict per individual
) %>%
  mutate(
    age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                     levels = c("Young", "Old")),
    age10 = (AgeYears - age_senesce) / 10,
    total_resights = mean(intrinsic_variables$total_resights, na.rm = TRUE),
    age_last_seen = mean(intrinsic_variables$age_last_seen, na.rm = TRUE),
    year_born_cont = mean(intrinsic_variables$year_born_cont, na.rm = TRUE)
  )

seal_years <- intrinsic_variables %>%
  count(season_fct)

##Make predictions
pred_grid <- pred_grid %>%
  mutate(pred_logit = predict(mod_piecewise_all, newdata = pred_grid, re.form = NULL),
         se_logit = predict(mod_piecewise_all, newdata = pred_grid, re.form = NA, se.fit = TRUE)$se.fit,
         pred_fixed = predict(mod_piecewise_all, newdata = pred_grid, re.form = NA)) %>%
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
  labs(x = "Age (Years)", y = "Probability of 1", color = "Age class", fill = "Age class") +
  theme_few()


############################ Piecewise gamma model 1996-2025 ####################################

##Trying multiple versions of the same model to compare AIC

##Added this to make intrinsic_variables_sub compatible for everything! 

## FROM MADDIE: Make sure to keep only the scale(year_born) as fixed effect model here too (adjust these models) - delete/comment out the rest

intrinsic_variables_sub <- intrinsic_variables_sub %>%
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old"))) %>%
  mutate(age10 = (AgeYears - age_senesce) / 10) %>%
  left_join(harem_assignment, by = c("animalID", "season")) ### we need avg_density here merged w/ intrinsic data to add avg_density here

#Make year_born continuous so we can scale and add as a fixed effect 
intrinsic_variables_sub <- intrinsic_variables_sub %>%
  mutate(year_born_cont = as.numeric(year_born_fct))

##Lowest AIC and best residuals
mod_piecewise_gamma_all <- glmmTMB(flipped_prop ~ age10 : age_cat + age_last_seen + scale(year_born_cont) + (1 | animalID_fct) + (1 | season_fct),
                               family = Gamma(link = "log"),
                               data = intrinsic_variables_sub) ##Aditi : See comment above, add avg_density
summary(mod_piecewise_gamma_all)
exp(fixef(mod_piecewise_gamma_all)$cond) 
DHARMa::simulateResiduals(mod_piecewise_gamma_all, plot = TRUE)
#residuals look not the best but fine
check_collinearity(mod_piecewise_gamma_all)


intrinsic_variables_sub <- intrinsic_variables_sub %>%
  mutate(
    flipped_prop_trans = car::logit(flipped_prop, adjust = 0.001)
  )

## I would recommend this one! It can handle 0s which gamma is bad at, but we can keep the gamma as a sensitivty test 
mod_gauss <- glmmTMB(
  flipped_prop_trans ~ age10:age_cat + age_last_seen + scale(year_born_cont) +
    (1 | animalID_fct) + (1 | season_fct),
  family = gaussian(),
  data = intrinsic_variables_sub
)
summary(mod_gauss)
exp(fixef(mod_gauss)$cond) 
DHARMa::simulateResiduals(mod_gauss, plot = TRUE)
check_collinearity(mod_gauss)

###I think either stick with the gaussian one above here or the gamma, and delete/comment out the ones below 

mod_piecewise_gamma_two <- glmmTMB(flipped_prop ~ age10 : age_cat + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                                   family = Gamma(link = "log"),
                                   data = intrinsic_variables_sub)
summary(mod_piecewise_gamma_two)
exp(fixef(mod_piecewise_gamma_two)$cond) 
DHARMa::simulateResiduals(mod_piecewise_gamma_two, plot = TRUE)
#residuals look a bit worse
check_collinearity(mod_piecewise_gamma_two)

mod_piecewise_gamma_one <- glmmTMB(flipped_prop ~ age10 : age_cat + age_last_seen + total_resights + (1 | animalID_fct),
                                   family = Gamma(link = "log"),
                                   data = intrinsic_variables_sub)
summary(mod_piecewise_gamma_one)
exp(fixef(mod_piecewise_gamma_one)$cond) 
DHARMa::simulateResiduals(mod_piecewise_gamma_one, plot = TRUE)
#residuals look a bit worse
check_collinearity(mod_piecewise_gamma_one)

##################################### Piecewise binomial 2016 - 2025 ########################################

##Subset the data for only seasons 2016 - 2025
intrinsic_2016_2025 <- intrinsic_variables %>%
  filter(season >= 2016)

##Drop levels from the model with all years!
intrinsic_2016_2025 <- intrinsic_2016_2025 %>% droplevels()

##Add harem density data for 2016-2025
intrinsic_2016_2025 <- intrinsic_2016_2025 %>%
  left_join(harem_assignment, by = c("animalID", "season")) %>%
  left_join(area_density, by = c("area", "season"))

##Histogram of filtered data's proportions
ggplot(data = intrinsic_2016_2025, aes(x = proportion)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "darkblue") +
  labs(title = "Histogram of Proportion MPA",
       x = "Proportion MPA", y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_few()

##Create the age variables for the senescence model age_cat and age10
intrinsic_2016_2025 <- intrinsic_2016_2025 %>%
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old"))) %>%
  mutate(age10 = (AgeYears - age_senesce) / 10)  %>%   #scaled numeric version of Age centered at senescence threshold 
  mutate(year_born_cont = as.numeric(year_born_fct))
  
##Piecewise model with season and animalID pre- and post- senescent threshold for 2016-2025- lowest AIC + good residuals
mod_piecewise_2016_2025 <- glmmTMB(is_one ~ age10 : age_cat + avg_density + age_last_seen + total_resights + scale(year_born_cont) + (1 | season_fct) + (1 | animalID_fct),
                                     family = binomial(link = "logit"),
                                     data = intrinsic_2016_2025)
summary(mod_piecewise_2016_2025)
exp(fixef(mod_piecewise_2016_2025)$cond)
DHARMa::simulateResiduals(mod_piecewise_2016_2025, plot = TRUE)
#residuals look good
check_collinearity(mod_piecewise_2016_2025)

##A table for raw observed proportion to overlay on plot
#FROM MADDIE: Please adjust prediction graphs w/ new fixed effect 
observed_data_2016_2025 <- intrinsic_2016_2025 %>%
  group_by(AgeYears) %>%
  summarize(n = n(),
            n_one = sum(is_one),
            perc_one = n_one / n,
            lwr = binom.test(n_one, n)$conf.int[1],
            upr = binom.test(n_one, n)$conf.int[2])

##Create ANNUAL (i.e. rand effect of season only) predictions
pred_season_2016_2025 <- ggpredict(
  mod_piecewise_2016_2025, 
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
  animalID_fct = levels(intrinsic_2016_2025$animalID_fct)) %>%
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old")),
         age10 = (AgeYears - age_senesce) / 10,
         total_resights = mean(intrinsic_2016_2025$total_resights),
         age_last_seen = mean(intrinsic_2016_2025$age_last_seen))

seal_years_2016_2025 <- intrinsic_2016_2025 %>%
  count(season_fct)

##Make predictions
pred_grid_2016_2025 <- pred_grid_2016_2025 %>%
  mutate(pred_logit = predict(mod_piecewise_2016_2025, newdata = pred_grid_2016_2025, re.form = NULL),
         se_logit = predict(mod_piecewise_2016_2025, newdata = pred_grid_2016_2025, re.form = NA, se.fit = TRUE)$se.fit,
         pred_fixed = predict(mod_piecewise_2016_2025, newdata = pred_grid_2016_2025, re.form = NA)) %>%
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
  labs(x = "Age (Years)", y = "Probability of 1", color = "Age class", fill = "Age class") +
  theme_few()

############################# Piecewise gamma model 2016-2025 #####################################

#Maddie's notes from here: 
#Check out gamma vs. gaussian here (I would once again recommend the transformed gaussian model if histogram looks good)
#Add year_born fixed effect here as done above!

##Subset the data for years 2016-2025
intrinsic_variables_sub_2016_2025 <- intrinsic_variables_sub %>%
  filter(season >= 2016) %>%
  left_join(harem_assignment, by = c("animalID", "season")) %>%
  left_join(area_density, by = c("area", "season"))

##Compare models with combinations of random effects for AIC

##Model with animalID and season
mod_piecewise_gamma_2016_2025_both <- glmmTMB(flipped_prop ~ age10 : age_cat + avg_density + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                                         family = Gamma(link = "log"),
                                         data = intrinsic_variables_sub_2016_2025)
summary(mod_piecewise_gamma_2016_2025_both)
exp(fixef(mod_piecewise_gamma_2016_2025_both)$cond) 
DHARMa::simulateResiduals(mod_piecewise_gamma_2016_2025_both, plot = TRUE)
#residuals look ok
check_collinearity(mod_piecewise_gamma_2016_2025_both)

##Model with only animalID, lower AIC
mod_piecewise_gamma_2016_2025_ID_only <- glmmTMB(flipped_prop ~ age10 : age_cat + age_last_seen + total_resights + (1 | animalID_fct),
                                              family = Gamma(link = "log"),
                                              data = intrinsic_variables_sub_2016_2025)
summary(mod_piecewise_gamma_2016_2025_ID_only)
exp(fixef(mod_piecewise_gamma_2016_2025_ID_only)$cond) 
DHARMa::simulateResiduals(mod_piecewise_gamma_2016_2025_ID_only, plot = TRUE)
#residuals look ok
check_collinearity(mod_piecewise_gamma_2016_2025_ID_only)




######### MADDIE PLAYING WITH CONSISTENCY STUFF ############
#Make sure package is loaded in (can add to the top in future)
library(rptR)

#Subset by age to avoid need to account for age/specific piecewise split 
#For all data
young_data <- intrinsic_variables %>% filter(age_cat == "Young") 
old_data   <- intrinsic_variables %>% filter(age_cat == "Old")

#For 0 - 1 (gamma model) data
young_data_sub <- intrinsic_variables_sub %>% filter(age_cat == "Young")
old_data_sub   <- intrinsic_variables_sub %>% filter(age_cat == "Old")

##Population-level repeatability calculations 

# Repeatability for ALL individuals from all years, for is_one binomial 
rpt_all <- rptBinary( #Use Binary version of rpt
  is_one ~ age10 + age_last_seen + total_resights + scale(year_born_cont) + 
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = intrinsic_variables,
  link     = "logit",
  nboot    = 1000,
  npermut  = 0
); summary(rpt_all) #summarize results


# Repeatability for young individuals from all years, for is_one binomial 
rpt_young <- rptBinary( #Use Binary version of rpt
  is_one ~ age10 + age_last_seen + total_resights + scale(year_born_cont) + 
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = young_data,
  link     = "logit",
  nboot    = 1000,
  npermut  = 0
); summary(rpt_young) #summarize results

# Repeatability for old individuals from all years, for is_one binomial 
rpt_old <- rptBinary( #Use Binary version of rpt
  is_one ~ age10 + age_last_seen + total_resights + scale(year_born_cont) + 
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = old_data,
  link     = "logit",
  nboot    = 1000,
  npermut  = 0
); summary(rpt_old) #summarize results

#Repeatability for ALL individuals, all years, for between 0 and 1 
rpt_sub_all <- rpt( #Use regular rpt that works for normally distributed data
  flipped_prop_trans ~ age10 + age_last_seen + total_resights + scale(year_born_cont) +
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = intrinsic_variables_sub,
  datatype = "Gaussian",
  nboot    = 1000,
  npermut  = 0
); summary(rpt_sub_all) #summarize results

#Repeatability for young individuals, all years, for between 0 and 1 
rpt_young_sub <- rpt( #Use regular rpt that works for normally distributed data
  flipped_prop_trans ~ age10 + age_last_seen + total_resights + scale(year_born_cont) +
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = young_data_sub,
  datatype = "Gaussian",
  nboot    = 1000,
  npermut  = 0
)

summary(rpt_young_sub) #summarize results

#Repeatability for old individuals, all years, for between 0 and 1 
rpt_old_sub <- rpt( #Use regular rpt that works for normally distributed data
  flipped_prop_trans ~ age10 + age_last_seen + total_resights + scale(year_born_cont) +
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = old_data_sub,
  datatype = "Gaussian",
  nboot    = 1000,
  npermut  = 0
); summary(rpt_old_sub) #summarize results


#Individual level consistency 

#Build a dataframe pre-merge for modelling to calculate individual-level consistency
season_level <- Proportion_MPA %>%                 # has animalID, season, count_1_pup, total_resights
  mutate(proportion = count_1_pup / total_resights) %>% #Includes ALL data
  left_join(metadata %>% select(animalID, season, AgeYears),
            by = c("animalID","season")) %>%
  mutate(
    animalID_fct = factor(animalID),
    season_fct   = factor(season)
  )

# helper for weighted variance calculations
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
  summarise(
    n_seasons   = n(),
    n_obs_total = sum(total_resights),
    mean_assoc  = weighted.mean(proportion, w = total_resights, na.rm = TRUE),     # preserves direction
    sd_assoc    = sqrt(wvar(proportion, w = total_resights)),                       # variability across seasons
    frac_zero   = mean(proportion == 0),
    frac_one    = mean(proportion == 1),
    .groups = "drop"
  ) %>%
  mutate(
  # Consistency: standardized so that SD = 0 (no variation) → 1, and SD = 0.5 (max) → 0
    consistency = pmax(pmin(1 - (sd_assoc / 0.5), 1), 0)   # 0–1, higher = more consistent
  ) %>%
  filter(n_seasons >= 2) # require at least two repeated breeding seasons per individual

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
# - Maternal specialist: consistently associated (mean > 0.8, consistency > 0.8)
# - Non-maternal specialist: consistently unassociated (mean < 0.2, consistency > 0.8)
# - Plastic: fluctuates in association (consistency < 0.8)
# - Intermediate specialist: edge cases with moderate consistency near thresholds

individual_consistency <- individual_consistency %>%
  mutate(strategy = case_when(
    mean_assoc <= 0.2 & consistency >= 0.8 ~ "Non-maternal specialist",
    mean_assoc >= 0.8 & consistency >= 0.8 ~ "Maternal specialist",
    consistency < 0.8 ~ "Plastic",
    TRUE ~ "Intermediate specialist"
  ))

#Plot these strategies with total observations (not season!) in consistency / MPA space 
ggplot(individual_consistency,
       aes(x = mean_assoc, y = consistency, color = strategy)) +
  geom_point(aes(size = n_obs_total), alpha = 0.8) +
  geom_vline(xintercept = c(0.2, 0.8), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c(
    "Maternal specialist" = "#1f78b4",
    "Non-maternal specialist" = "#e31a1c",
    "Plastic" = "#ff7f00",
    "Intermediate specialist" = "#66C2A5"
  )) +
  scale_size_continuous(name = "Total observations", range = c(2,10)) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x = "Mean maternal association (0–1)",
    y = "Consistency (0–1)",
    color = "Behavioral strategy",
    title = "Individual consistency and maternal behavior strategy"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.title = element_text(face = "bold")
  )
