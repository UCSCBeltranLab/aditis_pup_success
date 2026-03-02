##run the source code from Data Processing
source("./DataCurationCode/Data_processing_MPA.R")

### test collinearity between model variables ###
cor.test(intrinsic_variables$n_extreme_both, intrinsic_variables$season, use = "complete.obs")

########################## Model with all predictors (2016-2023) ###########################

# 1) Subset the data for only seasons 2016 - 2025
intrinsic_2016_2023 <- intrinsic_variables %>%
  filter(season %in% 2016 : 2023) %>%
  droplevels() ##Drop levels from the model with all years!

### Piecewise model for is_one 2016-2023 ###
mod_binom_2016_2023 <- glmmTMB(proportion ~ age_cat : age10 + n_extreme_both + avg_density + (1 | animalID_fct) + (1 | season_fct),
                               weights = total_resights,
                               family = binomial(link = "logit"),
                               data = intrinsic_2016_2023); summary(mod_binom_2016_2023) #model summary
ranef(mod_binom_2016_2023)
exp(fixef(mod_binom_2016_2023)$cond) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_binom_2016_2023, plot = TRUE) #plot all residuals 
testDispersion(simulateResiduals(mod_binom_2016_2023)) #check for overdispersion
check_collinearity(mod_binom_2016_2023) #check predictor VIFs

# 1) A table for raw observed proportion to overlay on plot
observed_data_2016_2023 <- intrinsic_2016_2023 %>%
  group_by(AgeYears) %>%
  summarize(n = n(),
            n_success = sum(count_1_pup),
            n_trials  = sum(total_resights),
            avg_prop  = mean(count_1_pup / total_resights),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(lwr = binom.test(n_success, n_trials)$conf.int[1],
         upr = binom.test(n_success, n_trials)$conf.int[2]) %>%
  ungroup()

# 2) Create ANNUAL (rand effect of season only) predictions
pred_season_2016_2023 <- ggpredict(mod_binom_2016_2023, 
                                   terms = c("age10 [all]", "age_cat", "season_fct"), 
                                   type = "random") %>%
  as_tibble() %>%
  mutate(AgeYears = x * 10 + age_senesce,
         age_cat = factor(group, levels = c("Young", "Old")),
         season_fct = factor(facet)) %>%
  filter((age_cat == "Young" & AgeYears < age_senesce) | (age_cat == "Old" & AgeYears >= age_senesce))

# 3) Create a prediction grid to visualize the Age effect
pred_grid_2016_2023 <- expand.grid(AgeYears = 3:21,
                                   season_fct = levels(intrinsic_2016_2023$season_fct),
                                   animalID_fct = NA) %>% #placeholder since we won't predict per individual
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old")),
         age10 = (AgeYears - age_senesce) / 10,
         total_resights = mean(intrinsic_2016_2023$total_resights, na.rm = TRUE),
         avg_density = mean(intrinsic_2016_2023$avg_density, na.rm = TRUE),
         n_extreme_both = mean(intrinsic_2016_2023$n_extreme_both, na.rm =TRUE)) #use means for predictors we want to hold constant

##identify number of observations per year
seal_years_2016_2023 <- intrinsic_2016_2023 %>%
  count(season_fct)

# 4) Make predictions for age effect
pred_grid_2016_2023 <- pred_grid_2016_2023 %>%
  mutate(pred_logit = predict(mod_binom_2016_2023, newdata = pred_grid_2016_2023, re.form = NULL), #linear predictor (logit scale) using fixed effects only
         se_logit = predict(mod_binom_2016_2023, newdata = pred_grid_2016_2023, re.form = NA, se.fit = TRUE)$se.fit) %>% #standard error of prediction on logit scale
  mutate(conf_lo = pred_logit - 1.96 * se_logit, 
         conf_hi = pred_logit + 1.96 * se_logit) %>% #Compute 95% confidence intervals (still on logit scale)
  left_join(seal_years_2016_2023, by = "season_fct") %>% #Add sample size weights for each season
  group_by(AgeYears, age_cat) %>%
  #Compute weighted summary statistics on probability scale
  summarise(predicted = weighted.mean(plogis(pred_logit), n), #avg predicted probability
            conf_lo = weighted.mean(plogis(conf_lo), n), #CI lower
            conf_hi = weighted.mean(plogis(conf_hi), n), #CI upper
            .groups = "drop")

# 5) Plot predictions with senescent threshold
ggplot(pred_grid_2016_2023, aes(x = AgeYears, y = predicted, color = age_cat, fill = age_cat)) +
  
  #Predictions for individual seasons
  geom_line(aes(group = interaction(age_cat, season_fct)),
            pred_season_2016_2023,
            alpha = 0.3, color = "#7EAAC1") +
  
  #Main ribbon for model w/ CI
  geom_ribbon(data = pred_grid_2016_2023,
              aes(x = AgeYears, ymin = conf_lo, ymax = conf_hi, fill = age_cat),
              alpha = 0.3, color = NA, inherit.aes = FALSE) +
  
  #Thick colored prediction line
  geom_line(data = pred_grid_2016_2023,
            aes(x = AgeYears, y = predicted, color = age_cat),
            linewidth = 1.2, inherit.aes = FALSE) +
  
  #Raw observed points
  geom_pointrange(data = observed_data_2016_2023,
                  aes(x = AgeYears, y = avg_prop, ymin = lwr, ymax = upr),
                  inherit.aes = FALSE, color = "black") +
  
  #Overlay sample size for each age
  geom_text(data = observed_data_2016_2023,
            aes(x = AgeYears, y = 1.01, label = n),
            inherit.aes = FALSE, size = 3) +
  
  #Formatting, colors
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_y_continuous(n.breaks = 5) + 
  scale_x_continuous(n.breaks = 20) +
  scale_color_manual(values = c("Young" = "#92BAEE", "Old" = "#EB99D2")) +
  scale_fill_manual(values = c("Young" = "#92BAEE", "Old" = "#EB99D2")) +
  labs(x = "Age (Years)", 
       y = "Proportion", 
       color = "Age class", 
       fill = "Age class") +
  theme_few()

# 1) predicted effect of density on flipped proportion
pred_density_2016_2023 <- ggpredict(mod_binom_2016_2023, terms = "avg_density [all]")

# 2) plot flipped proportion and average density
ggplot(data = pred_density_2016_2023, aes(x = x, y = predicted)) +
  
  #Raw points
  geom_point(data = intrinsic_2016_2023,
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
  theme_minimal()

######################## Model with age only (1996 - 2025) ##################################

### Piecewise model for is_one - full dataset ###
mod_binom_1996_2025 <- glmmTMB(proportion ~ age_cat : age10 + (1 | animalID_fct) + (1 | season_fct),
                               weights = total_resights,
                               family = binomial(link= "logit"),
                               data = intrinsic_variables); summary(mod_binom_1996_2025) #model summary
ranef(mod_binom_1996_2025)
exp(fixef(mod_binom_1996_2025)$cond) #converts fixed-effect log-odds to odds ratios
resid <- simulateResiduals(mod_binom_1996_2025, plot = TRUE) #plot all residuals testDispersion(simulateResiduals(mod_binom_1996_2025)) #check for overdispersion
plotResiduals(resid, form = intrinsic_variables$age_cat)
check_collinearity(mod_binom_1996_2025) #check predictor VIFs

observed_data <- intrinsic_variables %>%
  group_by(AgeYears) %>%
  summarise(
    n = n(),
    n_success = sum(count_1_pup, na.rm=TRUE),
    n_trials  = sum(total_resights, na.rm=TRUE),
    avg_prop  = n_success / n_trials,
    lwr = binom.test(n_success, n_trials)$conf.int[1],
    upr = binom.test(n_success, n_trials)$conf.int[2],
    .groups="drop")

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
pred_grid <- expand.grid(AgeYears = 3:21,
                         season_fct = levels(intrinsic_variables$season_fct),
                         animalID_fct = NA) %>% #placeholder since we won't predict per individual
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old")),
         age10 = (AgeYears - age_senesce) / 10,
         total_resights = mean(intrinsic_variables$total_resights, na.rm = TRUE)) #use means for predictors we want to hold constant

seal_years <- intrinsic_variables %>%
  group_by(season_fct) %>%
  summarise(trials = sum(total_resights, na.rm=TRUE), .groups="drop")

pred_grid <- pred_grid %>%
  mutate(pred_logit = predict(mod_binom_1996_2025, newdata=pred_grid, re.form= NULL, type="link"),
         se_logit = predict(mod_binom_1996_2025, newdata=pred_grid, re.form=NA, se.fit=TRUE, type="link")$se.fit) %>%
  mutate(conf_lo = pred_logit - 1.96*se_logit,
         conf_hi = pred_logit + 1.96*se_logit) %>%
  left_join(seal_years, by="season_fct") %>%
  group_by(AgeYears, age_cat) %>%
  summarise(predicted = weighted.mean(plogis(pred_logit), trials),
            conf_lo   = weighted.mean(plogis(conf_lo),   trials),
            conf_hi   = weighted.mean(plogis(conf_hi),   trials),
            .groups="drop")

# 5) Plot predictions with senescent threshold
age_plot <- ggplot(pred_grid, aes(x = AgeYears, y = predicted, color = age_cat, fill = age_cat)) +
  
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
                  aes(x = AgeYears, y = avg_prop, ymin = lwr, ymax = upr),
                  inherit.aes = FALSE, color = "black") +
  
  geom_text(data = observed_data,
            aes(x = AgeYears, y = 1.01, label = n),
            inherit.aes = FALSE, size = 3) +
  
  #Formatting, colors
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_y_continuous(n.breaks = 5) + 
  scale_x_continuous(n.breaks = 20) +
  scale_color_manual(values = c("Young" = "#92BAEE", "Old" = "#EB99D2")) +
  scale_fill_manual(values = c("Young" = "#92BAEE", "Old" = "#EB99D2")) +
  labs(x = "Age (Years)", 
       y = "Proportion", 
       color = "Age class", 
       fill = "Age class") +
  theme_few(); age_plot

##################################### Model with age + density (2016 - 2025) ########################################

# 1) Subset the data for only seasons 2016 - 2025
intrinsic_2016_2025 <- intrinsic_variables %>%
  filter(season >= 2016) %>%
  droplevels() ##Drop levels from the model with all years!

### Piecewise model for is_one 2016-2025 ###
mod_binom_2016_2025 <- glmmTMB(proportion ~ age_cat : age10 + avg_density + (1 | animalID_fct) + (1 | season_fct),
                               weights = total_resights,
                               family = binomial(link = "logit"),
                               data = intrinsic_2016_2025); summary(mod_binom_2016_2025) #model summary
ranef(mod_binom_2016_2025)
exp(fixef(mod_binom_2016_2025)$cond) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_binom_2016_2025, plot = TRUE) #plot all residuals 
testDispersion(simulateResiduals(mod_binom_2016_2025)) #check for overdispersion
check_collinearity(mod_binom_2016_2025) #check predictor VIFs


# 1) A table for raw observed proportion to overlay on plot
observed_data_2016_2025 <- intrinsic_2016_2025 %>%
  group_by(AgeYears) %>%
  summarize(n = n(),
            n_success = sum(count_1_pup),
            n_trials  = sum(total_resights),
            avg_prop  = mean(count_1_pup / total_resights),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(lwr = binom.test(n_success, n_trials)$conf.int[1],
         upr = binom.test(n_success, n_trials)$conf.int[2]) %>%
  ungroup()

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
pred_grid_2016_2025 <- expand.grid(AgeYears = 3:21,
                                   season_fct = levels(intrinsic_2016_2025$season_fct),
                                   animalID_fct = NA) %>% #placeholder since we won't predict per individual
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old")),
         age10 = (AgeYears - age_senesce) / 10,
         total_resights = mean(intrinsic_2016_2025$total_resights, na.rm = TRUE),
         avg_density = mean(intrinsic_2016_2025$avg_density, na.rm = TRUE)) #use means for predictors we want to hold constant

##identify number of observations per year
seal_years <- intrinsic_2016_2025 %>%
  count(season_fct)

# 4) Make predictions for age effect
pred_grid_2016_2025 <- pred_grid_2016_2025 %>%
  mutate(pred_logit = predict(mod_binom_2016_2025, newdata = pred_grid_2016_2025, re.form = NULL), #linear predictor (logit scale) using fixed effects only
         se_logit = predict(mod_binom_2016_2025, newdata = pred_grid_2016_2025, re.form = NA, se.fit = TRUE)$se.fit) %>% #standard error of prediction on logit scale
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
                  aes(x = AgeYears, y = avg_prop, ymin = lwr, ymax = upr),
                  inherit.aes = FALSE, color = "black") +
  
  #Overlay sample size for each age
  geom_text(data = observed_data_2016_2025,
            aes(x = AgeYears, y = 1.01, label = n),
            inherit.aes = FALSE, size = 3) +
  
  #Formatting, colors
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_y_continuous(n.breaks = 5) + 
  scale_x_continuous(n.breaks = 20) +
  scale_color_manual(values = c("Young" = "#92BAEE", "Old" = "#EB99D2")) +
  scale_fill_manual(values = c("Young" = "#92BAEE", "Old" = "#EB99D2")) +
  labs(x = "Age (Years)", 
       y = "Proportion", 
       color = "Age class", 
       fill = "Age class") +
  theme_few()

# 1) predicted effect of density on flipped proportion
pred_density <- ggpredict(mod_binom_2016_2025, terms = "avg_density [all]")

# 2) plot flipped proportion and average density
ggplot(data = pred_density, aes(x = x, y = predicted)) +
  
  #Raw points
  geom_point(data = intrinsic_2016_2025,
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
  theme_minimal()

################# Model with age + extreme wave/tide events (1996-2023) ##################
# 1) Subset the data for only seasons 2016 - 2025
intrinsic_1996_2023 <- intrinsic_variables %>%
  filter(season %in% 1996 : 2023) %>%
  droplevels() ##Drop levels from the model with all years!

### Piecewise binomial model with per-year extreme wave/tide events ###
mod_binom_1996_2023 <- glmmTMB(proportion ~ age_cat : age10 + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                           weights = total_resights,
                           data = intrinsic_1996_2023,
                           family = binomial(link = "logit")); summary(mod_binom_1996_2023) #model summary
ranef(mod_binom_1996_2023)
exp(fixef(mod_binom_1996_2023)$cond) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_binom_1996_2023, plot = TRUE) #plot all residuals 
testDispersion(simulateResiduals(mod_binom_1996_2023)) #check for overdispersion
check_collinearity(mod_binom_1996_2023) #check predictor VIFs

# 1) A table for raw observed proportion to overlay on plot
observed_data_1996_2023 <- intrinsic_1996_2023 %>%
  group_by(AgeYears) %>%
  summarize(n = n(),
            n_success = sum(count_1_pup),
            n_trials  = sum(total_resights),
            avg_prop  = mean(count_1_pup / total_resights),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(lwr = binom.test(n_success, n_trials)$conf.int[1],
         upr = binom.test(n_success, n_trials)$conf.int[2]) %>%
  ungroup()

# 2) Create ANNUAL (rand effect of season only) predictions
pred_season_1996_2023 <- ggpredict(mod_binom_1996_2023, 
                                   terms = c("age10 [all]", "age_cat", "season_fct"), 
                                   type = "random") %>%
  as_tibble() %>%
  mutate(AgeYears = x * 10 + age_senesce,
         age_cat = factor(group, levels = c("Young", "Old")),
         season_fct = factor(facet)) %>%
  filter((age_cat == "Young" & AgeYears < age_senesce) | (age_cat == "Old" & AgeYears >= age_senesce))

# 3) Create a prediction grid to visualize the Age effect
pred_grid_1996_2023 <- expand.grid(AgeYears = 3:21,
                                   season_fct = levels(intrinsic_1996_2023$season_fct),
                                   animalID_fct = NA) %>% #placeholder since we won't predict per individual
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old")),
         age10 = (AgeYears - age_senesce) / 10,
         total_resights = mean(intrinsic_1996_2023$total_resights, na.rm = TRUE),
         n_extreme_both = mean(intrinsic_1996_2023$n_extreme_both, na.rm = TRUE)) #use means for predictors we want to hold constant

##identify number of observations per year
seal_years <- intrinsic_1996_2023 %>%
  count(season_fct)

# 4) Make predictions for age effect
pred_grid_1996_2023 <- pred_grid_1996_2023 %>%
  mutate(pred_logit = predict(mod_binom_1996_2023, newdata = pred_grid_1996_2023, re.form = NULL), #linear predictor (logit scale) using fixed effects only
         se_logit = predict(mod_binom_1996_2023, newdata = pred_grid_1996_2023, re.form = NA, se.fit = TRUE)$se.fit) %>% #standard error of prediction on logit scale
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
ggplot(pred_grid_1996_2023, aes(x = AgeYears, y = predicted, color = age_cat, fill = age_cat)) +
  
  #Predictions for individual seasons
  geom_line(aes(group = interaction(age_cat, season_fct)),
            pred_season_1996_2023,
            alpha = 0.3, color = "#7EAAC1") +
  
  #Main ribbon for model w/ CI
  geom_ribbon(data = pred_grid_1996_2023,
              aes(x = AgeYears, ymin = conf_lo, ymax = conf_hi, fill = age_cat),
              alpha = 0.3, color = NA, inherit.aes = FALSE) +
  
  #Thick colored prediction line
  geom_line(data = pred_grid_1996_2023,
            aes(x = AgeYears, y = predicted, color = age_cat),
            linewidth = 1.2, inherit.aes = FALSE) +
  
  #Raw observed points
  geom_pointrange(data = observed_data_1996_2023,
                  aes(x = AgeYears, y = avg_prop, ymin = lwr, ymax = upr),
                  inherit.aes = FALSE, color = "black") +
  
  #Overlay sample size for each age
  geom_text(data = observed_data_1996_2023,
            aes(x = AgeYears, y = 1.01, label = n),
            inherit.aes = FALSE, size = 3) +
  
  #Formatting, colors
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_y_continuous(n.breaks = 5) + 
  scale_x_continuous(n.breaks = 20) +
  scale_color_manual(values = c("Young" = "#92BAEE", "Old" = "#EB99D2")) +
  scale_fill_manual(values = c("Young" = "#92BAEE", "Old" = "#EB99D2")) +
  labs(x = "Age (Years)", 
       y = "Proportion", 
       color = "Age class", 
       fill = "Age class") +
  theme_few()

### Plot effects of group (NP or SP) on MOA ###

# 1) predicted effect of group on MOA
effect_group <- ggpredict(mod_binom_1996_2023, terms = "group_fct")

# 2) plot group vs. predicted MOA
ggplot(effect_group, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  geom_point(size = 3) +
  labs(x = "Group",
       y = "Mother-offspring association") +
  scale_y_continuous(n.breaks = 8, limits = c(0.9,1)) +
  theme_minimal()

# 1) predicted effect of n_extreme_both on MOA
weather_effect <- ggpredict(mod_binom_1996_2023, terms = "n_extreme_both")

# 2) plot weather vs. predicted MOA
ggplot(weather_effect, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.15, fill = "black") +
  geom_line(linewidth = 1.1, color = "black") +
  geom_point(data = intrinsic_1996_2023,
             aes(x = n_extreme_both,
                 y = proportion),
             inherit.aes = FALSE,
             alpha = 0.4, size = 2) +
  labs(x = "Number of extreme wave/tide events (per year)",
       y = "Mother–offspring association (proportion)",
       color = "Season") +
  theme_minimal()
