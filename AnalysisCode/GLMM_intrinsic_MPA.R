
######################## Piecewise binomial 1996 - 2025 ##################################

##Setting senescence threshold at 11 (Allison paper!)
age_senesce <- 11

##Setting threshold in data
intrinsic_variables <- intrinsic_variables %>%
mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
       levels = c("Young", "Old"))) %>%
mutate(age10 = (AgeYears - age_senesce) / 10) #scaled numeric version of Age centered at senescence threshold

##trying multiple versions of random effects to compare AIC

##Piecewise model with only season- lowest AIC
mod_piecewise_season <- glmmTMB(is_one ~ age10 : age_cat + age_last_seen + total_resights + (1 | season_fct),
                                     family = binomial(link = "logit"),
                                     data = intrinsic_variables)
summary(mod_piecewise_season)
exp(fixef(mod_piecewise_season)$cond) 
DHARMa::simulateResiduals(mod_piecewise_season, plot = TRUE)
#residuals look good here
check_collinearity(mod_piecewise_season)

##Piecewise model with year born and season
mod_piecewise_yrborn_season <- glmmTMB(is_one ~ age10 : age_cat + age_last_seen + total_resights + (1 | year_born_fct) + (1 | season_fct), 
                                   family = binomial(link = "logit"), 
                                   data = intrinsic_variables)
summary(mod_piecewise_yrborn_season)
exp(fixef(mod_piecewise_yrborn_season)$cond) 
DHARMa::simulateResiduals(mod_piecewise_yrborn_season, plot = TRUE)
#residuals look a bit worse than above model
check_collinearity(mod_piecewise_yrborn_season)

##Piecewise model with animalID, season, year_born
mod_piecewise_all <- glmmTMB(is_one ~ age10 : age_cat + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct) + (1 | year_born_fct),
                             family = binomial(link = "logit"),
                             data = intrinsic_variables)
summary(mod_piecewise_all)
exp(fixef(mod_piecewise_all)$cond) 
DHARMa::simulateResiduals(mod_piecewise_all, plot = TRUE)
#residuals look good here too
check_collinearity(mod_piecewise_all)

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
pred_season <- ggpredict(
  mod_piecewise_season, 
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
  season_fct = levels(intrinsic_variables$season_fct)) %>%
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old")),
         age10 = (AgeYears - age_senesce) / 10,
         total_resights = mean(intrinsic_variables$total_resights),
         age_last_seen = mean(intrinsic_variables$age_last_seen))

seal_years <- intrinsic_variables %>%
  count(season_fct)

##Make predictions
pred_grid <- pred_grid %>%
  mutate(pred_logit = predict(mod_piecewise_only_season, newdata = pred_grid, re.form = NULL),
         se_logit = predict(mod_piecewise_only_season, newdata = pred_grid, re.form = NA, se.fit = TRUE)$se.fit,
         pred_fixed = predict(mod_piecewise_only_season, newdata = pred_grid, re.form = NA)) %>%
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

##Lowest AIC and best residuals
mod_piecewise_gamma_all <- glmmTMB(flipped_prop ~ age10 : age_cat + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct) + (1 | year_born_fct),
                               family = Gamma(link = "log"),
                               data = intrinsic_variables_sub)
summary(mod_piecewise_gamma_all)
exp(fixef(mod_piecewise_gamma_all)$cond) 
DHARMa::simulateResiduals(mod_piecewise_gamma_all, plot = TRUE)
#residuals look not the best but fine
check_collinearity(mod_piecewise_gamma_all)

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
  mutate(age10 = (AgeYears - age_senesce) / 10) #scaled numeric version of Age centered at senescence threshold

##Piecewise model with season and animalID pre- and post- senescent threshold for 2016-2025- lowest AIC + good residuals
mod_piecewise_2016_2025 <- glmmTMB(is_one ~ age10 : age_cat + age_last_seen + total_resights + (1 | season_fct) + (1 | animalID_fct),
                                     family = binomial(link = "logit"),
                                     data = intrinsic_2016_2025)
summary(mod_piecewise_2016_2025)
exp(fixef(mod_piecewise_2016_2025)$cond)
DHARMa::simulateResiduals(mod_piecewise_2016_2025, plot = TRUE)
#residuals look good
check_collinearity(mod_piecewise_2016_2025)

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

##Subset the data for years 2016-2025
intrinsic_variables_sub_2016_2025 <- intrinsic_variables_sub %>%
  filter(season >= 2016)

##Compare models with combinations of random effects for AIC

##Model with animalID and season
mod_piecewise_gamma_2016_2025_both <- glmmTMB(flipped_prop ~ age10 : age_cat + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
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
