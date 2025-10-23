
source("./DataCurationCode/Data_processing_MPA.R") ##do the source code from Data Processing

######################## Piecewise binomial 1996 - 2025 ##################################

#NOTE FROM MADDIE: I would stick with the third model (year_born as a numeric fixed effect), which seems to allow for year born & animalID to be controlled for 

##Piecewise model with animalID, season, year_born
mod_binom_1996_2025 <- glmmTMB(is_one ~ age10 : age_cat + age_last_seen + total_resights + scale(year_born_num) + (1 | animalID_fct) + (1 | season_fct),
                             family = binomial(link = "logit"),
                             data = intrinsic_variables)
summary(mod_binom_1996_2025)
exp(fixef(mod_binom_1996_2025)$cond) 
DHARMa::simulateResiduals(mod_binom_1996_2025, plot = TRUE) #Slight deviation in residuals (check dispersion)
DHARMa::testDispersion(DHARMa::simulateResiduals(mod_binom_1996_2025)) #No effect of over/under dispersion, which means we can ignore slight deviance in model quantiles!
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
  labs(x = "Age (Years)", y = "Probability of 1", color = "Age class", fill = "Age class") +
  theme_few()

# Get model summary table
tab <- tidy(mod_binom_1996_2025, effects = "fixed") %>%
  # Drop effect/component columns (if any remain)
  select(term, estimate, std.error, statistic, p.value) %>%
  # Add significance markers
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.10  ~ "·",
      TRUE ~ ""
    ))

# Then feed it into a table display
sjPlot::tab_df(
  tab,
  title = "Binomial Results (1996-2025)",
  show.rownames = FALSE,
  digits = 3
)


############################ Piecewise gamma model 1996-2025 ####################################

##FROM MADDIE: Make sure to keep only the scale(year_born) as fixed effect model here too (adjust these models) - delete/comment out the rest

mod_gamma_1996_2025 <- glmmTMB(flipped_prop ~ age10 : age_cat + age_last_seen + total_resights + scale(year_born_num) + (1 | animalID_fct) + (1 | season_fct),
                                   family = Gamma(link = "log"),
                                   data = intrinsic_variables_sub) ##Aditi : See comment above, add avg_density
summary(mod_gamma_1996_2025)
exp(fixef(mod_gamma_1996_2025)$cond) 
DHARMa::simulateResiduals(mod_gamma_1996_2025, plot = TRUE)
DHARMa::testDispersion(DHARMa::simulateResiduals(mod_binom_1996_2025))
check_collinearity(mod_gamma_1996_2025)


##FROM MADDIE: I would recommend this one! It can handle 0s which gamma is bad at, but we can keep the gamma as a sensitivty test 
mod_gauss_1996_2025 <- glmmTMB(prop_trans ~ age10:age_cat + age_last_seen + total_resights + scale(year_born_num) + (1 | animalID_fct) + (1 | season_fct),
                     family = gaussian(),
                     data = intrinsic_variables_sub)
summary(mod_gauss_1996_2025)
exp(fixef(mod_gauss_1996_2025)$cond) 
DHARMa::simulateResiduals(mod_gauss_1996_2025, plot = TRUE)
check_collinearity(mod_gauss_1996_2025)

library(broom)
library(broom.mixed)
library(dplyr)

# Get model summary table
tab <- tidy(mod_gamma_1996_2025, effects = "fixed") %>%
  # Drop effect/component columns (if any remain)
  select(term, estimate, std.error, statistic, p.value) %>%
  # Add significance markers
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.10  ~ "·",
      TRUE ~ ""
    )
  )

# Then feed it into a table display
sjPlot::tab_df(
  tab,
  title = "Gamma Results (1996-2025)",
  show.rownames = FALSE,
  digits = 3
)

##A table for raw observed proportion to overlay on plot
observed_data_1996_2025_sub <- intrinsic_variables_sub %>%
  group_by(AgeYears) %>%
  summarize(
    n = n(),
    mean_prop = mean(flipped_prop, na.rm = TRUE),
    sd_prop = sd(flipped_prop, na.rm = TRUE),
    se = sd_prop / sqrt(n),
    lwr = mean_prop - 1.96 * se,
    upr = mean_prop + 1.96 * se
  )

##Create ANNUAL (i.e. rand effect of season only) predictions
pred_season_1996_2025_sub <- ggpredict(
  mod_gamma_1996_2025, 
  terms = c("age10 [all]", "age_cat", "season_fct"),
  type = "random"
) %>% 
  as_tibble() %>% 
  mutate(
    AgeYears = x * 10 + age_senesce,
    age_cat = factor(group, levels = c("Young", "Old")),
    season_fct = factor(facet)
  ) %>%
  filter(
    (age_cat == "Young" & AgeYears < age_senesce) | 
      (age_cat == "Old" & AgeYears >= age_senesce)
  )

##Need a prediction grid to visualize the Age effect
pred_grid_1996_2025_sub <- expand.grid(
  AgeYears = 3:22
) %>%
  mutate(
    age_cat = factor(
      ifelse(AgeYears < age_senesce, "Young", "Old"),
      levels = c("Young", "Old")
    ),
    age10 = (AgeYears - age_senesce) / 10,
    total_resights = mean(intrinsic_variables_sub$total_resights, na.rm = TRUE),
    age_last_seen  = mean(intrinsic_variables_sub$age_last_seen, na.rm = TRUE),
    year_born_num  = mean(intrinsic_variables_sub$year_born_num, na.rm = TRUE),
    animalID_fct = factor(
      levels(model.frame(mod_gamma_1996_2025)$animalID_fct)[1],
      levels = levels(model.frame(mod_gamma_1996_2025)$animalID_fct)),
    season_fct = factor(
      levels(model.frame(mod_gamma_1996_2025)$season_fct)[1],
      levels = levels(model.frame(mod_gamma_1996_2025)$season_fct))
  )

seal_years_1996_2025_sub <- intrinsic_variables_sub %>%
  count(season_fct)

##Make predictions
pred_out <- predict(
  mod_gamma_1996_2025,
  newdata = pred_grid_1996_2025_sub,
  type = "link",
  re.form = NA,
  se.fit = TRUE,
  allow.new.levels = TRUE
)

pred_grid_1996_2025_sub <- pred_grid_1996_2025_sub %>%
  mutate(pred_link = pred_out$fit,
    se_link = pred_out$se.fit,
    # use plogis if model uses a logit link, exp() if log link
    predicted = plogis(pred_link),
    conf_lo = plogis(pred_link - 1.96 * se_link),
    conf_hi = plogis(pred_link + 1.96 * se_link)
  ) %>%
  left_join(seal_years_1996_2025_sub, by = "season_fct") %>%
  group_by(AgeYears, age_cat) %>%
  summarise(
    predicted = weighted.mean(predicted, n),
    conf_lo = weighted.mean(conf_lo, n),
    conf_hi = weighted.mean(conf_hi, n),
    .groups = "drop"
  )

##Plot predictions with senescent threshold
ggplot(pred_grid_1996_2025_sub, aes(x = AgeYears, y = predicted, color = age_cat, fill = age_cat)) +
  
  #Predictions for individual seasons
  geom_line(aes(group = interaction(age_cat, season_fct)),
            pred_season_1996_2025_sub,
            alpha = 0.3, color = "#7EAAC1") +
  
  #Main ribbon for model w/ CI
  geom_ribbon(data = pred_grid_1996_2025_sub,
              aes(x = AgeYears, ymin = conf_lo, ymax = conf_hi, fill = age_cat),
              alpha = 0.3, color = NA, inherit.aes = FALSE) +
  
  #Thick colored prediction line
  geom_line(data = pred_grid_1996_2025_sub,
            aes(x = AgeYears, y = predicted, color = age_cat),
            linewidth = 1.2, inherit.aes = FALSE) +
  
  #Raw observed points
  geom_pointrange(data = observed_data_1996_2025_sub,
                  aes(x = AgeYears, y = mean_prop, ymin = lwr, ymax = upr),
                  inherit.aes = FALSE, color = "black") +
  
  #Overlay sample size for each age
  geom_text(data = observed_data_1996_2025_sub,
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
  
##Piecewise model for 2016-2025
mod_binom_2016_2025 <- glmmTMB(is_one ~ age10 : age_cat + avg_density + age_last_seen + total_resights + scale(year_born_num) + (1 | season_fct) + (1 | animalID_fct),
                                     family = binomial(link = "logit"),
                                     data = intrinsic_2016_2025)
summary(mod_binom_2016_2025)
exp(fixef(mod_binom_2016_2025)$cond)
DHARMa::simulateResiduals(mod_binom_2016_2025, plot = TRUE)
#residuals look good
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
         avg_density = mean(intrinsic_2016_2025$avg_density, na.rm = TRUE),
         year_born_num = mean(intrinsic_2016_2025$year_born_num, na.rm = TRUE))

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

# Get model summary table
tab <- tidy(mod_binom_2016_2025, effects = "fixed") %>%
  # Drop effect/component columns (if any remain)
  select(term, estimate, std.error, statistic, p.value) %>%
  # Add significance markers
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.10  ~ "·",
      TRUE ~ ""
    )
  )

# Then feed it into a table display
sjPlot::tab_df(
  tab,
  title = "Binomial Results (2016-2025)",
  show.rownames = FALSE,
  digits = 3
)

############################# Piecewise gamma model 2016-2025 #####################################

#Maddie's notes from here: 
#Check out gamma vs. gaussian here (I would once again recommend the transformed gaussian model if histogram looks good)
#Add year_born fixed effect here as done above!

##Subset the data for years 2016-2025
intrinsic_sub_2016_2025 <- intrinsic_variables_sub %>%
  filter(season >= 2016) %>%
  left_join(harem_assignment, by = c("animalID", "season")) %>%
  left_join(area_density, by = c("area", "season"))

##Compare models with combinations of random effects for AIC

##Model with animalID and season
mod_gamma_2016_2025 <- glmmTMB(flipped_prop ~ age_cat : age10 + avg_density + age_last_seen + scale(year_born_num) + total_resights + (1 | animalID_fct) + (1 | season_fct),
                                         family = Gamma(link = "log"),
                                         data = intrinsic_sub_2016_2025)
summary(mod_gamma_2016_2025)
exp(fixef(mod_gamma_2016_2025)$cond) 
DHARMa::simulateResiduals(mod_gamma_2016_2025, plot = TRUE)
#residuals look ok
check_collinearity(mod_gamma_2016_2025)

#Gaussian model 2016-2025
mod_gauss_2016_2025 <- glmmTMB(prop_trans ~ age_cat : age10 + avg_density + age_last_seen + scale(year_born_num) + total_resights + (1 | animalID_fct) + (1 | season_fct),
                               family = gaussian(),
                               data = intrinsic_sub_2016_2025)
summary(mod_gauss_2016_2025)
exp(fixef(mod_gauss_2016_2025)$cond) 
DHARMa::simulateResiduals(mod_gauss_2016_2025, plot = TRUE)
check_collinearity(mod_gauss_2016_2025)

# Get model summary table
tab <- tidy(mod_gamma_2016_2025, effects = "fixed") %>%
  # Drop effect/component columns (if any remain)
  select(term, estimate, std.error, statistic, p.value) %>%
  # Add significance markers
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.10  ~ "·",
      TRUE ~ ""
    )
  )

# Then feed it into a table display
sjPlot::tab_df(
  tab,
  title = "Gamma Results (2016-2025)",
  show.rownames = FALSE,
  digits = 3
)

##A table for raw observed proportion to overlay on plot
observed_data_2016_2025_sub <- intrinsic_sub_2016_2025 %>%
  group_by(AgeYears) %>%
  summarize(n = n(),
            mean_prop = mean(flipped_prop, na.rm = TRUE),
            sd_prop = sd(flipped_prop, na.rm = TRUE),
            se = sd_prop / sqrt(n),
            lwr = mean_prop - 1.96 * se,
            upr = mean_prop + 1.96 * se)

##Create ANNUAL (i.e. rand effect of season only) predictions
pred_season_2016_2025_sub <- ggpredict(
  mod_gamma_2016_2025, 
  terms = c("age10 [all]", "age_cat", "season_fct"),
  type = "random") %>% 
  as_tibble() %>% 
  mutate(AgeYears = x * 10 + age_senesce,
         age_cat = factor(group, levels = c("Young", "Old")),
         season_fct = factor(facet)) %>%
  filter((age_cat == "Young" & AgeYears < age_senesce) | (age_cat == "Old" & AgeYears >= age_senesce))

##Need a prediction grid to visualize the Age effect
pred_grid_2016_2025_sub <- expand.grid(
  AgeYears = 3:22) %>%
  mutate(
    age_cat = factor(
      ifelse(AgeYears < age_senesce, "Young", "Old"),
      levels = c("Young", "Old")
    ),
    age10 = (AgeYears - age_senesce) / 10,
    total_resights = mean(intrinsic_sub_2016_2025$total_resights, na.rm = TRUE),
    age_last_seen  = mean(intrinsic_sub_2016_2025$age_last_seen, na.rm = TRUE),
    avg_density    = mean(intrinsic_sub_2016_2025$avg_density, na.rm = TRUE),
    year_born_num = mean(intrinsic_sub_2016_2025$year_born_num, na.rm = TRUE),
    animalID_fct = factor(levels(model.frame(mod_gamma_2016_2025)$animalID_fct)[1],
                          levels = levels(model.frame(mod_gamma_2016_2025)$animalID_fct)),
    season_fct   = factor(levels(model.frame(mod_gamma_2016_2025)$season_fct)[1],
                          levels = levels(model.frame(mod_gamma_2016_2025)$season_fct))
  )

seal_years_2016_2025_sub <- intrinsic_variables_sub_2016_2025 %>%
  count(season_fct)

##Make predictions
pred_grid_2016_2025_sub <- pred_grid_2016_2025_sub %>%
  mutate(pred_link = predict(mod_gamma_2016_2025,
                             newdata = ., type = "link",
                             re.form = NA, allow.new.levels = TRUE),
         se_link = predict(
           mod_gamma_2016_2025,
           newdata = ., type = "link",
           re.form = NA, se.fit = TRUE, allow.new.levels = TRUE)$se.fit,
         predicted = exp(pred_link),
         conf_lo = exp(pred_link - 1.96 * se_link),
         conf_hi = exp(pred_link + 1.96 * se_link)) %>%
  left_join(seal_years_2016_2025_sub, by = "season_fct") %>%
  group_by(AgeYears, age_cat) %>%
  summarise(predicted = weighted.mean(predicted, n),
            conf_lo = weighted.mean(conf_lo, n),
            conf_hi = weighted.mean(conf_hi, n),
            .groups = "drop")

##Plot predictions with senescent threshold
ggplot(pred_grid_2016_2025_sub, aes(x = AgeYears, y = predicted, color = age_cat, fill = age_cat)) +
  
  #Predictions for individual seasons
  geom_line(aes(group = interaction(age_cat, season_fct)),
            pred_season_2016_2025_sub,
            alpha = 0.3, color = "black") +
  
  #Main ribbon for model w/ CI
  geom_ribbon(data = pred_grid_2016_2025_sub,
              aes(x = AgeYears, ymin = conf_lo, ymax = conf_hi, fill = age_cat),
              alpha = 0.3, color = NA, inherit.aes = FALSE) +
  
  #Thick colored prediction line
  geom_line(data = pred_grid_2016_2025_sub,
            aes(x = AgeYears, y = predicted, color = age_cat),
            linewidth = 1.2, inherit.aes = FALSE) +
  
  #Raw observed points
  geom_pointrange(data = observed_data_2016_2025_sub,
                  aes(x = AgeYears, y = mean_prop, ymin = lwr, ymax = upr),
                  inherit.aes = FALSE, color = "black") +
  
  #Overlay sample size for each age
  geom_text(data = observed_data_2016_2025_sub,
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




gge_overall <- ggpredict(
  mod_piecewise_gamma_2016_2025_both,
  terms = "avg_density",  # only one term -> overall marginal effect
  type = "fixed"             # fixed effects only (population-level)
)

plot(gge_overall) +
  labs(
    x = "Average Density",
    y = "Predicted flipped proportion",
    title = "Overall Effect of Average Density (Gamma GLMM)",
    subtitle = "Marginal effect averaged across all other variables"
  )






##Density effect predictions
gge_overall <- ggpredict(
  mod_gamma_2016_2025,
  terms = "avg_density",  # only one term -> overall marginal effect
  type = "fixed"             # fixed effects only (population-level)
)

# Generate marginal predictions for avg_density
gge_overall <- ggpredict(
  mod_gamma_2016_2025,
  terms = "avg_density",
  type = "fixed"
) %>%
  mutate(
    predicted_prop = mean(predicted),
    conf.low_prop  = sd(conf.low),
    conf.high_prop = sd(conf.high)
  )

# Plot with observed data
ggplot(gge_overall, aes(x = x, y = predicted_prop)) +
  # Confidence ribbon
  geom_ribbon(aes(ymin = conf.low_prop, ymax = conf.high_prop),
              fill = "#7EAAC1", alpha = 0.3) +
  # Model line
  geom_line(color = "#2C7BB6", linewidth = 1.2) +
  # Overlay observed points
  geom_point(data = intrinsic_sub_2016_2025,
             aes(x = avg_density, y = mean(prop_trans)),
             alpha = 0.5, color = "#0B7348", size = 1.5) +
  theme_few() +
  labs(
    x = "Average Density",
    y = "Proportion (<1)"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))








##Gaussian plot
## 1. Observed (raw) data summary
observed_data_gauss_2016_2025 <- intrinsic_sub_2016_2025 %>%
  group_by(AgeYears) %>%
  summarise(
    n = n(),
    mean_prop = mean(prop_trans, na.rm = TRUE),
    sd_prop   = sd(prop_trans, na.rm = TRUE),
    se        = sd_prop / sqrt(n),
    lwr       = mean_prop - 1.96 * se,
    upr       = mean_prop + 1.96 * se,
    .groups   = "drop"
  )

## 2. Per-season random-effect predictions

pred_season_gauss_2016_2025 <- tryCatch({
  ggpredict(
    mod_gauss_2016_2025,
    terms = c("age10 [all]", "age_cat", "season_fct"),
    type  = "random"
  ) %>%
    as_tibble() %>%
    mutate(
      AgeYears   = x * 10 + age_senesce,
      age_cat    = factor(group, levels = c("Young", "Old")),
      season_fct = factor(facet)
    ) %>%
    filter(
      (age_cat == "Young" & AgeYears < age_senesce) |
        (age_cat == "Old"   & AgeYears >= age_senesce)
    )
}, error = function(e) {
  message("⚠️ ggpredict() failed — returning empty tibble.")
  tibble(AgeYears = numeric(), predicted = numeric(), age_cat = factor(), season_fct = factor())
})

## 3. Population-level prediction grid (fixed)

model_levels <- list(
  season_fct  = levels(mod_gauss_2016_2025$frame$season_fct),
  animalID_fct = levels(mod_gauss_2016_2025$frame$animalID_fct)
)

pred_grid_gauss_2016_2025 <- expand.grid(
  AgeYears     = 3:22,
  season_fct   = model_levels$season_fct,
  animalID_fct = model_levels$animalID_fct[1]  # one valid level
) %>%
  mutate(
    age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                     levels = c("Young", "Old")),
    age10 = (AgeYears - age_senesce) / 10,
    avg_density   = mean(intrinsic_sub_2016_2025$avg_density, na.rm = TRUE),
    age_last_seen = mean(intrinsic_sub_2016_2025$age_last_seen, na.rm = TRUE),
    year_born_num = mean(intrinsic_sub_2016_2025$year_born_num, na.rm = TRUE)
  ) %>%
  filter(
    (age_cat == "Young" & AgeYears < age_senesce) |
      (age_cat == "Old"   & AgeYears >= age_senesce)
  )


## 4. Predictions
seal_years_gauss_2016_2025 <- intrinsic_sub_2016_2025 %>% count(season_fct)

pred_grid_gauss_2016_2025 <- pred_grid_gauss_2016_2025 %>%
  mutate(
    pred_fix = predict(mod_gauss_2016_2025, newdata = ., re.form = NA),
    se_val   = tryCatch(
      predict(mod_gauss_2016_2025, newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
      error = function(e) rep(sd(residuals(mod_gauss_2016_2025)), nrow(.))
    ),
    conf_lo = pred_fix - 1.96 * se_val,
    conf_hi = pred_fix + 1.96 * se_val,
    pred_val = predict(mod_gauss_2016_2025, newdata = ., re.form = NULL)
  ) %>%
  left_join(seal_years_gauss_2016_2025, by = "season_fct") %>%
  group_by(AgeYears, age_cat) %>%
  summarise(
    predicted     = weighted.mean(pred_val, n, na.rm = TRUE),
    predicted_pop = weighted.mean(pred_fix, n, na.rm = TRUE),
    conf_lo       = weighted.mean(conf_lo, n, na.rm = TRUE),
    conf_hi       = weighted.mean(conf_hi, n, na.rm = TRUE),
    .groups = "drop"
  )


# Back-transform predictions to proportion scale
invlogit <- function(x) 1 / (1 + exp(-x))

pred_grid_gauss_2016_2025 <- pred_grid_gauss_2016_2025 %>%
  mutate(
    predicted_prop     = invlogit(predicted_pop),
    conf_lo_prop        = invlogit(conf_lo),
    conf_hi_prop        = invlogit(conf_hi)
  )

pred_season_gauss_2016_2025 <- pred_season_gauss_2016_2025 %>%
  mutate(predicted_prop = invlogit(predicted))

observed_data_gauss_2016_2025 <- observed_data_gauss_2016_2025 %>%
  mutate(
    mean_prop_back = invlogit(mean_prop),
    lwr_back       = invlogit(lwr),
    upr_back       = invlogit(upr)
  )

## PLOT on proportion scale
y_offset <- 1.05 * max(observed_data_gauss_2016_2025$upr_back, na.rm = TRUE)

ggplot() +
  
  # per-season random lines (on proportion scale)
  geom_line(
    data = pred_season_gauss_2016_2025,
    aes(x = AgeYears, y = predicted_prop,
        group = interaction(season_fct, age_cat)),
    color = "gray60", alpha = 0.6, linewidth = 0.6
  ) +
  
  # confidence ribbon (proportion)
  geom_ribbon(
    data = pred_grid_gauss_2016_2025,
    aes(x = AgeYears, ymin = conf_lo_prop, ymax = conf_hi_prop, fill = age_cat),
    alpha = 0.3, color = NA
  ) +
  
  # main line
  geom_line(
    data = pred_grid_gauss_2016_2025,
    aes(x = AgeYears, y = predicted_prop, color = age_cat),
    linewidth = 1.3
  ) +
  
  # observed data (means ± SE)
  geom_pointrange(
    data = observed_data_gauss_2016_2025,
    aes(x = AgeYears, y = mean_prop_back, ymin = lwr_back, ymax = upr_back),
    color = "black",
    linewidth = 0.4
  ) +
  
  geom_text(
    data = observed_data_gauss_2016_2025,
    aes(x = AgeYears, label = n),
    y = y_offset,   # <- outside aes()
    size = 3
  ) +
  
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_y_continuous(
    name = "Predicted proportion (0–1 scale)",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_x_continuous("Age (Years)", n.breaks = 20) +
  scale_color_manual(values = c("Young" = "#66C2A5", "Old" = "#D5A5C9")) +
  scale_fill_manual(values  = c("Young" = "#66C2A5", "Old" = "#D5A5C9")) +
  theme_few() +
  labs(
    color = "Age class",
    fill  = "Age class",
    title = "Gaussian GLMM: Age effect on proportion (2016–2025)",
    subtitle = "Back-transformed to proportion scale (0–1)"
  )



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
  is_one ~ age10 + age_last_seen + total_resights + scale(year_born_num) + 
    (1 | animalID_fct) + (1 | season_fct),
  grname   = "animalID_fct", #grouping by individual to get influence of animalID
  data     = intrinsic_variables,
  link     = "logit",
  nboot    = 1000,
  npermut  = 0
); summary(rpt_all) #summarize results


# Repeatability for young individuals from all years, for is_one binomial 
rpt_young <- rptBinary( #Use Binary version of rpt
  is_one ~ age10 + age_last_seen + total_resights + scale(year_born_num) + 
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
  filter(n_seasons >= 3) # require at least two repeated breeding seasons per individual

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
  geom_point(aes(size = n_obs_total, alpha = 0.8)) +
  geom_vline(xintercept = c(0.2, 0.8), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c(
    "Maternal specialist" = "#1f78b4",
    "Non-maternal specialist" = "#B74F6F",
    "Plastic" = "#6F73D2",
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


n <- nrow(season_level)
season_level <- season_level %>%
  mutate(prop_adj = (proportion * (n - 1) + 0.5) / n)
library(brms)
fit_beta <- brm(
  bf(proportion ~ 1 + (1|animalID),          # mean (mu)
     phi ~ 1 + (1|animalID)),                # precision (inverse of variance)
  data   = season_level,
  family = zero_one_inflated_beta(),                   # NOT gaussian()
  cores = 4, chains = 4
)

