##GLMM modeling for intrinsic variables
library(glmmTMB)
library(splines)
library(DHARMa)
library(ggeffects)
library(dplyr)
library(performance)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

################################## Old Models ########################################

##GLMM Gamma model with polynomial Age
mod_01_intrinsic_glmm <- glmmTMB(flipped_prop ~ poly(AgeYears, 3) + year_born_fct + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
  data = intrinsic_variables,
  family = Gamma(link = "log"))
summary(mod_01_intrinsic_glmm)
DHARMa::simulateResiduals(mod_01_intrinsic_glmm, plot = TRUE)

##Trying flipped proportion to make it more gamma-friendly
flipped_prop <- max(intrinsic_variables$proportion) - intrinsic_variables$proportion #subtract prop from max prop
flipped_prop <- flipped_prop - min(flipped_prop) + 0.001  #ensure strictly positive
intrinsic_variables$flipped_prop <- flipped_prop 

##Improved GLMM with natural splines Age
mod_02_intrinsic_glmm <- glmmTMB(flipped_prop ~ ns(AgeYears, 3) + year_born_fct + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
  data = intrinsic_variables,
  family = Gamma(link = "log"))
summary(mod_02_intrinsic_glmm)
DHARMa::simulateResiduals(mod_02_intrinsic_glmm, plot = TRUE)
#even with flipped proportion, the model's residuals aren't great
check_collinearity(mod_02_intrinsic_glmm)
#moderate co-linearity in age

######################## Binomial model for proportion == 1 ################################

##First, modify intrinsic_variables so it has a column for proportion = 1 (1) OR not 1 (0)
intrinsic_variables <- intrinsic_variables %>%
  mutate(is_one = as.numeric(proportion == 1))


##Binomial model for proportion == 1
mod_binom_1_only <- glmmTMB(is_one ~ ns(AgeYears, 3) + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
  data = intrinsic_variables,
  family = binomial(link = "logit"))
summary(mod_binom_1_only)
exp(fixef(mod_binom_1_only)$cond) 
DHARMa::simulateResiduals(mod_binom_1_only, plot = TRUE)
#residuals look great
check_collinearity(mod_binom_1_only)

##Use ggeffects to predict the likelihood that proportion is 1 based on age
pred_age_mod_binom <- ggpredict(mod_binom_1_only,
                                terms = "AgeYears [all]")

pred_total_resights_mod_binom <- ggpredict(mod_binom_1_only,
                                           terms = "total_resights")

##Plot the effect of age on probability of proportion == 1
ggplot(pred_age_mod_binom, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1.2, color = "#2A5EA7") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "#7EAAC1") + #confidence intervals
  labs(title = "Predicted Probability of Perfect 1-Pup Association vs. Age",
    x = "Age (Years)",
    y = "Predicted Probability of Proportion = 1") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 20) +
  theme_few()
#looks like a peak/senescent effect but very wide error at older ages

##Plot the effect of total resights on probability of prop == 1
ggplot(pred_total_resights_mod_binom, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1.2, color = "#2A5EA7") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "#7EAAC1") + #confidence intervals
  labs(title = "Effect of Total Resights on Probability of Perfect Proportion",
       x = "Total Resights",
       y = "Predicted Probability of Proportion = 1") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 20) +
  theme_few()
#more resights = lower probability that prop is 1

######################### Gamma model for all proportions < 1 ###############################

#Create a version of intrinsic_variables with proportion < 1
intrinsic_variables_sub <- intrinsic_variables %>%
filter(proportion < 1)

##Transform proportion for this filtered data to make it suitable for gamma distribution 
flipped_prop <- max(intrinsic_variables_sub$proportion) - intrinsic_variables_sub$proportion #subtract all proportions from the maximum
flipped_prop <- flipped_prop - min(flipped_prop) + 0.0000001 #ensure no 0s by adding small constant
intrinsic_variables_sub$flipped_prop <- flipped_prop

##Histogram of flipped proportion
ggplot(data = intrinsic_variables_sub, aes(x = flipped_prop)) +
geom_histogram(binwidth = 0.01, fill = "skyblue", color = "darkblue") +
  labs(title = "Histogram of Flipped Proportion MPA <1",
       x = "Flipped Proportion MPA", y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_few()

##Gamma model <1
mod_gamma_non1 <- glmmTMB(flipped_prop ~ ns(AgeYears, 3) + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
  data = intrinsic_variables_sub,
  family = Gamma(link = "log"))
summary(mod_gamma_non1)
exp(fixef(mod_gamma_non1)$cond) 
DHARMa::simulateResiduals(mod_gamma_non1, plot = TRUE)
#residuals look good
check_collinearity(mod_gamma_non1)

#Predict proportion across age (see if it makes sense)
pred_gamma <- ggpredict(mod_gamma_non1,
                        terms = "AgeYears [all]")

##Plot the effect of age on proportion
ggplot(pred_gamma, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1.2, color = "#2A5EA7") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "#7EAAC1") +
  labs(title = "Predicted Flipped Proportion by Age",
       x = "Age (Years)",
       y = "Predicted Flipped Proportion") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 20) +
  theme_few()

######################## Piecewise segmented regression ##################################

age_senesce <- 11 #setting senescence threshold at 11 (Allison paper!)

##Setting threshold in data
intrinsic_variables <- intrinsic_variables %>%
mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
       levels = c("Young", "Old"))) %>%
mutate(age10 = (AgeYears - age_senesce) / 10) #Scaled numeric version of Age centered at senescence threshold

##A table for raw observed proportion to overlay on plot
observed_data <- intrinsic_variables %>%
  group_by(AgeYears) %>%
  summarize(n = n(),
            n_one = sum(is_one),
            perc_one = n_one / n,
            lwr = binom.test(n_one, n)$conf.int[1],
            upr = binom.test(n_one, n)$conf.int[2])

##Piecewise model for age pre- and post- senescent threshold
mod_piecewise <- glmmTMB(is_one ~ age10 : age_cat + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
  family = binomial(link = "logit"),
  data = intrinsic_variables)
summary(mod_piecewise)
exp(fixef(mod_piecewise)$cond) 
DHARMa::simulateResiduals(mod_piecewise, plot = TRUE)
#residuals look good here too
check_collinearity(mod_piecewise)

##Create ANNUAL (i.e. rand effect of season only) predictions
pred_season <- ggpredict(
  mod_piecewise, 
  terms = c("age10 [all]", "age_cat", "season_fct"),
  type = "random"
) %>% 
  as_tibble() %>% 
  mutate(AgeYears = x * 10 + age_senesce,
         age_cat = factor(group, levels = c("Young", "Old")),
         season_fct = factor(facet)) %>% 
  filter((age_cat == "Young" & AgeYears < age_senesce) |
           (age_cat == "Old" & AgeYears >= age_senesce))

##Need a prediction grid to visualize the Age effect
pred_grid <- expand.grid(
  AgeYears = 3:22,
  season_fct = levels(intrinsic_variables$season_fct)
) %>%
  mutate(
    age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                     levels = c("Young", "Old")),
    age10 = (AgeYears - age_senesce) / 10,
    total_resights = mean(intrinsic_variables$total_resights),
    age_last_seen = mean(intrinsic_variables$age_last_seen),
    animalID_fct = factor(NA, levels = levels(intrinsic_variables$animalID_fct))
  )

seal_years <- intrinsic_variables %>%
  count(season_fct)

##Make predictions
pred_grid <- pred_grid %>%
  mutate(pred_logit = predict(mod_piecewise, newdata = pred_grid, re.form = NULL),
         se_logit = predict(mod_piecewise, newdata = pred_grid, re.form = NA, se.fit = TRUE)$se.fit,
         pred_fixed = predict(mod_piecewise, newdata = pred_grid, re.form = NA)) %>%
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
  
  # Predictions for individual seasons
  geom_line(aes(group = interaction(age_cat, season_fct)),
            pred_season,
            alpha = 0.3, color = "#7EAAC1") +
  
  # Main ribbon for model w/ CI
  geom_ribbon(data = pred_grid,
              aes(x = AgeYears, ymin = conf_lo, ymax = conf_hi, fill = age_cat),
              alpha = 0.3, color = NA, inherit.aes = FALSE) +
  
  # Thick colored prediction line
  geom_line(data = pred_grid,
            aes(x = AgeYears, y = predicted, color = age_cat),
            linewidth = 1.2, inherit.aes = FALSE) +
  
  # Raw observed points
  geom_pointrange(data = observed_data,
                  aes(x = AgeYears, y = perc_one, ymin = lwr, ymax = upr),
                  inherit.aes = FALSE, color = "black") +
  
  # Overlay sample size for each age
  geom_text(data = observed_data,
            aes(x = AgeYears, y = 1.03, label = n),
            inherit.aes = FALSE, size = 3) +
  
  # Formatting, colors
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_y_continuous(n.breaks = 5) +
  scale_x_continuous(n.breaks = 20) +
  scale_color_manual(values = c("Young" = "#66C2A5", "Old" = "#D5A5C9")) +
  scale_fill_manual(values = c("Young" = "#66C2A5", "Old" = "#D5A5C9")) +
  labs(x = "Age (Years)", y = "Probability of 1", color = "Age class", fill = "Age class") +
  theme_few()


######################## Age class (Young, Prime, Senescent) model ###############################

#First binning age classes as "young", "prime", "senescent"
intrinsic_variables <- intrinsic_variables %>%
  mutate(AgeClass = case_when(
    AgeYears >= 3 & AgeYears <= 6 ~ "Young",
    AgeYears >= 7 & AgeYears <= 11 ~ "Prime",
    AgeYears >= 12 ~ "Senescent"),
    AgeClass = factor(AgeClass, levels = c("Young", "Prime", "Senescent")))

##Model age class interaction with age
mod_age_class <- glmmTMB(
  is_one ~ AgeYears : AgeClass + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
  family = binomial(link = "logit"),
  data = intrinsic_variables)
summary(mod_age_class)
exp(fixef(mod_age_class)$cond) 
DHARMa::simulateResiduals(mod_age_class, plot = TRUE)
#residuals look good

##Generate predicted age and age class effects on probability of prop == 1
pred_mod_age_class <- ggpredict(mod_age_class, 
                                terms = c("AgeYears", "AgeClass"))

##Restrict so that it only shows predictions within each age class
pred_mod_age_class_filtered <- pred_mod_age_class %>%
  filter((group == "Young" & x >= 3 & x <= 6) |
           (group == "Prime" & x >= 7 & x <= 11) |
           (group == "Senescent" & x >= 12))

##Plot effect of age and age class on prob of prop == 1
ggplot(pred_mod_age_class_filtered, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(title = "Predicted Perfect 1-Pup Association by Age and Age Class",
    x = "Age (Years)",
    y = "Predicted Probability",
    color = "Age Class",
    fill = "Age Class") +
  theme_minimal()


##################################### Models subset 2016 - 2025 ########################################

#First, subset the data for only seasons 2016 - 2025
intrinsic_2016_2025 <- intrinsic_variables %>%
  filter(season >= 2016)

#Histogram of filtered data's proportions
ggplot(data = intrinsic_2016_2025, aes(x = proportion)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Proportion MPA",
       x = "Proportion MPA", y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  theme_few()

##Binomial model 2016-2025
mod_binom_1_only_2016_2025 <- glmmTMB(is_one ~ ns(AgeYears, 3) + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                                      data = intrinsic_2016_2025,
                                      family = binomial(link = "logit"))
summary(mod_binom_1_only_2016_2025)
exp(fixef(mod_binom_1_only_2016_2025)$cond)
DHARMa::simulateResiduals(mod_binom_1_only_2016_2025, plot = TRUE)

#Predict the effects of age using the model
pred_age_mod_binom_2016_2025 <- ggpredict(mod_binom_1_only_2016_2025, 
                                  terms = "AgeYears [all]")

##Plot the effect of age on probability of proportion == 1
ggplot(pred_age_mod_binom_2016_2025, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) + #confidence intervals
  labs(title = "Predicted Probability of Perfect 1-Pup Association vs. Age",
       x = "Age (Years)",
       y = "Predicted Probability of Proportion = 1") +
  theme_minimal()

##Gamma model 2016-2025
intrinsic_sub_2016_2025 <- intrinsic_variables_sub %>%
  filter(season >= 2016)

mod_gamma_non1_2016_2025 <- glmmTMB(flipped_prop ~ ns(AgeYears, 3) + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                                    data = intrinsic_sub_2016_2025,
                                    family = Gamma(link = "log"))
summary(mod_gamma_non1_2016_2025)
exp(fixef(mod_gamma_non1_2016_2025)$cond) 
DHARMa::simulateResiduals(mod_gamma_non1_2016_2025, plot = TRUE)

##Piecewise model for age pre- and post- senescent threshold

#Make raw data table to overlay points on the plot

observed_data_2016_2025 <- intrinsic_2016_2025 %>%
  group_by(AgeYears) %>%
  summarize(n = n(), 
    n_one = sum(is_one),
    perc_one = n_one / n, #Summarize the proportion of 1s by AgeYears
    lwr = binom.test(n_one, n)$conf.int[1], #Confidence intervals
    upr = binom.test(n_one, n)$conf.int[2])

mod_piecewise_2016_2025 <- glmmTMB(is_one ~ age10 : age_cat + total_resights + (1 | animalID_fct) + (1 | season_fct),
                                  family = binomial(link = "logit"),
                                  data = intrinsic_2016_2025)
summary(mod_piecewise_2016_2025)
exp(fixef(mod_piecewise_2016_2025)$cond) 
DHARMa::simulateResiduals(mod_piecewise_2016_2025, plot = TRUE)
#residuals look good here too

##Need a prediction grid to visualize the effect
pred_grid_2016_2025 <- expand.grid(AgeYears = 3:22) %>%
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old")),
         age10 = (AgeYears - age_senesce) / 10,
         total_resights = mean(intrinsic_2016_2025$total_resights),
         animalID_fct = factor(NA, levels = levels(intrinsic_2016_2025$animalID)),
         season_fct = factor(NA, levels = levels(intrinsic_2016_2025$season)))

#make predictions
pred_grid_2016_2025$predicted <- predict(mod_piecewise_2016_2025, 
                                     newdata = pred_grid_2016_2025, 
                                     type = "response", 
                                     re.form = NA)

##Plot predictions with senescent threshold
ggplot(pred_grid_2016_2025, aes(x = AgeYears, y = predicted, color = age_cat)) +
  geom_line(linewidth = 1.2) +
  
  #Add senescence threshold line
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  
  # Add observed points and confidence intervals
  geom_pointrange(data = observed_data_2016_2025,
                  aes(x = AgeYears, y = perc_one, ymin = lwr, ymax = upr),
                  inherit.aes = FALSE,
                  size = 0.4) +
  
  # Add sample size labels above the points
  geom_text(data = observed_data,
            aes(x = AgeYears, y = 1.02, label = n),
            inherit.aes = FALSE,
            size = 3,
            vjust = -0.5) +
  scale_y_continuous("Predicted probability of 1") +
  scale_x_continuous("Age (Years)") +
  theme_minimal() +
  labs(color = "Age class")


##Need a prediction grid to visualize the effect
pred_grid_2016_2025 <- expand.grid(AgeYears = 3:22,
                                   season_fct = levels(mod_piecewise_2016_2025$frame$season_fct)) %>%
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"), levels = c("Young", "Old")),
         age10 = (AgeYears - age_senesce) / 10,
         age_last_seen = mean(intrinsic_2016_2025$age_last_seen),
         total_resights = mean(intrinsic_2016_2025$total_resights),
         animalID_fct = factor(NA, levels = levels(intrinsic_2016_2025$animalID_fct)))

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

pred_grid_season_2016_2025 <- expand.grid(AgeYears = 3:22,
                                          season_fct = levels(mod_piecewise_2016_2025$frame$season_fct)) %>%
  mutate(season_fct = factor(season_fct, levels = levels(intrinsic_variables$season_fct)),
         age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"), levels = levels(intrinsic_variables$age_cat)),
         age10 = (AgeYears - age_senesce) / 10,
         age_last_seen = mean(intrinsic_2016_2025$age_last_seen, na.rm = TRUE),
         total_resights = mean(intrinsic_2016_2025$total_resights, na.rm = TRUE),
         animalID_fct = factor(NA, levels = levels(intrinsic_2016_2025$animalID_fct)))

pred_grid_season_2016_2025 <- pred_grid_season_2016_2025 %>%
  mutate(predicted = predict(mod_piecewise_2016_2025, newdata = pred_grid_season_2016_2025, re.form = NULL)) %>%
  mutate(predicted_prob = plogis(predicted))

##Plot predictions with senescent threshold
ggplot(pred_grid_2016_2025, aes(x = AgeYears, y = predicted, color = age_cat, fill = age_cat)) +
  
  # Per-season prediction lines (thin gray)
  geom_line(data = pred_grid_season_2016_2025,
            aes(x = AgeYears, y = predicted_prob, group = season_fct),
            color = "gray70", alpha = 0.2, linewidth = 0.4, inherit.aes = FALSE) +
  
  # Main ribbon for weighted model CI
  geom_ribbon(data = pred_grid_2016_2025,
              aes(x = AgeYears, ymin = conf_lo, ymax = conf_hi, fill = age_cat),
              alpha = 0.3, color = NA, inherit.aes = FALSE) +
  
  # Thick colored prediction line (weighted)
  geom_line(data = pred_grid_2016_2025,
            aes(x = AgeYears, y = predicted, color = age_cat),
            linewidth = 1.2, inherit.aes = FALSE) +
  
  # Dotted line: fixed effects only
  geom_line(data = pred_grid_2016_2025,
            aes(x = AgeYears, y = predicted_pop),
            linetype = "dotted", color = "gray40", linewidth = 1) +
  
  # Raw observed points
  geom_pointrange(data = observed_data_2016_2025,
                  aes(x = AgeYears, y = perc_one, ymin = lwr, ymax = upr),
                  inherit.aes = FALSE, color = "black") +
  geom_text(data = observed_data_2016_2025,
            aes(x = AgeYears, y = 1.03, label = n),
            inherit.aes = FALSE, size = 3) +
  
  # Formatting, colors
  geom_vline(xintercept = age_senesce, linetype = "dashed") +
  scale_y_continuous("Probability of 1", limits = c(0, 1.05)) +
  scale_color_manual(values = c("Young" = "#66c2a5", "Old" = "#fc8d62")) +
  scale_fill_manual(values = c("Young" = "#66c2a5", "Old" = "#fc8d62")) +
  labs(x = "Age (Years)", color = "Age class", fill = "Age class") +
  theme_classic()
  












######################### Model with Birth Date Timing #########################

library(lubridate)
intrinsic_variables <- intrinsic_variables %>%
  mutate(BirthDOY = yday(as.Date(paste0("2023-", BirthDate))))

early_end <- yday(as.Date("2024-01-15"))    # 5
early_start <- yday(as.Date("2024-12-13"))  # 348
peak_start <- yday(as.Date("2024-01-15"))   # 6
peak_end <- yday(as.Date("2024-02-05"))     # 20

intrinsic_variables <- intrinsic_variables %>%
  mutate(birth_group = case_when(
    BirthDOY <= 10 ~ "Early",    # Jan 1-10
    BirthDOY <= 30 ~ "Peak",     # Jan 11-20
    TRUE ~ "Late"                # Jan 21 onwards
  ))

table(intrinsic_variables$birth_group)

mod_binom <- glmmTMB(is_one ~ ns(AgeYears, 3) + birth_group + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                          data = intrinsic_variables,
                          family = binomial())
summary(mod_binom)
exp(fixef(mod_gamma_non1)$cond) 
DHARMa::simulateResiduals(mod_gamma_non1, plot = TRUE)

################################# 2+ pups model ####################################


intrinsic_variables_2 <- metadata %>%
  select(animalID, season, AgeYears, withpup) %>%
  filter(withpup != '?') %>%
  group_by(animalID) %>%
  mutate(pup_more_than_2 = ifelse(any(withpup >= 2), 1, 0)) %>%  #Check if any year for this animal has withpup >= 2
  ungroup() %>%
  distinct(animalID, .keep_all = TRUE) %>% #Keep only one row per animalID
  mutate(animalID_fct = factor(animalID), season_fct = factor(season))

mod_binom_greater_2 <- glmmTMB(pup_more_than_2 ~ AgeYears + (1 | animalID_fct) + (1 | season_fct),
                           data = intrinsic_variables_2,
                           family = binomial(link = "logit"))
summary(mod_binom_greater_2)


pred_age_greater_2 <- ggpredict(mod_binom_greater_2,
                                terms = "AgeYears [all]")


ggplot(pred_age_greater_2, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1.2, color = "#2A5EA7") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "#7EAAC1") + #confidence intervals
  labs(title = "Predicted probability of at least 1 2+ Pup Sighting vs. Age",
       x = "Age (Years)",
       y = "Predicted Probability of 2+ Sighting") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 20) +
  theme_few()






