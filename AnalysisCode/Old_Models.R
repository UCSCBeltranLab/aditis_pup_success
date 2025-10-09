
######################## Spline binomial model for proportion == 1 ################################

#Trying different versions with random effects for AIC and residuals

mod_binom_season_ID <- glmmTMB(is_one ~ ns(AgeYears, 3) + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                               data = intrinsic_variables,
                               family = binomial(link = "logit"))
summary(mod_binom_season_ID)
exp(fixef(mod_binom_season_ID)$cond) 
DHARMa::simulateResiduals(mod_binom_season_ID, plot = TRUE)
#residuals look fine
check_collinearity(mod_binom_season_ID)

##trying multiple versions with different combinations of random effects

mod_binom_yr_born_ID <- glmmTMB(is_one ~ ns(AgeYears, 3) + age_last_seen + total_resights + (1 | animalID_fct) + (1 | year_born_fct),
                                data = intrinsic_variables,
                                family = binomial(link = "logit"))
summary(mod_binom_yr_born_ID)
exp(fixef(mod_binom_yr_born_ID)$cond) 
DHARMa::simulateResiduals(mod_binom_yr_born_ID, plot = TRUE)
#residuals look great
check_collinearity(mod_binom_1_yr_born_ID)

mod_binom_1_all <- glmmTMB(is_one ~ ns(AgeYears, 3) + age_last_seen + total_resights + (1 | animalID_fct) + (1 | year_born_fct) + (1 | season_fct),
                           data = intrinsic_variables,
                           family = binomial(link = "logit"))
summary(mod_binom_1_all)
exp(fixef(mod_binom_1_all)$cond) 
DHARMa::simulateResiduals(mod_binom_1_all, plot = TRUE)
#residuals look great
check_collinearity(mod_binom_1_all)

##Use ggeffects to predict the likelihood that proportion is 1 based on age
pred_age_mod_binom <- ggpredict(mod_binom_all,
                                terms = "AgeYears [all]")

pred_total_resights_mod_binom <- ggpredict(mod_binom_all,
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
#more resights = lower probability that MPA is 1

######################### Spline gamma model for all proportions < 1 ###############################

#Trying different versions with random effects for AIC and residuals

mod_gamma_all <- glmmTMB(flipped_prop ~ ns(AgeYears, 3) + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct) + (1 | year_born_fct),
                         data = intrinsic_variables_sub,
                         family = Gamma(link = "log"))
summary(mod_gamma_all)
exp(fixef(mod_gamma_all)$cond) 
DHARMa::simulateResiduals(mod_gamma_all, plot = TRUE)
#residuals look good
check_collinearity(mod_gamma_all)

mod_gamma_season_ID <- glmmTMB(flipped_prop ~ ns(AgeYears, 3) + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                               data = intrinsic_variables_sub,
                               family = Gamma(link = "log"))
summary(mod_gamma_season_ID)
exp(fixef(mod_gamma_season_ID)$cond) 
DHARMa::simulateResiduals(mod_gamma_season_ID, plot = TRUE)
#residuals look good
check_collinearity(mod_gamma_season_ID)

mod_gamma_ID <- glmmTMB(flipped_prop ~ ns(AgeYears, 3) + age_last_seen + total_resights + (1 | animalID_fct),
                        data = intrinsic_variables_sub,
                        family = Gamma(link = "log"))
summary(mod_gamma_ID)
exp(fixef(mod_gamma_ID)$cond) 
DHARMa::simulateResiduals(mod_gamma_ID, plot = TRUE)
#residuals look good
check_collinearity(mod_gamma_ID)

#Predict proportion across age (just to see if it makes sense)
pred_gamma <- ggpredict(mod_gamma_all,
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

################################## Old Models ########################################

##Trying flipped proportion to make it more gamma-friendly
flipped_prop <- max(intrinsic_variables$proportion) - intrinsic_variables$proportion #subtract prop from max prop
flipped_prop <- flipped_prop - min(flipped_prop) + 0.001  #ensure strictly positive
intrinsic_variables$flipped_prop <- flipped_prop 

##GLMM Gamma model with polynomial Age
mod_01_intrinsic_glmm <- glmmTMB(flipped_prop ~ poly(AgeYears, 3) + year_born_fct + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                                 data = intrinsic_variables,
                                 family = Gamma(link = "log"))
summary(mod_01_intrinsic_glmm)
DHARMa::simulateResiduals(mod_01_intrinsic_glmm, plot = TRUE)

##Improved GLMM with natural splines Age
mod_02_intrinsic_glmm <- glmmTMB(flipped_prop ~ ns(AgeYears, 3) + year_born_fct + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
                                 data = intrinsic_variables,
                                 family = Gamma(link = "log"))
summary(mod_02_intrinsic_glmm)
DHARMa::simulateResiduals(mod_02_intrinsic_glmm, plot = TRUE)
#even with flipped proportion, the model's residuals aren't great
check_collinearity(mod_02_intrinsic_glmm)
#moderate co-linearity in age


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

## First clean and convert withpup
metadata_clean <- metadata %>%
  filter(!grepl("\\?", withpup)) %>%
  mutate(withpup = as.numeric(withpup)) %>%
  filter(!is.na(withpup))

## Second make pup_more_than_2 per animalID and season
intrinsic_variables_2 <- metadata_clean %>%
  select(animalID, season, AgeYears, withpup) %>%
  group_by(animalID, season) %>%
  summarise(AgeYears = first(AgeYears),
            pup_more_than_2 = ifelse(any(withpup >= 2), 1, 0)) %>%
  ungroup() %>%
  mutate(animalID_fct = factor(animalID),
         season_fct = factor(season))

mod_binom_greater_2 <- glmmTMB(pup_more_than_2 ~ AgeYears + (1 | animalID_fct) + (1 | season_fct),
                               data = intrinsic_variables_2,
                               family = binomial(link = "logit"))
summary(mod_binom_greater_2)

pred_age_greater_2 <- ggpredict(mod_binom_greater_2,
                                terms = "AgeYears [all]")


ggplot(pred_age_greater_2, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1.2, color = "#2A5EA7") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "#7EAAC1") +
  labs(title = "Predicted probability of at least 1 2+ Pup Sighting vs. Age",
       x = "Age (Years)",
       y = "Predicted Probability of 2+ Sighting") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 20) +
  theme_few()