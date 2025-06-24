##GLMM modeling for intrinsic variables
library(glmmTMB)
library(splines)
library(DHARMa)
library(ggeffects)
library(dplyr)

##Original GLMM Gamma distribution with polynomial Age
mod_01_intrinsic_glmm <- glmmTMB(proportion ~ poly(AgeYears, 3) + year_born_fct + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
  data = intrinsic_variables,
  family = Gamma(link = "log"))
summary(mod_01_intrinsic_glmm)

##Trying flipped proportion to make it more gamma-friendly
flipped_prop <- max(intrinsic_variables$proportion) - intrinsic_variables$proportion
flipped_prop <- flipped_prop - min(flipped_prop) + 0.001  # ensure strictly positive
intrinsic_variables$flipped_prop <- flipped_prop

##Improved GLMM with natural splines Age
mod_02_intrinsic_glmm <- glmmTMB(flipped_prop ~ ns(AgeYears, df = 3) + year_born_fct + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
  data = intrinsic_variables,
  family = Gamma(link = "log"))
summary(mod_02_intrinsic_glmm)
##Residuals
DHARMa::simulateResiduals(mod_02_intrinsic_glmm, plot = TRUE)
#even with flipped proportion, the model's residuals aren't great

####Trying a two-model approach; binomial model using bernoulli prop == 1, and gamma model for proportions under 1

##First, modify intrinsic_variables so it has a column for proportion = 1 (1) OR not 1 (0)
intrinsic_variables$is_one <- as.numeric(intrinsic_variables$proportion == 1)

##Binomial model for proportion == 1
mod_binom_1_only <- glmmTMB(is_one ~ ns(AgeYears, 3) + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
  data = intrinsic_variables,
  family = binomial(link = "logit"))
summary(mod_binom_1_only)
DHARMa::simulateResiduals(mod_binom_1_only, plot = TRUE)
#residuals look great

##Use ggeffects to predict the likelihood that proportion is 1 based on age
pred_age_mod_binom_1 <- ggpredict(mod_binom_1_only, terms = "AgeYears [all]")

##Plot the effect of age on probability of proportion == 1
ggplot(pred_age_mod_binom_1, aes(x = x, y = predicted)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) + #confidence intervals
  labs(title = "Predicted Probability of Always 1-Pup Association vs. Age",
    x = "Age (Years)",
    y = "Predicted Probability (is_one = 1)") +
  theme_minimal()
#looks like a peak/senescent effect but very wide error at older ages

##A gamma model for all proportions < 1
#Create a version of intrinsic_variables with proportion < 1
intrinsic_variables_sub <- intrinsic_variables %>%
filter(proportion < 1)

##Transform proportion for this filtered data to make it suitable for gamma distribution 
flipped_prop <- max(intrinsic_variables_sub$proportion) - intrinsic_variables_sub$proportion
flipped_prop <- flipped_prop - min(flipped_prop) + 0.001
intrinsic_variables_sub$flipped_prop <- flipped_prop

##Gamma model <1
mod_gamma_non1 <- glmmTMB(flipped_prop ~ ns(AgeYears, 3) + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
  data = intrinsic_variables_sub,
  family = Gamma(link = "log"))
summary(mod_gamma_non1)
DHARMa::simulateResiduals(mod_gamma_non1, plot = TRUE)
#residuals look ok

##Exploring senescence thresholds and age class effects

##Trying an age effect threshold piecewise model
age_senesce <- 13 #setting senescence threshold at 11 (Allison paper!)

##Setting threshold in data
intrinsic_variables <- intrinsic_variables %>%
mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
       levels = c("Young", "Old"))) %>%
mutate(age10 = (AgeYears - age_senesce) / 10) #Scaled numeric version of Age centered at senescence threshold

##Piecewise model for age pre- and post- senescent threshold
mod_piecewise <- glmmTMB(is_one ~ age10 : age_cat + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
  family = binomial(link = "logit"),
  data = intrinsic_variables)
summary(mod_piecewise)
DHARMa::simulateResiduals(mod_piecewise, plot = TRUE)
#residuals look good here too

##Need a prediction grid to visualize the effect
prediction_grid <- expand.grid(
  AgeYears = seq(4, 20, 0.1)) %>%
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old")),
    age10 = (AgeYears - age_senesce) / 10,
    age_last_seen = mean(intrinsic_variables$age_last_seen),
    total_resights = mean(intrinsic_variables$total_resights),
    animalID = factor(NA, levels = levels(intrinsic_variables$animalID)),
    season = factor(NA, levels = levels(intrinsic_variables$season)))

##Make predictions
prediction_grid$predicted <- predict(mod_piecewise,
                             newdata = prediction_grid,
                             type = "response",
                             re.form = NA)

##Plot predictions with senescent threshold
ggplot(prediction_grid, aes(x = AgeYears, y = predicted, color = age_cat)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_y_continuous("Predicted probability of 1") +
  scale_x_continuous("Age (Years)") +
  theme_minimal() +
  labs(color = "Age class")

##Age class interaction model
#First binning age classes as "young", "prime", "senescent"
intrinsic_variables <- intrinsic_variables %>%
  mutate(AgeClass = case_when(
      AgeYears >= 3 & AgeYears <= 6 ~ "Young",
      AgeYears >= 7 & AgeYears <= 13 ~ "Prime",
      AgeYears >= 14 ~ "Senescent"),
    AgeClass = factor(AgeClass, levels = c("Young", "Prime", "Senescent")))


##Model age class interaction with age
mod_age_class <- glmmTMB(
  is_one ~ AgeYears * AgeClass + age_last_seen + total_resights + (1 | animalID_fct) + (1 | season_fct),
  family = binomial(link = "logit"),
  data = intrinsic_variables)
summary(mod_age_class)
DHARMa::simulateResiduals(mod_age_class, plot = TRUE)
#residuals look good

##Generate predicted age and age class effects on probability of prop == 1
pred_mod_age_class <- ggpredict(mod_age_class, terms = c("AgeYears", "AgeClass"))

##Restrict so that it only shows predictions within each age class
pred_mod_age_class_filtered <- pred_mod_age_class %>%
  filter((group == "Young" & x >= 3 & x <= 6) |
           (group == "Prime" & x >= 7 & x <= 13) |
           (group == "Senescent" & x >= 14))

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
