library(glmmTMB)
library(splines)

##Here I'm going to use raw count data rather than proportion

ggplot(intrinsic_count_version, aes(x = count_1_pup / total_resights)) +
  geom_histogram(binwidth = 0.05) +
  labs(title = "Distribution of Proportion of Days with 1 Pup",
       x = "Proportion of days seen with exactly 1 pup",
       y = "Number of observations")

intrinsic_count_version <- metadata %>%
  select(animalID, season, AgeYears, BirthDate, year_born, count_1_pup, total_resights) %>%
  distinct() %>% #eliminate the metadata multiple rows; we only need one per animalID per season
  filter(!is.na(count_1_pup)) %>%
  filter(!is.na(total_resights)) %>%
  group_by(animalID) %>%  
  mutate(year_born_fct = factor(year_born), #make year_born a factor in separate column
         animalID_fct = factor(animalID), #make animalID a factor in separate column
         season_fct = factor(season),
         age_last_seen = max(AgeYears, na.rm = TRUE)) %>% #calculate age at last observation
  ungroup()

#modeling using raw counts
mod_binom_count_version <- glmmTMB(
  cbind(count_1_pup, total_resights - count_1_pup) ~ 
    ns(AgeYears, df = 3) + age_last_seen + (1 | animalID_fct) + (1 | season_fct),
  data = intrinsic_count_version,
  family = binomial()
)
summary(mod_binom_count_version)
DHARMa::simulateResiduals(mod_binom_count_version, plot = TRUE)
performance::check_overdispersion(mod_binom_count_version)

icc(mod_binom_count_version)



















##Ignore this until you talk to Maddie
intrinsic_count_version$is_one <- as.numeric(
  (intrinsic_count_version$count_1_pup / intrinsic_count_version$total_resights) == 1
)

#Using a two-hurdle approach to model all 1s separately
mod_bernoulli <- glmmTMB(
  is_one ~ ns(AgeYears, df =3) + year_born_fct + age_last_seen + 
    (1 | animalID_fct) + (1 | season_fct),
  data = intrinsic_count_version,
  family = binomial()
)
summary(mod_bernoulli)
DHARMa::simulateResiduals(mod_bernoulli, plot = TRUE)

mod_binom_non1 <- glmmTMB(
  cbind(count_1_pup, total_resights - count_1_pup) ~ 
    ns(AgeYears, df = 3) + year_born_fct + age_last_seen +
    (1 + AgeYears | animalID_fct) + 
    (1 | season_fct),
  data = intrinsic_count_version[intrinsic_count_version$count_1_pup < intrinsic_count_version$total_resights, ],
  family = binomial()
)
summary(mod_binom_non1)
DHARMa::simulateResiduals(mod_binom_non1, plot = TRUE)

##Plotting the effect of age on association
library(ggeffects)
library(ggplot2)

pred_age_bernoulli <- ggpredict(mod_bernoulli, terms = "AgeYears [all]")
ggplot(pred_age_bernoulli, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(
    title = "Predicted Probability of 1 Pup vs. Age",
    x = "Age (Years)",
    y = "Predicted Probability"
  ) +
  theme_minimal()

pred_age_binom_non1 <- ggpredict(mod_binom_non1, terms = "AgeYears [all]")
ggplot(pred_age_binom_non1, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(
    title = "Daily Probability of 1 Pup Association by Age (Excludes Always-1 Mothers)",
    x = "Age (Years)",
    y = "Predicted Probability"
  ) +
  theme_minimal()

##CONSISTENCY METRICS through this model
##Re-creating the binomial model for consistency

library(performance)

filtered_intrinsic_count_data <- intrinsic_count_version %>%
  group_by(animalID) %>%
  filter(n() >= 3) %>%  
  ungroup()

mod_binom_non1_filtered <- glmmTMB(
  cbind(count_1_pup, total_resights - count_1_pup) ~ 
    ns(AgeYears, df = 3) + year_born_fct + age_last_seen +
    (1 | animalID_fct) + (1 | season_fct),
  data = filtered_intrinsic_count_data[filtered_intrinsic_count_data$count_1_pup < filtered_intrinsic_count_data$total_resights, ],
  family = binomial()
)
summary(mod_binom_non1_filtered)
DHARMa::simulateResiduals(mod_binom_non1_filtered, plot = TRUE)

##Now calculate ICC, which will be more robust
icc(mod_binom_non1_filtered)

