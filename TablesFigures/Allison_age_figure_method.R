################# trying Allison's figure method (doesn't work because of total resight weighting) ################

library(lme4)

mod_binom_1996_2025 <- glmer(proportion ~ age_cat : age10 + (1 | animalID_fct) + (1 | season_fct),
                             weights = total_resights,
                             family = binomial(link= "logit"),
                             control = glmerControl(optimizer = "bobyqa"),
                             data = intrinsic_variables); summary(mod_binom_1996_2025) #model summary
ranef(mod_binom_1996_2025)
exp(fixef(mod_binom_1996_2025)) #converts fixed-effect log-odds to odds ratios
resid <- simulateResiduals(mod_binom_1996_2025, plot = TRUE) #plot all residuals testDispersion
plotResiduals(resid, form = intrinsic_variables$AgeYears)
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

#Create sample size data frame
n_breed <- intrinsic_variables %>% 
  count(AgeYears)

# Create ANNUAL (i.e. rand effect of year only) predictions
breed_pred_years <- ggpredict(
  mod_binom_1996_2025, 
  terms = c("age10", "age_cat", "season_fct"),
  type = "random") %>% 
  as_tibble() %>% 
  mutate(AgeYears = x * 10 + age_senesce,
         age_cat = factor(group, levels = c("Young", "Old")),
         season_fct = factor(facet)) %>% 
  filter((age_cat == "Young" & AgeYears < age_senesce) |
           (age_cat == "Old" & AgeYears >= age_senesce))

age_rng <- range(intrinsic_variables$AgeYears, na.rm = TRUE)

# Create WEIGHTED-AVERAGE predictions using years' random effects
# ggpredict can't do that for confidence intervals
# create prediction grid
pred_grid <- expand_grid(
  AgeYears   = seq(age_rng[1], age_rng[2]),
  season_fct = levels(intrinsic_2016_2023$season_fct),
) %>% 
  mutate(age10 = (AgeYears - age_senesce) / 10,
         age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old"))) %>%
  mutate(avg_density = mean(intrinsic_variables$avg_density),
         n_extreme_both = mean(intrinsic_variables$n_extreme_both))

seal_years <- count(intrinsic_variables, season_fct)

pred_wgt <- pred_grid %>% 
  mutate(predicted = predict(mod_binom_1996_2025, 
                             newdata = pred_grid,
                             re.form = ~ (1 | season_fct)),
         predicted_se = suppressWarnings(
           # warning: "se.fit computation uses an approximation to estimate the sampling distribution of the parameters"
           predict(mod_binom_1996_2025, 
                   newdata = pred_grid,
                   re.form = NA,
                   se.fit = TRUE)
         )$se.fit,
         predicted_pop = predict(mod_binom_1996_2025, 
                                 newdata = pred_grid,
                                 re.form = NA),
         conf_lo = predicted - 1.96 * predicted_se,
         conf_hi = predicted + 1.96 * predicted_se) %>% 
  left_join(seal_years, by = "season_fct") %>% 
  group_by(AgeYears, age_cat) %>% 
  summarize(across(c(predicted, predicted_pop, conf_lo, conf_hi), 
                   \(x) weighted.mean(x, n)),
            .groups = "drop") %>%
  # invert logit
  mutate(across(c(predicted, predicted_pop, conf_lo, conf_hi), 
                family(mod_binom_1996_2025)$linkinv)) %>% 
  left_join(n_breed, by = "AgeYears") 

pred_wgt %>% filter(AgeYears > 6, 
                    AgeYears <14) %>% 
  summarize(mean_prop = mean(predicted))

pred_wgt %>% 
  mutate(perc_pop = n/sum(n), 
         pred_pups = n * predicted)

breed_plot <- ggplot(pred_wgt, aes(AgeYears, predicted)) +
  # Predictions for individual years
  geom_line(aes(group = interaction(age_cat, season_fct)),
            breed_pred_years,
            alpha = 0.1) +
  # Weighted model (CI and mean response)
  geom_ribbon(aes(fill = age_cat, ymin = conf_lo, ymax = conf_hi), 
              alpha = 0.5) +
  geom_line(aes(color = age_cat),
            linewidth = 1.2) +
  # Population-level mean response
  geom_line(aes(y = predicted_pop, group = age_cat),
            linetype = "dotted",
            linewidth = 1.2,
            alpha = 0.5) +
  # Summarized raw data
  geom_pointrange(aes(y = avg_prop,
                      ymin = lwr,
                      ymax = upr),
                  observed_data) +
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_y_continuous("proportion", 
                     labels = scales::percent,
                     breaks = seq(0.4, 1.0, by = 0.1)) +
  scale_color_manual(values = c("#7fbc41", "#de77ae")) +
  scale_fill_manual(values = c("#7fbc41", "#de77ae")) +
  coord_cartesian(ylim = c(0.8, 1.02)) +
  labs(x = "Age (years)") +
  theme(legend.position = "none") + 
  geom_text(aes(label = n), y = 1, vjust = -0.5)

breed_plot