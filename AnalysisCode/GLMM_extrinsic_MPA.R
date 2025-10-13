
######################### Extrinsic model with density and extreme events (2016-2023) ##########################
##Setting threshold in data
extrinsic_variables <- extrinsic_variables %>%
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"),
                          levels = c("Young", "Old"))) %>%
  mutate(age10 = (AgeYears - age_senesce) / 10) #scaled numeric version of Age centered at senescence threshold

##binomial
mod_extrinsic_1a <- glmmTMB(is_one ~ AgeYears * avg_density + n_extreme_tide + total_resights + (1 | season_fct) + (1 | animalID_fct),
                            data = extrinsic_variables,
                            family = binomial(link = "logit"))
summary(mod_extrinsic_1a)
exp(fixef(mod_extrinsic_1a)$cond) 
DHARMa::simulateResiduals(mod_extrinsic_1a, plot = TRUE)
#residuals look fine
check_collinearity(mod_extrinsic_1a)

##gamma
mod_extrinsic_1b <- glmmTMB(flipped_prop ~ age_cat : avg_density + n_extreme_both + total_resights + (1 | season_fct) + (1 | animalID_fct) + (1 | area_fct),
                            data = extrinsic_variables_sub,
                            family = Gamma(link = "log"))
summary(mod_extrinsic_1b)
exp(fixef(mod_extrinsic_1b)$cond) 
DHARMa::simulateResiduals(mod_extrinsic_1b, plot = TRUE)
check_collinearity(mod_extrinsic_1b)


#try high, medium, low density areas or NP/SP 
mod_extrinsic_3 <- glmmTMB(is_one ~ avg_density + area_fct + total_resights + (1 | season_fct),
                           data = extrinsic_variables,
                           family = binomial(link = "logit"))
summary(mod_extrinsic_3)
exp(fixef(mod_extrinsic_3)$cond) 
DHARMa::simulateResiduals(mod_extrinsic_3, plot = TRUE)
#residuals look fine
check_collinearity(mod_extrinsic_3)

#################################### Interaction models  ####################################

mod_extrinsic_int <- glmmTMB(is_one ~ avg_density * AgeYears + n_extreme_both * AgeYears + total_resights + (1 | season_fct) + (1 | animalID_fct),
                            data = extrinsic_variables,
                            family = binomial(link = "logit"))
summary(mod_extrinsic_int)
exp(fixef(mod_extrinsic_int)$cond) 
DHARMa::simulateResiduals(mod_extrinsic_int, plot = TRUE)
#residuals look fine
check_collinearity(mod_extrinsic_int)

mod_extrinsic_int2 <- glmmTMB(flipped_prop ~ avg_density * AgeYears + n_extreme_both * AgeYears + total_resights + (1 | season_fct) + (1 | animalID_fct) + (1 | area_fct), 
                            data = extrinsic_variables_sub,
                            family = Gamma(link = "log"))
summary(mod_extrinsic_int2)

exp(fixef(mod_extrinsic_int2)$cond) 
DHARMa::simulateResiduals(mod_extrinsic_int2, plot = TRUE)
check_collinearity(mod_extrinsic_int2)

range(extrinsic_variables$avg_density) 

pred_mod_extrinsic_int <- ggpredict(mod_extrinsic_int, terms = c("AgeYears [all]", "avg_density [5,20]"))

ggplot(pred_mod_extrinsic_int, aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_line(linewidth = 1.3) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.25, colour = NA) +
  scale_color_manual(values = c("lightblue", "navy")) +
  scale_fill_manual(values = c("lightblue", "navy")) +
  labs(x = "Age (years)",
       y = "Probability of perfect proportion (0–1)",
       colour = "Density", fill = "Density") +
  theme_few(base_size = 14) +
  theme(legend.position = "top",
        legend.title = element_blank())

