
######################### Extrinsic model with density and extreme events (2016-2023) ##########################

##binomial
mod_extrinsic_1a <- glmmTMB(is_one ~ avg_density + n_extreme_both + total_resights + (1 | season_fct) + (1 | animalID_fct) + (1 | area_fct),
                               data = extrinsic_variables,
                               family = binomial(link = "logit"))
summary(mod_extrinsic_1a)
exp(fixef(mod_extrinsic_1a)$cond) 
DHARMa::simulateResiduals(mod_extrinsic_1a, plot = TRUE)
#residuals look fine
check_collinearity(mod_extrinsic_1a)

##gamma
mod_extrinsic_1b <- glmmTMB(flipped_prop ~ avg_density + n_extreme_both + total_resights + (1 | season_fct) + (1 | animalID_fct) + (1| area_fct), 
                            data = extrinsic_variables_sub,
                            family = Gamma(link = "log"))
summary(mod_extrinsic_1b)
exp(fixef(mod_extrinsic_1b)$cond) 
DHARMa::simulateResiduals(mod_extrinsic_1b, plot = TRUE)
check_collinearity(mod_extrinsic_1b)



#try high, medium, low density areas or NP/SP 
mod_extrinsic_3 <- glmmTMB(is_one ~ n_extreme_both + total_resights + (1 | season_fct),
                           data = extrinsic_weather,
                           family = binomial(link = "logit"))
summary(mod_extrinsic_3)
exp(fixef(mod_extrinsic_3)$cond) 
DHARMa::simulateResiduals(mod_extrinsic_3, plot = TRUE)
#residuals look fine
check_collinearity(mod_extrinsic_3)

#################################### Figure development! ####################################

range(extrinsic_variables$avg_density) 

pred_mod_extrinsic_1 <- ggpredict(mod_extrinsic_1, terms = c("AgeYears [all]", "avg_density [5,20]"))

ggplot(pred_mod_extrinsic_1, aes(x = x, y = predicted, colour = group, fill = group)) +
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

