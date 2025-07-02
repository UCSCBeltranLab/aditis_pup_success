##Modeling intrinsic variables!
library(gamlss)
library(gamlss.dist)

#Linear model- this won't work as proportion is skewed, not normally distributed
mod_Age_by_Proportion <- lm(proportion ~ AgeYears, data = metadata)
summary(mod_Age_by_Proportion)

#The qqplot to prove it in case u don't believe me
qqnorm(metadata$proportion, pch = 1, frame = FALSE)
qqline(metadata$proportion, col = "steelblue", lwd = 2)

##Poisson model for age and proportion
mod_Poisson_Age_by_Proportion <- glm(proportion ~ AgeYears, data = metadata, family = "poisson")
summary(mod_Poisson_Age_by_Proportion)

library(MASS)
##Negative Binomial for age and proportion
mod_NB_Age_by_Proportion <- glm.nb(proportion ~ AgeYears, data = metadata)
summary(mod_NB_Age_by_Proportion)

##Zero-One Inflated Beta Regression model for intrinsic variables

##First model with age and year born as independent effects
mod_01_BR_intrinsic <- gamlss(proportion ~ pb(AgeYears), family = "BEINF", data = intrinsic_variables) #the model
summary(mod_01_BR_intrinsic) #summary statistics
exp(coef(mod_01_BR_Age_YearBorn)) #test for covariance
wp(mod_01_BR_intrinsic, ylim.all = 4) #worm plot
wp(mod_01_BR_intrinsic, xvar = intrinsic_variables$AgeYears)

#Better model with age and year born, made both smooth/non-linear
mod_02_BR_intrinsic_smooth <- gamlss(proportion ~ pb(AgeYears) + pb(year_born) + age_last_seen + random(animalID_fct) + random(season_fct),
  sigma.formula = ~ pb(AgeYears) + pb(year_born) + age_last_seen,
  tau.formula   = ~ pb(AgeYears) + pb(year_born) + age_last_seen,
  family = BEINF,
  data = intrinsic_variables)
summary(mod_02_BR_intrinsic_smooth) #summary stats
exp(coef(mod_02_BR_intrinsic_smooth)) #test for covariance
wp(mod_02_BR_intrinsic_smooth, ylim.all = 4) #worm plot
wp(mod_02_BR_intrinsic_smooth, xvar = intrinsic_variables$AgeYears, ylim.worm = 4) #worm plot for Age
plot(mod_02_BR_intrinsic_smooth) #model diagnostics

#Visualize tau, probability of getting proportion = 1 with age and year_born
term.plot(mod_02_BR_intrinsic_smooth, what = "tau", ask = FALSE) 

##try this for assessing co-linearity
cor(intrinsic_variables$AgeYears, intrinsic_variables$year_born)
cor(intrinsic_variables$AgeYears, intrinsic_variables$age_last_seen)

#try cohorts to see if it fixes adding cohort to model
intrinsic_variables$cohort <- cut(
  intrinsic_variables$year_born,
  breaks = c(1980, 1990, 2000, 2010, 2015, 2022),
  labels = c("1981–89", "1990–99", "2000–09", "2010–14", "2015–21"),
  right = FALSE)

##Cohorts work for this model
mod_03_BR_intrinsic_smooth <- gamlss(proportion ~ pb(AgeYears) + as.factor(cohort) + random(as.factor(animalID)) + random(as.factor(season)),
  sigma.formula = ~ pb(AgeYears)
  tau.formula   = ~ pb(AgeYears)
  family = BEINF,
  data = intrinsic_variables,
  trace = TRUE)
summary(mod_03_BR_intrinsic_smooth) #summary stats
exp(coef(mod_03_BR_intrinsic_smooth)) #test for covariance
wp(mod_03_BR_intrinsic_smooth, ylim.all = 4) #worm plot
wp(mod_03_BR_intrinsic_smooth, xvar = intrinsic_variables$AgeYears, ylim.worm = 2) #worm plot for Age
plot(mod_03_BR_intrinsic_smooth)


##I HATE THIS DO NOT USE BRMS
#Bayesian regression models using Stan
##brms version of beta 1-inflated model
library(brms)

intrinsic_variables$is_one <- as.numeric(intrinsic_variables$proportion == 1)

mod_brms_BEINF <- brm(
  formula = bf(
    proportion ~ s(AgeYears) + s(year_born) + age_last_seen + (1 | animalID_fct) + (1 | season_fct),
    phi ~ s(AgeYears) + s(year_born) + age_last_seen, #like sigma
    zoi ~ s(AgeYears) + s(year_born) + age_last_seen), #like tau
  family = zero_one_inflated_beta(),
  data = intrinsic_variables,
  chains = 4,
  cores = 4,
  iter = 4000,
  control = list(adapt_delta = 0.95))

summary(mod_brms_BEINF)

#term plots for each effect
plot(conditional_effects(mod_brms_BEINF, "AgeYears")) ##mu
plot(conditional_effects(mod_brms_BEINF, "AgeYears", dpar = "zoi")) #1-inflation
plot(conditional_effects(mod_brms_BEINF, "year_born"))

pp_check(mod_brms_BEINF)

res <- residuals(mod_brms_BEINF, summary = TRUE)
fitted_vals <- fitted(mod_brms_BEINF)[, "Estimate"]

plot(fitted_vals, res[, "Estimate"],
     xlab = "Fitted", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2)

pp_check(mod_brms_BEINF, type = "dens_overlay")    # Distribution match
pp_check(mod_brms_BEINF, type = "hist")            # Histogram comparison
pp_check(mod_brms_BEINF, type = "scatter_avg")     # Fitted vs observed
pp_check(mod_brms_BEINF, type = "mean")            # Predicted vs observed means

##Run two-part model with only one and beta components
mod_bernoulli_brms <- brm(
  is_one ~ s(AgeYears) + s(year_born) + age_last_seen + (1 | animalID_fct) + (1 | season_fct),
  family = bernoulli(link = "logit"),
  data = intrinsic_variables,
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)
summary(mod_bernoulli_brms)
pp_check(mod_bernoulli_brms, type = "bars")
pp_check(mod_bernoulli_brms, type = "scatter_avg")
bayesplot::ppc_bars(y = intrinsic_variables$is_one, yrep = yrep_small)

plot(conditional_effects(mod_bernoulli_brms, "AgeYears"))

mod_beta_brms <- brm(
  proportion ~ s(AgeYears) + s(year_born) + age_last_seen + (1 | animalID_fct) + (1 | season_fct),
  family = beta(),
  data = intrinsic_variables[intrinsic_variables$proportion < 1, ],
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)
summary(mod_beta_brms)
