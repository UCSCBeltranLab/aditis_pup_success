##run the source code from Data Processing
source("./DataCurationCode/Data_processing_MPA.R")

#### These are the main models:

##The first model includes all predictors with a smaller subset of data (2016-2023)
##The second model allows us to test for senescence/non-linearity of the age predictor using the full dataset (1996-2025)
##Interaction models with either age*extreme events or age*density using the subset (2016-2023)

########## 1a) 2016-2023 subset model with all predictors ###############

#make dataframe for modeling and figures
intrinsic_2016_2023 <- intrinsic_variables %>%
  filter(season <= 2023, season >= 2016) %>%
  filter(!is.na(avg_density)) %>%
  mutate(season_fct = droplevels(season_fct))

## model for 2016-2023 subset to test all predictors together
## age as a linear predictor (best supported, see 2b)
## random effects of animalID and year

# 2) all predictors model
mod_binom_2016_2023 <- glmer(proportion ~ AgeYears + n_extreme_both + avg_density + (1 | animalID_fct) + (1 | season_fct),
                             weights = total_resights,
                             family = binomial(link= "logit"),
                             control = glmerControl(optimizer = "bobyqa"),
                             data = intrinsic_2016_2023); summary(mod_binom_2016_2023)
ranef(mod_binom_2016_2023)
exp(fixef(mod_binom_2016_2023)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_binom_2016_2023, plot = TRUE) #plot residuals
check_collinearity(mod_binom_2016_2023) #check predictor VIFs

############ 1b) 2016-2023 model with pupping experience instead of age ###################

## experience prior = number of prior breeding seasons observed at least once with a pup
## experience is a weaker predictor, AIC difference from age model = +8
## random effects of animalID and year

mod_binom_exp_2016_2023 <- glmer(proportion ~ experience_prior + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                 weights = total_resights,
                                 family = binomial(link= "logit"),
                                 control = glmerControl(optimizer = "bobyqa"),
                                 data = intrinsic_2016_2023); summary(mod_binom_exp_2016_2023) #model summary
ranef(mod_binom_exp_2016_2023)
exp(fixef(mod_binom_exp_2016_2023)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_binom_exp_2016_2023, plot = TRUE) #plot residuals
check_collinearity(mod_binom_exp_2016_2023) #check predictor VIFs

######### Linear experience figure ###########

# 1) Define x-axis experience values and seasons used for season-specific curves
exp_vals_2016_2023 <- sort(unique(intrinsic_2016_2023$experience_prior)) #unique experience values present (x values)
seasons_2016_2023 <- levels(intrinsic_2016_2023$season_fct) #season factor levels used in predictions

# 2) Compute pooled observed proportion by experience (black points + binomial CIs)
observed_data_exp_2016_2023 <- intrinsic_2016_2023 %>%
  group_by(experience_prior) %>% #group by rows within each experience level
  summarise(n_experience = n(), #sample size at each experience value for labels
            n_success = sum(count_1_pup, na.rm = TRUE), #total successes
            n_trials = sum(total_resights, na.rm = TRUE), #total effort
            avg_prop = n_success/n_trials, #observed proportion
            lwr = binom.test(n_success, n_trials)$conf.int[1], #binomial CI lower
            upr = binom.test(n_success, n_trials)$conf.int[2], #binomial CI upper
            .groups = "drop")

# 3) Build standardized newdata for each experience value for post-stratification
nd_by_experience_2016_2023 <- lapply(exp_vals_2016_2023, \(E) transform(intrinsic_2016_2023, #retain covariate distribution
                                                           experience_prior = E)) #set all seals to experience value E

# 4) Post-stratified prediction curve (row predictions -> weighted mean)
post_curve_exp_2016_2023 <- function(fit, season = NULL){
  vapply(seq_along(nd_by_experience_2016_2023), function(e){
    nd <- nd_by_experience_2016_2023[[e]] #dataset for one experience value
    if(!is.null(season)) nd$season_fct <- factor(season, levels = seasons_2016_2023) #force season level
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL) #predict including all RE
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE) #total resight weighted mean
  }, numeric(1))
}

# 5) Thin season curves (one post-stratified curve per season)
season_lines_exp_2016_2023 <- bind_rows(lapply(seasons_2016_2023, \(S)
                                           tibble(season_fct = S, #season grouping
                                                  experience_prior = exp_vals_2016_2023, #x-axis
                                                  pred = post_curve_exp_2016_2023(mod_binom_exp_2016_2023, season = S)))) #season-specific predictions

# 6) Overall curve + parametric bootstrap CI
set.seed(123) #set for reproducibility
nsim <- 100 #bootstrap draws
post_curve_exp_overall_2016_2023 <- function(fit) post_curve_exp_2016_2023(fit, season = NULL) #overall curve

overall_exp_2016_2023 <- tibble(experience_prior = exp_vals_2016_2023,
                            pred = post_curve_exp_overall_2016_2023(mod_binom_exp_2016_2023)) %>% #mean prediction
  left_join(observed_data_exp_2016_2023 %>% select(experience_prior, n_experience), by = "experience_prior") #add sample sizes

boot_exp_2016_2023 <- bootMer(mod_binom_exp_2016_2023, FUN = post_curve_exp_overall_2016_2023,
                          nsim = nsim, type = "parametric", use.u = FALSE) #parametric bootstrap

# 7) Convert bootstrap draws into 95% CI ribbon
ci_exp_2016_2023 <- overall_exp_2016_2023 %>%
  mutate(lo = apply(boot_exp_2016_2023$t, 2, quantile, 0.025, na.rm = TRUE), #lower CI
         hi = apply(boot_exp_2016_2023$t, 2, quantile, 0.975, na.rm = TRUE)) #upper CI

# 8) Plot: thin season curves + bootstrap ribbon + mean curve + observed data
plot_exp_2016_2023 <- ggplot() +
  geom_line(data = season_lines_exp_2016_2023, #thin season-specific curves
            aes(experience_prior, pred, group = season_fct),
            color = "grey40",
            alpha = 0.3) +
  geom_ribbon(data = ci_exp_2016_2023, #bootstrap CI band
              aes(experience_prior, ymin = lo, ymax = hi),
              fill = "#2BB295", alpha = 0.28) +
  geom_line(data = ci_exp_2016_2023, #overall predicted curve
            aes(experience_prior, pred),
            color = "#2BB295", linewidth = 1.5) +
  geom_pointrange(data = observed_data_exp_2016_2023, #observed pooled proportions
                  aes(experience_prior, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_text(data = observed_data_exp_2016_2023, #sample size labels
            aes(experience_prior, 1.01, label = n_experience),
            vjust = -0.5) +
  scale_x_continuous(breaks = exp_vals_2016_2023) +
  coord_cartesian(ylim = c(0.65, 1.02), clip = "off") +
  theme_few() +
  labs(x = "Previous pupping experience", y = "Mother-offspring association"); plot_exp_2016_2023

plot_exp_jitter_2016_2023 <- ggplot() +
  geom_jitter(data = intrinsic_2016_2023, #raw data weighted by total resights
              aes(experience_prior, proportion, size = total_resights),
              color = "#2BB295",
              alpha = 0.3,
              width = 0.15,
              height = 0.02) +
  geom_line(data = season_lines_exp_2016_2023, #thin season-specific curves
            aes(experience_prior, pred, group = season_fct),
            color = "grey40",
            alpha = 0.4) +
  geom_ribbon(data = ci_exp_2016_2023, #bootstrap CI band
              aes(experience_prior, ymin = lo, ymax = hi),
              fill = "#2BB295", alpha = 0.28) +
  geom_line(data = ci_exp_2016_2023, #overall predicted curve
            aes(experience_prior, pred),
            color = "#2BB295", linewidth = 1.5) +
  geom_text(data = observed_data_exp_2016_2023, #sample size labels
            aes(experience_prior, 1.03, label = n_experience),
            vjust = -0.5) +
  scale_x_continuous(breaks = exp_vals_2016_2023) +
  coord_cartesian(ylim = c(0, 1.03), clip = "off") +
  theme_few() +
  labs(x = "Experience Prior", y = "Mother-offspring association"); plot_exp_jitter_2016_2023

########## Linear age figure ##########

# 1) Define x-axis ages used in model
ages_2016_2023 <- sort(unique(intrinsic_2016_2023$AgeYears)) #unique ages present (x values)

# 2) Compute pooled observed proportion by age (black points + binomial CIs)
observed_data_age_2016_2023 <- intrinsic_2016_2023 %>%
  group_by(AgeYears) %>% #group by rows within each age
  summarise(n_age = n(), #sample size at each age for labels
            n_success = sum(count_1_pup, na.rm = TRUE), #total successes
            n_trials = sum(total_resights, na.rm = TRUE), #total effort
            avg_prop = n_success/n_trials, #observed proportion
            lwr = binom.test(n_success, n_trials)$conf.int[1], #binomial CI lower
            upr = binom.test(n_success, n_trials)$conf.int[2], #binomial CI upper
            .groups = "drop")

# 3) Build standardized newdata for each age for post-stratification
nd_by_age_2016_2023 <- lapply(ages_2016_2023, \(A) transform(intrinsic_2016_2023, #retain covariate distribution
                                         AgeYears = A)) #set all seals to age A

# 4) Post-stratified prediction curve (row predictions -> weighted mean)
post_curve_age_2016_2023 <- function(fit, season = NULL){
  vapply(seq_along(nd_by_age_2016_2023), function(a){
    nd <- nd_by_age_2016_2023[[a]] #dataset for one age
    if(!is.null(season)) nd$season_fct <- factor(season, levels = seasons_2016_2023) #force season level
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL) #predict including all RE
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE) #total resight weighted mean
  }, numeric(1))
}

# 5) Thin season curves (one post-stratified curve per season)
season_lines_age_2016_2023 <- bind_rows(lapply(seasons_2016_2023, \(S)
                                           tibble(season_fct = S, #season grouping
                                                  AgeYears = ages, #x-axis
                                                  pred = post_curve_age_2016_2023(mod_binom_2016_2023, season = S)))) #season-specific predictions

# 6) Overall curve + parametric bootstrap CI
set.seed(123) #set for reproducibility
nsim <- 100 #bootstrap draws
post_curve_age_overall_2016_2023 <- function(fit) post_curve_age_2016_2023(fit, season = NULL) #overall curve

overall_age_2016_2023 <- tibble(AgeYears = ages_2016_2023,
                            pred = post_curve_age_overall_2016_2023(mod_binom_2016_2023)) %>% #mean prediction
  left_join(observed_data_age_2016_2023 %>% select(AgeYears, n_age), by = "AgeYears") #add sample sizes

boot_age_016_2023 <- bootMer(mod_binom_2016_2023, FUN = post_curve_age_overall_2016_2023,
                          nsim = nsim, type = "parametric", use.u = FALSE) #parametric bootstrap

# 7) Convert bootstrap draws into 95% CI ribbon
ci_age_2016_2023 <- overall_age_2016_2023 %>%
  mutate(lo = apply(boot_2016_2023$t, 2, quantile, 0.025, na.rm = TRUE), #lower CI
         hi = apply(boot_2016_2023$t, 2, quantile, 0.975, na.rm = TRUE)) #upper CI

# 8) Plot: thin season curves + bootstrap ribbon + mean curve + observed data
plot_age_2016_2023 <- ggplot() +
  geom_line(data = season_lines_age_2016_2023, #thin season-specific curves
            aes(AgeYears, pred, group = season_fct),
            color = "grey40",
            alpha = 0.4) +
  geom_ribbon(data = ci_age_2016_2023, #bootstrap CI band
              aes(AgeYears, ymin = lo, ymax = hi),
              fill = "#D295E3", alpha = 0.28) +
  geom_line(data = ci_age_2016_2023, #overall predicted curve
            aes(AgeYears, pred),
            color = "#D295E3", linewidth = 1.5) +
  geom_pointrange(data = observed_data_age_2016_2023, #observed pooled proportions
                  aes(AgeYears, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_text(data = observed_data_age_2016_2023, #sample size labels
            aes(AgeYears, 1.01, label = n_age),
            vjust = -0.5) +
  scale_x_continuous(breaks = ages) +
  coord_cartesian(ylim = c(0.75, 1.02), clip = "off") +
  theme_few() +
  labs(x = "Age (Years)", y = "Mother-offspring association"); plot_age_2016_2023

plot_age_jitter_2016_2023 <- ggplot() +
  geom_jitter(data = intrinsic_2016_2023, #raw data weighted by total resights
              aes(AgeYears, proportion, size = total_resights),
              color = "#D295E3",
              alpha = 0.2,
              width = 0.15,
              height = 0.02) +
  geom_line(data = season_lines_age_2016_2023, #thin season-specific curves
            aes(AgeYears, pred, group = season_fct),
            color = "grey40",
            alpha = 0.4) +
  geom_ribbon(data = ci_age_2016_2023, #bootstrap CI band
              aes(AgeYears, ymin = lo, ymax = hi),
              fill = "#D295E3", alpha = 0.4) +
  geom_line(data = ci_age_2016_2023, #overall predicted curve
            aes(AgeYears, pred),
            color = "#D295E3", linewidth = 1.5) +
  geom_text(data = observed_data_age_2016_2023, #sample size labels
            aes(AgeYears, 1.03, label = n_age),
            vjust = -0.5) +
  scale_x_continuous(breaks = ages) +
  coord_cartesian(ylim = c(0, 1.02), clip = "off") +
  theme_few() +
  labs(x = "Age (Years)", y = "Mother-offspring association"); plot_age_jitter_2016_2023

################### Density figure (raw) #####################

# 1) Define x-axis density values across the observed range
density_vals <- seq(min(intrinsic_2016_2023$avg_density, na.rm = TRUE), #minimum observed density
                    max(intrinsic_2016_2023$avg_density, na.rm = TRUE), #maximum observed density
                    length.out = 100) #number of x-axis values for smooth curve

# 2) Build standardized newdata for each density value D (for post-stratification)
# For each D, keep the full predictors of intrinsic_2016_2023 but overwrite avg_density as if everyone had density D
nd_by_density <- lapply(density_vals, \(D) transform(intrinsic_2016_2023, #preserve covariate distribution and weights
                                                     avg_density = D)) #overwrite avg_density so every row is evaluated at density D

# 3) Post-stratified curve: predict for each density, then compute weighted average using total_resights
post_curve_density <- function(fit) {
  vapply(seq_along(nd_by_density), function(d){ #loop over density values
    nd <- nd_by_density[[d]] #newdata for one density value
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL) #predict probability including all random effects
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE) #post-stratified mean, weighted by total resights
  }, numeric(1)) #numeric vector of predictions across density values
}

# 4) Bootstrap CI for the overall density curve using parametric bootstrapping of the fitted GLMM
set.seed(123)
nsim <- 100 #number of bootstrap simulations (can increase if needed, but it takes longer)

overall_density_2016_2023 <- tibble(avg_density = density_vals, #continuous density values
                                    pred = post_curve_density(mod_binom_2016_2023)) #overall post-stratified predictions

boot_density_2016_2023 <- bootMer(mod_binom_2016_2023, FUN = post_curve_density, nsim = nsim, #bootstrap curves
                                  type = "parametric", use.u = FALSE) #parametric bootstrap = simulate ranef each time

# 5) Turn bootstrap curves into 95% CI bands
ci_density_2016_2023 <- overall_density_2016_2023 %>%
  mutate(lo = apply(boot_density_2016_2023$t, 2, quantile, 0.025, na.rm = TRUE), #2.5% quantile at each density value
         hi = apply(boot_density_2016_2023$t, 2, quantile, 0.975, na.rm = TRUE)) #97.5% quantile at each density value

intrinsic_2016_2023 <- intrinsic_2016_2023 %>%
  mutate(season_fct = factor(season_fct, levels = c("2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")))

# 6) Make density plot: bootstrap ribbon + weighted mean line + raw points colored by season
plot_density_jitter <- ggplot() +
  geom_jitter(data = intrinsic_2016_2023, #raw data weighted by total resights
              aes(avg_density, proportion, size = total_resights),
              color = "#6D8DD6",
              alpha = 0.3,
              width = 0.15,
              height = 0.02) +
  geom_ribbon(data = ci_density_2016_2023, #bootstrap CI band for overall curve
              aes(avg_density, ymin = lo, ymax = hi),
              fill = "#1E46A4",
              alpha = 0.4) +
  geom_line(data = ci_density_2016_2023, #overall post-stratified prediction line
            aes(avg_density, pred),
            color = "#1E46A4",
            linewidth = 1.5) +
  scale_color_manual(values = pal, name = "Year") +
  coord_cartesian(ylim = c(0, 1.01), clip = "off") +
  theme_few() +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x  = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)) +
  labs(x = "Location density (average number of seals/10m radius)",
       y = "Mother-offspring association"); plot_density_jitter

ggsave("./TablesFigures/plot_raw_density.png", plot_density_jitter, width = 10, height = 8, dpi = 600)

########### Density figure (binned) ###########

# 1) Define density bins used for observed points
intrinsic_2016_2023$density_bin <- cut_number(intrinsic_2016_2023$avg_density,
                                              n = 8, labels = FALSE)

# 2) Compute pooled observed proportion by density bin (black points + binomial CIs)
observed_density_2016_2023 <- intrinsic_2016_2023 %>%
  group_by(density_bin) %>%
  summarise(avg_density = mean(avg_density, na.rm = TRUE),
            n_bin = n(),
            n_success = sum(count_1_pup, na.rm = TRUE),
            n_trials = sum(total_resights, na.rm = TRUE),
            avg_prop = n_success / n_trials,
            lwr = binom.test(n_success, n_trials)$conf.int[1],
            upr = binom.test(n_success, n_trials)$conf.int[2],
            .groups = "drop")

# 8) Make density plot: bootstrap ribbon + weighted mean line + binned observed data by season
plot_density <- ggplot() +
  geom_ribbon(data = ci_density_2016_2023, #bootstrap CI band for overall curve
              aes(avg_density, ymin = lo, ymax = hi),
              fill = "#598E9E",
              alpha = 0.3) +
  geom_line(data = ci_density_2016_2023, #overall post-stratified prediction line
            aes(avg_density, pred),
            color = "#598E9E",
            linewidth = 1.3) +
  geom_pointrange(data = observed_density_2016_2023, #observed binned proportions by season
                  aes(avg_density, avg_prop, ymin = lwr, ymax = upr, color = avg_density),
                  size = 1.2,
                  linewidth = 1.2) +
  geom_text(data = observed_density_2016_2023, #sample size labels
            aes(avg_density, 1, label = n_bin),
            vjust = -0.5) +
  scale_color_viridis_c(option = "mako", begin = 0.1, end = 0.7, direction = -1, name = "Density") +
  coord_cartesian(ylim = c(0.8, 1), clip = "off") + 
  scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(n.breaks = 3) +
  theme_few() +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x  = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)) +
  labs(x = "Location density (average number of seals in a 10m radius)",
       y = "Mother-offspring association"); plot_density

ggsave("./TablesFigures/plot_density.png", plot_density, width = 10, height = 8, dpi = 600)

########### Extreme events figure (binned) ############

# Set plotting colors for seasons
pal <- c("#E57373","#E6A64C","#E5C76B","#7FBF7B","#2E7D32","#5C7EE5","#6DAEDB","#D77CC8")

# 1) Define x-axis extreme event values and seasons used for plotting
extreme_vals <- sort(unique(intrinsic_2016_2023$n_extreme_both)) #unique extreme event counts present in data

# 2) Compute observed proportion by extreme event count and season (colored points + binomial CI whiskers)
observed_extreme_2016_2023 <- intrinsic_2016_2023 %>%
  group_by(n_extreme_both, season_fct) %>%
  summarise(n_obs = n(), #sample size at each extreme event count within season
            n_success = sum(count_1_pup, na.rm = TRUE), #total successes at each extreme event count within season
            n_trials = sum(total_resights, na.rm = TRUE), #total resights at each extreme event count within season
            avg_prop = n_success / n_trials, #pooled observed proportion at each extreme event count within season
            lwr = binom.test(n_success, n_trials)$conf.int[1], #CI lower
            upr = binom.test(n_success, n_trials)$conf.int[2], #CI upper
            .groups = "drop") 

# make sure years are in order
observed_extreme_2016_2023 <- observed_extreme_2016_2023 %>%
  mutate(season_fct = factor(season_fct, levels = c("2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")))

# 3) Build standardized newdata for each extreme event count E (for post-stratification)
# For each E, keep the full covariate mix of intrinsic_2016_2023 but overwrite n_extreme_both as if everyone experienced E events
nd_by_extreme <- lapply(extreme_vals, \(E) transform(intrinsic_2016_2023, #start from full dataset to preserve covariate distribution and weights
                                                     n_extreme_both = E)) #overwrite n_extreme_both so every row is evaluated at event count E

# 4) Post-stratified curve: predict for each extreme event count, then compute weighted average using total_resights
post_curve_extreme <- function(fit) {
  vapply(seq_along(nd_by_extreme), function(e){ #loop over extreme event counts
    nd <- nd_by_extreme[[e]] #newdata for one extreme event count
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL) #predict probability including all random effects
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE) #post-stratified mean, weighted by total resights
  }, numeric(1)) #numeric vector of predictions across extreme-event counts
}

# 5) Bootstrap CI for the overall extreme event curve using parametric bootstrapping of the fitted GLMM
set.seed(123)
nsim <- 100 #number of bootstrap simulations (can increase if needed, but it takes longer)

overall_extreme_2016_2023 <- tibble(n_extreme_both = extreme_vals, #discrete extreme event counts
                                    pred = post_curve_extreme(mod_binom_2016_2023)) #overall post-stratified predictions

boot_extreme_2016_2023 <- bootMer(mod_binom_2016_2023, FUN = post_curve_extreme, nsim = nsim, #bootstrap curves
                                  type = "parametric", use.u = FALSE) #parametric bootstrap = simulate ranefs each time

# 6) Turn bootstrap curves into 95% CI bands
ci_extreme_2016_2023 <- overall_extreme_2016_2023 %>%
  mutate(lo = apply(boot_extreme_2016_2023$t, 2, quantile, 0.025, na.rm = TRUE), #2.5% quantile at each extreme event count
         hi = apply(boot_extreme_2016_2023$t, 2, quantile, 0.975, na.rm = TRUE)) #97.5% quantile at each extreme event count

intrinsic_2016_2023$season_fct <- factor(intrinsic_2016_2023$season_fct,
                                         levels = c("2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023"))

observed_extreme_2016_2023$season_fct <- factor(observed_extreme_2016_2023$season_fct,
                                                levels = c("2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023"))

# 8) Make extreme events plot: observed seasonal points + bootstrap ribbon + weighted mean line + raw points by season
plot_n_extreme <- ggplot() +
  geom_ribbon(data = ci_extreme_2016_2023, #bootstrap CI band for overall curve
              aes(n_extreme_both, ymin = lo, ymax = hi),
              fill = "#2484A7",
              alpha = 0.4) +
  geom_line(data = ci_extreme_2016_2023, #overall post-stratified prediction line
            aes(n_extreme_both, pred),
            color = "#2484A7",
            linewidth = 1.2) +
  geom_pointrange(data = observed_extreme_2016_2023, #observed points with exact binomial CI by season
                  aes(n_extreme_both, avg_prop, ymin = lwr, ymax = upr, color = season_fct),
                  size = 1.2,
                  linewidth = 1.3) +
  geom_text(data = observed_extreme_2016_2023, #sample size labels
            aes(n_extreme_both, 0.98, label = n_obs),
            vjust = -0.5, hjust = 0.4) +
  scale_color_manual(values = pal, name = "Year") +
  scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(n.breaks = 3) +
  coord_cartesian(ylim = c(0.8, 1), clip = "off") +
  theme_few() +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x  = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)) +
  labs(x = "Number of extreme wave and tide events",
       y = "Mother-offspring association"); plot_n_extreme

ggsave("./TablesFigures/plot_extreme_events.png", plot_n_extreme, width = 10, height = 8, dpi = 600)

############# Extreme events figure (raw) #######################

# 9) Make extreme events plot with jittered raw points
plot_n_extreme_jitter <- ggplot() +
  geom_jitter(data = intrinsic_2016_2023, #raw observations weighted by resights
              aes(n_extreme_both, proportion, size = total_resights, color = season_fct),
              alpha = 0.5,
              width = 0.15,
              height = 0.02) +
  geom_ribbon(data = ci_extreme_2016_2023, #bootstrap CI band for overall curve
              aes(n_extreme_both, ymin = lo, ymax = hi),
              fill = "#2484A7",
              alpha = 0.5) +
  geom_line(data = ci_extreme_2016_2023, #overall post-stratified prediction line
            aes(n_extreme_both, pred),
            color = "#2484A7",
            linewidth = 1.5) +
  geom_pointrange(data = observed_extreme_2016_2023, #observed points with exact binomial CI by season
                  aes(n_extreme_both, avg_prop, ymin = lwr, ymax = upr),
                  color = "black",
                  size = 0.5,
                  linewidth = 1.3) +
  scale_color_manual(values = pal, name = "Year") +
  scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(n.breaks = 5) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  theme_few() +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x  = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)) +
  labs(x = "Number of extreme wave and tide events",
       y = "Mother-offspring association"); plot_n_extreme_jitter

ggsave("./TablesFigures/plot_extreme_events_raw.png", plot_n_extreme_jitter, width = 10, height = 8, dpi = 600)

##### Age threshold test #######

# 1) test all possible thresholds
age_cutoff <- 5:14

# 2) For each cutoff:
# a) define Young/Old at a, and center age at a
# b) fit the same model

age_threshold_test <- function(a) {
  
  d <- intrinsic_variables %>%
    mutate(age_cat = factor(if_else(AgeYears >= a, "Old", "Young"),
                            levels = c("Young","Old")), 
           age10 = (AgeYears - a) / 10)
  
  m <- glmer(proportion ~ age_cat : age10 + (1 | animalID_fct) + (1 | season_fct),
             weights = total_resights,
             family = binomial(link="logit"),
             control = glmerControl(optimizer="bobyqa"),
             data = d)
  
  # 3) Return model fit statistic
  tibble(Threshold = a,
         AIC = AIC(m),
         logLik = as.numeric(logLik(m)))
}

# 4) Fit all thresholds automatically
aic_age_threshold_comparison <- map_dfr(age_cutoff, age_threshold_test) %>%
  mutate(delta_AIC = AIC - min(AIC),
         AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(Threshold)

age_threshold_comparison_tbl <- aic_age_threshold_comparison %>%
  transmute(Threshold = Threshold,
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(AIC_weight, 3)) %>%
  arrange("Threshold")

# identify best threshold
best_age_threshold <- age_threshold_comparison_tbl %>%
  slice_min(AIC, n = 1) %>%
  pull(Threshold)

# bold the lowest AIC threshold
age_threshold_comparison_ft <- flextable(age_threshold_comparison_tbl) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(age_threshold_comparison_tbl$Threshold == best_age_threshold), part = "body"); age_threshold_comparison_ft

save_as_docx(age_threshold_comparison_ft, path = "./TablesFigures/AIC_Age_Threshold_Comparison.docx")

############# 2a) 1996-2025 full dataset model for age ###############

## model for full dataset to estimate age effects
## age with piecewise threshold (best supported, see above section)
## age_cat = "Young", "Old" based on age_senesce
## age10 = (AgeYears - age_senesce) / 10) scaled numeric version of age centered at senescence threshold
## random effects of animalID and year

mod_binom_1996_2025 <- glmer(proportion ~ age_cat : age10 + (1 | animalID_fct) + (1 | season_fct),
                             weights = total_resights,
                             family = binomial(link= "logit"),
                             control = glmerControl(optimizer = "bobyqa"),
                             data = intrinsic_variables); summary(mod_binom_1996_2025)
ranef(mod_binom_1996_2025)
icc(mod_binom_1996_2025, by_group = TRUE)
exp(fixef(mod_binom_1996_2025)) #converts fixed-effect log-odds to odds ratios
resid <- simulateResiduals(mod_binom_1996_2025, plot = TRUE) #plot residuals
plotResiduals(resid, form = intrinsic_variables$age_cat) #age residuals

############## Age figure (piecewise) ################

# 1) Set plotting colors for young/old
COL <- c(Young = "#92BAEE", Old = "#EB99D2")

# 2) Define x-axis ages and the set of seasons to draw thin season-specific curves
ages <- sort(unique(intrinsic_variables$AgeYears)) #unique ages present in data
seasons <- levels(intrinsic_variables$season_fct) #factor levels for seasons (used in predictions)

# 3) Compute pooled observed proportion by age (black points + binomial CI whiskers)
observed_data_1996_2025 <- intrinsic_variables %>%
  group_by(AgeYears) %>% 
  summarise(n_age = n(), #sample size at each age
            n_success = sum(count_1_pup, na.rm = TRUE), #total successes at each age
            n_trials = sum(total_resights, na.rm = TRUE), #total trials/effort at each age
            avg_prop = n_success/n_trials, #pooled observed proportion at each age
            lwr = binom.test(n_success, n_trials)$conf.int[1], #binomial CI lower
            upr = binom.test(n_success, n_trials)$conf.int[2], #binomial CI upper
            .groups = "drop") #return ungrouped tibble

# 4) Build standardized newdata for each age A (for post-stratification)
# For each A, keep the full covariate mix of intrinsic_variables but overwrite age variables as if everyone were age A
nd_by_age <- lapply(ages, \(A) transform(intrinsic_variables, #start from full dataset to preserve covariate distribution and weights
                                         AgeYears = A, #overwrite AgeYears so every row is evaluated at age A
                                         age_cat = factor(ifelse(A < age_senesce, "Young", "Old"), levels = c("Young","Old")), 
                                         age10 = (A - age_senesce)/10)) 

# 5) Post-stratified curve: predict for each age, then compute weighted average using total_resights
# season = NULL gives overall curve, season ="" forces season RE for thin curves
post_curve <- function(fit, season = NULL) {
  vapply(seq_along(nd_by_age), function(a){ #loop over ages
    nd <- nd_by_age[[a]] #newdata for one age
    if(!is.null(season)) nd$season_fct <- factor(season, levels = seasons) #force all rows to a single season level
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL) #predict probability including all random effects
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE) #post-stratified mean, weighted by total resights
  }, numeric(1)) #numeric vector of predictions across ages
}

# 6) Thin season lines: for each season, compute the post-stratified curve
season_lines_1996_2025 <- bind_rows(lapply(seasons, \(S)
                                           tibble(season_fct = S, #season for grouping/plotting
                                                  AgeYears = ages, #x-axis ages
                                                  pred = post_curve(mod_binom_1996_2025, season = S)))) #predicted curve for each season

# 7) Bootstrap CI for the overall curve using parametric bootstrapping of the fitted GLMM
set.seed(123)
nsim <- 100 #number of bootstrap simulations (can increase if needed, but it takes longer)
post_curve_overall <- function(fit) post_curve(fit, season = NULL) #overall curve only

overall_1996_2025 <- tibble(AgeYears = ages, #continuous age in years
                            pred = post_curve_overall(mod_binom_1996_2025)) %>% #overall post-stratified predictions
  left_join(observed_data_1996_2025 %>% select(AgeYears, n_age), by = "AgeYears") #add labels for sample size per age

boot_1996_2025 <- bootMer(mod_binom_1996_2025, FUN = post_curve_overall, nsim = nsim, #bootstrap curves
                          type = "parametric", use.u = FALSE) #parametric bootstrap simulate ranefs each time

# 8) Turn bootstrap curves into 95% CI bands; define Young/Old segment for coloring
ci_1996_2025 <- overall_1996_2025 %>%
  mutate(lo = apply(boot_1996_2025$t, 2, quantile, 0.025, na.rm = TRUE), #2.5% quantile at each age
         hi = apply(boot_1996_2025$t, 2, quantile, 0.975, na.rm = TRUE), #97.5% quantile at each age
         age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"), levels = c("Young","Old")), #segment label for each age
         seg = age_cat) #used to group lines on either side of the threshold in plot

# 9) Add young/old labels to the season-specific lines for correct grouping across the threshold
season_lines_1996_2025 <- season_lines_1996_2025 %>%
  mutate(age_cat = factor(ifelse(AgeYears < age_senesce, "Young", "Old"), levels = c("Young","Old")), #segment label
         seg = age_cat) #keep seg consistent

# 10) Make age plot: thin season curves + bootstrap ribbon + weighted mean line + observed avg points for each age + threshold line + sample size labels
plot_age_1996_2025 <- ggplot() +
  geom_line(data = season_lines_1996_2025, #thin season-specific curves for context
            aes(AgeYears, pred, group = interaction(season_fct, seg)),
            color = "grey40",
            alpha = 0.2) +
  geom_ribbon(data = ci_1996_2025, #bootstrap CI band for overall curve
              aes(AgeYears, ymin = lo, ymax = hi, fill = age_cat, group = seg), #separate ribbons on each side of threshold
              alpha = 0.3) +
  geom_line(data = subset(ci_1996_2025, AgeYears < age_senesce),
            aes(AgeYears, pred, color = age_cat, group = seg),
            linewidth = 1.3) +
  geom_line(data = subset(ci_1996_2025, AgeYears >= age_senesce),
            aes(AgeYears, pred, color = age_cat, group = seg),
            linewidth = 1.3,
            linetype = "dashed") +
  geom_pointrange(data = observed_data_1996_2025, #pooled observed points with binomial CI
                  aes(AgeYears, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  annotate("text", 
           x = age_senesce - 0.5, y = 0.8, 
           label = "Age threshold (9 years)",
           angle = 90, vjust = -0.8, size = 3.5) +
  geom_text(data = observed_data_1996_2025, aes(AgeYears, 1.01, label = n_age), vjust = -0.5) + #sample size labels per age
  coord_cartesian(ylim = c(0.75, 1.02), clip = "off") +
  scale_color_manual(name = "Age class", values = COL) +
  scale_fill_manual(name = "Age class", values = COL) +
  scale_x_continuous(breaks = ages) +
  theme_few() +
  labs(x = "Age (Years)", y = "Mother-offspring association"); plot_age_1996_2025

plot_age_jitter_1996_2025 <- ggplot() +
  geom_jitter(data = intrinsic_variables, #raw data points by age segment
              aes(AgeYears, proportion, color = age_cat, size = total_resights),
              width = 0.18,
              height = 0.02,
              alpha = 0.1,
              show.legend = FALSE) +
  geom_line(data = season_lines_1996_2025, #thin season-specific curves for context
            aes(AgeYears, pred, group = interaction(season_fct, seg)),
            color = "grey40",
            alpha = 0.3) +
  geom_ribbon(data = ci_1996_2025, #bootstrap CI band for overall curve
              aes(AgeYears, ymin = lo, ymax = hi, fill = age_cat, group = seg), #separate ribbons on each side of threshold
              alpha = 0.5) +
  geom_line(data = subset(ci_1996_2025, AgeYears < age_senesce),
            aes(AgeYears, pred, color = age_cat, group = seg),
            linewidth = 1.6) +
  geom_line(data = subset(ci_1996_2025, AgeYears >= age_senesce),
            aes(AgeYears, pred, color = age_cat, group = seg),
            linewidth = 1.6,
            linetype = "dashed") +
  geom_pointrange(data = observed_data_1996_2025, #pooled observed points with binomial CI
                  aes(AgeYears, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  annotate("text", 
           x = age_senesce - 0.5, y = 0.8, 
           label = "Age threshold (9 years)",
           angle = 90, vjust = -0.8, size = 3.5) +
  geom_text(data = observed_data_1996_2025, aes(AgeYears, 1.01, label = n_age), vjust = -0.5) + #sample size labels per age
  coord_cartesian(ylim = c(0, 1.02), clip = "off") +
  scale_color_manual(name = "Age class", values = COL) +
  scale_fill_manual(name = "Age class", values = COL) +
  scale_x_continuous(breaks = ages) +
  theme_few() +
  labs(x = "Age (Years)", y = "Mother-offspring association"); plot_age_jitter_1996_2025

################### Experience threshold test ###################

# 1) test all possible experience thresholds
exp_cutoff <- 1:11

# 2) For each cutoff:
# a) define inexperienced vs. experienced at a, and center experience at a
# b) fit the same model

threshold_test_exp <- function(a) {
  
  d <- intrinsic_variables %>%
    mutate(experience_cat = factor(if_else(experience_prior >= a, "Late", "Early"),
                                   levels = c("Early", "Late")),
           exp10 = (experience_prior - a) / 10)
  
  m <- glmer(proportion ~ experience_cat : exp10 + (1 | animalID_fct) + (1 | season_fct),
             weights = total_resights,
             family = binomial(link = "logit"),
             control = glmerControl(optimizer = "bobyqa"),
             data = d)
  
  tibble(Threshold = a,
         AIC = AIC(m),
         logLik = as.numeric(logLik(m)))
}

# Fit all thresholds automatically
aic_experience_threshold_comparison <- map_dfr(exp_cutoff, threshold_test_exp) %>%
  mutate(delta_AIC = AIC - min(AIC),
         AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(Threshold)

exp_threshold_comparison_tbl <- aic_experience_threshold_comparison %>%
  transmute(Threshold = Threshold,
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(AIC_weight, 3)) %>%
  arrange("Threshold")

# identify best threshold
best_threshold <- aic_experience_threshold_comparison %>%
  slice_min(AIC, n = 1) %>%
  pull(Threshold)

# bold the lowest AIC threshold
exp_threshold_comparison_ft <- flextable(exp_threshold_comparison_tbl) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(exp_threshold_comparison_tbl$Threshold == best_threshold), part = "body"); exp_threshold_comparison_ft

save_as_docx(exp_threshold_comparison_ft, path = "./TablesFigures/AIC_Experience_Threshold_Comparison.docx")

################# 2b) 1996-2025 full dataset model for pupping experience ###############

# 1) Set the threshold for experience_prior to 5
intrinsic_variables <- intrinsic_variables %>%
  mutate(experience_cat = case_when(experience_prior < 5  ~ "inexperienced",
                                    experience_prior >= 5 ~ "experienced"),
         experience_cat = factor(experience_cat, levels = c("inexperienced", "experienced")),
         exp10 = (experience_prior - 5) / 10) %>%
  ungroup()

## model for full dataset using experience instead of age
## piecewise threshold for experience
## experience_cat = "inexperienced" or "experienced", based on <5 or >= 5
## age10 = (experience_prior - 5) / 10) scaled numeric version of experience centered at the threshold
## random effects of animalID and year

# 2) Fit pupping experience model
mod_exp_1996_2025 <- glmer(proportion ~ experience_cat : exp10 + (1 | animalID_fct) + (1 | season_fct),
                           weights = total_resights,
                           family = binomial(link = "logit"),
                           control = glmerControl(optimizer = "bobyqa"),
                           data = intrinsic_variables); summary(mod_exp_1996_2025)

############ Experience figure (piecewise)  ###############

# 1) Colors for inexperienced/experienced
COL <- c(inexperienced = "#9FD46C", experienced = "#B787DA")

# 2) Define x-axis experience and seasons
exp_vals <- sort(unique(intrinsic_variables$experience_prior))
seasons <- levels(intrinsic_variables$season_fct)

# 3) Observed proportions by experience
observed_exp_data <- intrinsic_variables %>%
  group_by(experience_prior) %>% 
  summarise(n_exp = n(),
            n_success = sum(count_1_pup, na.rm = TRUE),
            n_trials = sum(total_resights, na.rm = TRUE),
            avg_prop = n_success/n_trials,
            lwr = binom.test(n_success, n_trials)$conf.int[1],
            upr = binom.test(n_success, n_trials)$conf.int[2],
            .groups = "drop")

# 4) Standardized newdata for each experience E
nd_by_exp <- lapply(exp_vals, \(E) transform(intrinsic_variables,
                                             experience_prior = E,
                                             experience_cat = factor(ifelse(E < 5, "inexperienced", "experienced"), levels = c("inexperienced","experienced")),
                                             exp10 = (E - 5)/10))

# 5) Post-stratified curve
post_curve_exp <- function(fit, season = NULL) {
  vapply(seq_along(nd_by_exp), function(i){
    nd <- nd_by_exp[[i]]
    if(!is.null(season)) nd$season_fct <- factor(season, levels = seasons)
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL)
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE)
  }, numeric(1))
}

# 6) Thin season lines
season_lines_exp <- bind_rows(lapply(seasons, \(S)
                                     tibble(season_fct = S,
                                            experience_prior = exp_vals,
                                            pred = post_curve_exp(mod_exp_1996_2025, season = S))))

# 7) Bootstrap CI
set.seed(123)
nsim <- 100
post_curve_exp_overall <- function(fit) post_curve_exp(fit, season = NULL)

overall_exp <- tibble(experience_prior = exp_vals,
                      pred = post_curve_exp_overall(mod_exp_1996_2025)) %>%
  left_join(observed_exp_data %>% select(experience_prior, n_exp), by = "experience_prior")

boot_exp <- bootMer(mod_exp_1996_2025, FUN = post_curve_exp_overall, nsim = nsim,
                    type = "parametric", use.u = FALSE)

# 8) CIs
ci_exp <- overall_exp %>%
  mutate(lo = apply(boot_exp$t, 2, quantile, 0.025, na.rm = TRUE),
         hi = apply(boot_exp$t, 2, quantile, 0.975, na.rm = TRUE),
         experience_cat = factor(ifelse(experience_prior < 5, "inexperienced", "experienced"), levels = c("inexperienced","experienced")),
         seg = experience_cat)

# 9) Label season lines
season_lines_exp <- season_lines_exp %>%
  mutate(experience_cat = factor(ifelse(experience_prior < 5, "inexperienced", "experienced"), levels = c("inexperienced","experienced")),
         seg = experience_cat)

# 10) Plot
plot_exp_1996_2025 <- ggplot() +
  geom_line(data = season_lines_exp,
            aes(experience_prior, pred, group = interaction(season_fct, seg)),
            color = "grey40",
            alpha = 0.2) +
  geom_ribbon(data = ci_exp,
              aes(experience_prior, ymin = lo, ymax = hi, fill = experience_cat, group = seg),
              alpha = 0.4) +
  geom_line(data = subset(ci_exp, experience_prior < 5),
            aes(experience_prior, pred, color = experience_cat, group = seg),
            linewidth = 1.3) +
  geom_line(data = subset(ci_exp, experience_prior >= 5),
            aes(experience_prior, pred, color = experience_cat, group = seg),
            linewidth = 1.3) +
  geom_pointrange(data = observed_exp_data,
                  aes(experience_prior, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_vline(xintercept = 5 - 0.5, linetype = "dashed") +
  annotate("text",
           x = 5 - 0.5,
           y = 0.8, label = "Experience threshold (5)",
           angle = 90, vjust = -0.8, size = 3.5) +
  coord_cartesian(ylim = c(0.7, 1.02), clip = "off") +
  geom_text(data = observed_exp_data,
            aes(experience_prior, 1.01, label = n_exp), vjust = -0.5) +
  scale_color_manual(name = "Experience class", values = COL) +
  scale_fill_manual(name = "Experience class", values = COL) +
  scale_x_continuous(breaks = exp_vals) +
  theme_few() +
  labs(x = "Prior pupping experience", y = "Mother-offspring association"); plot_exp_1996_2025

#Plot with jitter
plot_exp_jitter_1996_2025 <- ggplot() +
  geom_jitter(data = intrinsic_variables, #raw data points by experience segment
              aes(experience_prior, proportion, color = experience_cat, size = total_resights),
              width = 0.18,
              height = 0.02,
              alpha = 0.3,
              show.legend = FALSE) +
  geom_line(data = season_lines_exp,
            aes(experience_prior, pred, group = interaction(season_fct, seg)),
            color = "grey40",
            alpha = 0.3) +
  geom_ribbon(data = ci_exp,
              aes(experience_prior, ymin = lo, ymax = hi, fill = experience_cat, group = seg),
              alpha = 0.3) +
  geom_line(data = subset(ci_exp, experience_prior < 5),
            aes(experience_prior, pred, color = experience_cat, group = seg),
            linewidth = 1.5) +
  geom_line(data = subset(ci_exp, experience_prior >= 5),
            aes(experience_prior, pred, color = experience_cat, group = seg),
            linewidth = 1.5) +
  geom_vline(xintercept = 5 - 0.5, linetype = "dashed") +
  annotate("text",
           x = 5 - 0.5,
           y = 0.8, label = "Experience threshold (5)",
           angle = 90, vjust = -0.8, size = 3.5) +
  coord_cartesian(ylim = c(0, 1.02), clip = "off") +
  geom_text(data = observed_exp_data,
            aes(experience_prior, 1.01, label = n_exp), vjust = -0.5) +
  scale_color_manual(name = "Experience class", values = COL) +
  scale_fill_manual(name = "Experience class", values = COL) +
  scale_x_continuous(breaks = exp_vals) +
  theme_few() +
  labs(x = "Prior pupping experience", y = "Mother-offspring association"); plot_exp_jitter_1996_2025

ggsave("./TablesFigures/experience_figure.png", plot = plot_exp_1996_2025, width = 14, height = 10, dpi = 600)

#################### 3) Interaction models with age, extreme events, density ##########################

intrinsic_variables <- intrinsic_variables %>%
  mutate(n_extreme_both_scaled = n_extreme_both/10) %>%
  mutate(AgeYears_scaled = AgeYears/10) %>%
  mutate(avg_density_scaled = avg_density/10)

### model with interaction between age*n_extreme_both with all predictors ###
mod_int_weather_2016_2023 <- glmer(proportion ~ age_cat : age10 * n_extreme_both + avg_density + (1 | animalID_fct) + (1 | season_fct),
                                   weights = total_resights,
                                   family = binomial(link = "logit"),
                                   control = glmerControl(optimizer = "bobyqa"),
                                   data = intrinsic_variables); summary(mod_int_weather_2016_2023)
ranef(mod_int_weather_2016_2023)
exp(fixef(mod_int_weather_2016_2023)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_int_weather_2016_2023, plot = TRUE) #plot residuals
check_collinearity(mod_int_weather_2016_2023) #check predictor VIFs

library(ggeffects)
library(ggplot2)

# predicted values across age10 for 3 storm levels
pred_weather <- ggpredict(
  mod_int_weather_2016_2023,
  terms = c("age10 [-8:15 by=1]", "n_extreme_both [all]"),
  condition = c(avg_density = mean(intrinsic_variables$avg_density, na.rm = TRUE))
)

ggplot(pred_weather,
       aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.18, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_y_continuous(
    limits = c(0, 1),
    name = "Predicted mother-offspring association"
  ) +
  labs(
    x = "Age centered at threshold (9 years = 0)",
    color = "Extreme wave+tide events",
    fill = "Extreme wave+tide events"
  ) +
  theme_classic(base_size = 14)

# Predictions
pred_mod_int <- ggpredict(mod_int_weather_2016_2023, 
                          terms = c("n_extreme_both [all]", "AgeYears [all]"))

# Convert age groups to numeric for gradient coloring
pred_mod_int$AgeYears_num <- parse_number(as.character(pred_mod_int$group))

ggplot(pred_mod_int,
       aes(x = x, y = predicted, group = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = AgeYears_num), alpha = 0.15) +
  geom_line(aes(color = AgeYears_num), linewidth = 1, alpha = 0.9) +
  labs(x = "Number of extreme events",
       y = "Predicted MOA") +
  theme_minimal()

### model with interaction between age*density with all predictors ###
mod_int_density_2016_2023 <- glmer(proportion ~ age_cat : age10 * avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                   weights = total_resights,
                                   family = binomial(link = "logit"),
                                   control = glmerControl(optimizer = "bobyqa"),
                                   data = intrinsic_variables); summary(mod_int_density_2016_2023)
ranef(mod_int_density_2016_2023)
exp(fixef(mod_int_density_2016_2023)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_int_density_2016_2023, plot = TRUE) #plot residuals
check_collinearity(mod_int_density_2016_2023) #check predictor VIFs

pred_int_density <- ggpredict(mod_int_density_2016_2023,
                              terms = c("density_c [all]", "AgeYears [all]"))

ggplot(pred_int_density, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  labs(x = "Average density",
       y = "Predicted MOA",
       color = "AgeYears",
       fill = "AgeYears") +
  theme_minimal()

### model with interactions between age*density*n_extreme_both ###
mod_int_all_2016_2023 <- glmer(proportion ~ scale(AgeYears) * scale(avg_density) * scale(n_extreme_both) + (1 | animalID_fct) + (1 | season_fct),
                               weights = total_resights,
                               family = binomial(link = "logit"),
                               control = glmerControl(optimizer = "bobyqa"),
                               data = intrinsic_variables); summary(mod_int_all_2016_2023)
ranef(mod_int_all_2016_2023)
exp(fixef(mod_int_all_2016_2023)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_int_all_2016_2023, plot = TRUE) #plot residuals
check_collinearity(mod_int_all_2016_2023) #check predictor VIFs

#################### 4a) 1996-2025 age AIC comparison ########################

### threshold ###

mod_binom_1996_2025 <- glmer(proportion ~ age_cat : age10 + (1 | animalID_fct) + (1 | season_fct),
                             weights = total_resights,
                             family = binomial(link= "logit"),
                             control = glmerControl(optimizer = "bobyqa"),
                             data = intrinsic_variables); summary(mod_binom_1996_2025) #model summary

### linear ###

mod_age_linear_1996_2025 <- glmmTMB(proportion ~ AgeYears + (1 | animalID_fct) + (1 | season_fct),
                      weights = total_resights,
                      family = binomial(link = "logit"),
                      data = intrinsic_variables); summary(mod_linear)

### quadratic ###

mod_age_quad_1996_2025 <- glmmTMB(proportion ~ age10 + I(age10^2) + (1 | animalID_fct) + (1 | season_fct),
                    weights = total_resights,
                    family = binomial(link = "logit"),
                    data = intrinsic_variables); summary(mod_quad)

## make table with model type comparisons

mod_age_comparisons_1996_2025 <- list(Linear = mod_age_linear_1996_2025,
                        Threshold = mod_binom_1996_2025,
                        Quadratic = mod_age_quad_1996_2025)

aic_table_age_1996_2025 <- tibble(Model = names(mod_age_comparisons_1996_2025),
                              AIC = sapply(mod_age_comparisons_1996_2025, AIC),
                              K = sapply(mod_age_comparisons_1996_2025, function(m) attr(logLik(m), "df"))) %>%
  mutate(delta_AIC = AIC - min(AIC),
         aic_weights = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(delta_AIC)

aic_table_age_1996_2025 <- aic_table_age_1996_2025 %>%
  transmute(Model = Model, 
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(aic_weights, 3)) %>%
  arrange("ΔAIC")

best_model <- aic_table_age_1996_2025$Model[1]

#bold lowest AIC model
aic_table_age_1996_2025 <- flextable(aic_table_age_1996_2025) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(aic_table_age_1996_2025$Model == best_model),
       part = "body"); aic_table_age_1996_2025

#save final table
save_as_docx(aic_table_age_1996_2025, path = "./TablesFigures/1996_2025_Age_Predictor_Comparison.docx")

########## 4b) 2016-2023 age AIC comparison ###################

### threshold ###

mod_age_binom_2016_2023 <- glmer(proportion ~ age_cat : age10 + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                             weights = total_resights,
                             family = binomial(link= "logit"),
                             control = glmerControl(optimizer = "bobyqa"),
                             data = intrinsic_variables); summary(mod_binom_2016_2023) #model summary

### linear ###

mod_age_linear_2016_2023 <- glmmTMB(proportion ~ AgeYears + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                weights = total_resights,
                                family = binomial(link = "logit"),
                                data = intrinsic_variables); summary(mod_linear_2016_2023)

### quadratic ###

mod_age_quad_2016_2023 <- glmmTMB(proportion ~ age10 + I(age10^2) + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                    weights = total_resights,
                    family = binomial(link = "logit"),
                    data = intrinsic_variables); summary(mod_quad_2016_2023)


## name models to compare
mod_age_comparisons_2016_2023 <- list(Linear = mod_age_linear_2016_2023,
                        Quadratic = mod_age_quad_2016_2023,
                        Threshold = mod_age_binom_2016_2023)

# 2) make table with AIC comparisons
aic_table_age_2016_2023 <- tibble(Model = names(mod_age_comparisons_2016_2023),
                              AIC = sapply(mod_age_comparisons_2016_2023, AIC)) %>%
  mutate(delta_AIC = AIC - min(AIC),
         aic_weights = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(delta_AIC)

aic_table_age_2016_2023 <- aic_table_age_2016_2023 %>%
  transmute(Model = Model,
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(aic_weights, 3)) %>%
  arrange("ΔAIC")

best_model <- aic_table_age_2016_2023$Model[1]

#bold lowest AIC model
aic_table_age_2016_2023 <- flextable(aic_table_age_2016_2023) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(aic_table_age_2016_2023$Model == best_model),
       part = "body"); aic_table_age_2016_2023

#save final table
save_as_docx(aic_table_age_2016_2023, path = "./TablesFigures/2016_2023_Age_Predictor_Comparison.docx")

############# 4c) age threshold significance comparison figure (senescence) #################

#test all possible thresholds
age_cutoff <- 5:14

# For each cutoff:
# 1) define Young/Old at a, and center age at a
# 2) fit the same model
# 3) pull the post-threshold slope term (Old × age10), its SE, and p-value
thr_res_age <- map_dfr(age_cutoff, \(a){
  d <- intrinsic_variables %>%
    mutate(age_cat = factor(if_else(AgeYears >= a, "Old", "Young"),
                            levels = c("Young","Old")),
           age10   = (AgeYears - a) / 10)
  
  m <- glmer(proportion ~ age_cat : age10 + (1 | animalID_fct) + (1 | season_fct),
             weights = total_resights,
             family = binomial(link="logit"),
             control = glmerControl(optimizer="bobyqa"),
             data = d)
  
  #coefficient table
  tab <- as.data.frame(summary(m)$coefficients) %>%
    rownames_to_column("term")
  
  #interaction term = slope of age10 for Old
  term_name <- "age_catOld:age10"
  
  row <- filter(tab, term == term_name)
  
  tibble(age_cutoff = a,
         coef = row$Estimate,
         se = row$`Std. Error`,
         p = 2 * pnorm(abs(row$`z value`), lower.tail = FALSE))
}) %>%
  mutate(sig = p < 0.05, # flag significance
         logp_cap = pmin(-log10(p), 5)) # point size weighted by p

# Plot: coefficient ± 95% CI; bigger points = smaller p; color = significance
age_threshold_comparison_figure <- ggplot(thr_res_age, aes(age_cutoff, coef)) +
  geom_hline(yintercept = 0, linetype = "dashed") +  # no post-threshold slope
  geom_pointrange(aes(ymin = coef - 1.96 * se,
                      ymax = coef + 1.96 * se,
                      color = sig,
                      size = logp_cap)) +
  scale_color_manual(values = c(`TRUE` = "deepskyblue", `FALSE` = "tomato"),
                     labels = c(`TRUE` = "p < 0.05", `FALSE` = "p > 0.05"),
                     name = "Significance") +
  scale_size_continuous(range = c(0.5, 2.2), guide = "none") +
  scale_x_continuous(breaks = cutoff) +
  labs(x = "Age senescence threshold",
       y = "Estimated old slope (log-odds)") +
  theme_classic(); age_threshold_comparison_figure

ggsave("./TablesFigures/age_threshold_comparison_figure.png", age_threshold_comparison_figure, width = 8, height = 6, dpi = 600)

############ 5a) 1996-2025 experience AIC comparison #################

### threshold ###

mod_exp_1996_2025 <- glmer(proportion ~ experience_cat : exp10 + (1 | animalID_fct) + (1 | season_fct),
                             weights = total_resights,
                             family = binomial(link= "logit"),
                             control = glmerControl(optimizer = "bobyqa"),
                             data = intrinsic_variables); summary(mod_exp_1996_2025) #model summary

### linear ###

mod_exp_linear_1996_2025 <- glmmTMB(proportion ~ experience_prior + (1 | animalID_fct) + (1 | season_fct),
                                    weights = total_resights,
                                    family = binomial(link = "logit"),
                                    data = intrinsic_variables); summary(mod_exp_linear_1996_2025)

### quadratic ###

mod_exp_quad_1996_2025 <- glmmTMB(proportion ~ exp10 + I(exp10^2) + (1 | animalID_fct) + (1 | season_fct),
                                  weights = total_resights,
                                  family = binomial(link = "logit"),
                                  data = intrinsic_variables); summary(mod_exp_quad_1996_2025)

## make table with model type comparisons

mod_exp_comparisons_1996_2025 <- list(Linear = mod_exp_linear_1996_2025,
                                      Threshold = mod_exp_1996_2025,
                                      Quadratic = mod_exp_quad_1996_2025)

aic_table_exp_1996_2025 <- tibble(Model = names(mod_exp_comparisons_1996_2025),
                                  AIC = sapply(mod_exp_comparisons_1996_2025, AIC),
                                  K = sapply(mod_exp_comparisons_1996_2025, function(m) attr(logLik(m), "df"))) %>%
  mutate(delta_AIC = AIC - min(AIC),
         aic_weights = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(delta_AIC)

aic_table_exp_1996_2025 <- aic_table_exp_1996_2025 %>%
  transmute(Model = Model, 
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(aic_weights, 3)) %>%
  arrange("ΔAIC")

best_model <- aic_table_exp_1996_2025$Model[1]

#bold lowest AIC model
aic_table_exp_1996_2025 <- flextable(aic_table_exp_1996_2025) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(aic_table_exp_1996_2025$Model == best_model),
       part = "body"); aic_table_exp_1996_2025

#save final table
save_as_docx(aic_table_exp_1996_2025, path = "./TablesFigures/1996_2025_Experience_Predictor_Comparison.docx")

############ 5b) 2016-2023 experience AIC comparison ##################

### threshold ###

mod_exp_2016_2023 <- glmer(proportion ~ experience_cat : exp10 + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                 weights = total_resights,
                                 family = binomial(link= "logit"),
                                 control = glmerControl(optimizer = "bobyqa"),
                                 data = intrinsic_variables); summary(mod_exp_2016_2023) #model summary

### linear ###

mod_exp_linear_2016_2023 <- glmmTMB(proportion ~ experience_prior + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                    weights = total_resights,
                                    family = binomial(link = "logit"),
                                    data = intrinsic_variables); summary(mod_exp_linear_2016_2023)

### quadratic ###

mod_exp_quad_2016_2023 <- glmmTMB(proportion ~ exp10 + I(exp10^2) + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                  weights = total_resights,
                                  family = binomial(link = "logit"),
                                  data = intrinsic_variables); summary(mod_exp_quad_2016_2023)


## name models to compare
mod_exp_comparisons_2016_2023 <- list(Linear = mod_exp_linear_2016_2023,
                                      Quadratic = mod_exp_quad_2016_2023,
                                      Threshold = mod_exp_2016_2023)

# 2) make table with AIC comparisons
aic_table_exp_2016_2023 <- tibble(Model = names(mod_exp_comparisons_2016_2023),
                                  AIC = sapply(mod_exp_comparisons_2016_2023, AIC)) %>%
  mutate(delta_AIC = AIC - min(AIC),
         aic_weights = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(delta_AIC)

aic_table_exp_2016_2023 <- aic_table_exp_2016_2023 %>%
  transmute(Model = Model,
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(aic_weights, 3)) %>%
  arrange("ΔAIC")

best_model <- aic_table_exp_2016_2023$Model[1]

#bold lowest AIC model
aic_table_exp_2016_2023 <- flextable(aic_table_exp_2016_2023) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(aic_table_exp_2016_2023$Model == best_model),
       part = "body"); aic_table_exp_2016_2023

#save final table
save_as_docx(aic_table_exp_2016_2023, path = "./TablesFigures/2016_2023_Experience_Predictor_Comparison.docx")

###################### 5c) experience threshold significance figure (senescence) ############################

#test all possible thresholds
exp_cutoff <- 1:11

# For each cutoff:
# 1) define Young/Old at a, and center age at a
# 2) fit the same model
# 3) pull the post-threshold slope term (Old × age10), its SE, and p-value
thr_res_exp <- map_dfr(exp_cutoff, \(e){
  d <- intrinsic_variables %>%
    mutate(experience_cat = factor(if_else(experience_prior >= e, "Experienced", "Inexperienced"),
                            levels = c("Inexperienced","Experienced")),
           exp10 = (experience_prior - e)/10)
  
  m <- glmer(proportion ~ experience_cat : exp10 + (1 | animalID_fct) + (1 | season_fct),
             weights = total_resights,
             family = binomial(link="logit"),
             control = glmerControl(optimizer="bobyqa"),
             data = d)
  
  #coefficient table
  tab <- as.data.frame(summary(m)$coefficients) %>%
    rownames_to_column("term")
  
  #interaction term = slope of age10 for Old
  term_name <- "experience_catExperienced:exp10"
  
  row <- filter(tab, term == term_name)
  
  tibble(exp_cutoff = e,
         coef = row$Estimate,
         se = row$`Std. Error`,
         p = 2 * pnorm(abs(row$`z value`), lower.tail = FALSE))
}) %>%
  mutate(sig = p < 0.05, # flag significance
         logp_cap = pmin(-log10(p), 5)) # point size weighted by p

# Plot: coefficient ± 95% CI; bigger points = smaller p; color = significance
experience_threshold_comparison_figure <- ggplot(thr_res_exp, aes(exp_cutoff, coef)) +
  geom_hline(yintercept = 0, linetype = "dashed") +  # no post-threshold slope
  geom_pointrange(aes(ymin = coef - 1.96 * se,
                      ymax = coef + 1.96 * se,
                      color = sig,
                      size = logp_cap)) +
  scale_color_manual(values = c(`TRUE` = "deepskyblue", `FALSE` = "tomato"),
                     labels = c(`TRUE` = "p < 0.05", `FALSE` = "p > 0.05"),
                     name = "Significance") +
  scale_size_continuous(range = c(0.5, 2.2), guide = "none") +
  scale_x_continuous(breaks = exp_cutoff) +
  labs(x = "Experience senescence threshold",
       y = "Estimated experienced slope (log-odds)") +
  theme_classic(); experience_threshold_comparison_figure

ggsave("./TablesFigures/experience_threshold_comparison_figure.png", experience_threshold_comparison_figure, width = 8, height = 6, dpi = 600)

####################### 6) Weaning mass analysis ######################

##read the weanling weight csv!
weaner_data <- read.csv("./RawData/weanling_weight_master.csv")

weaner_data <- weaner_data %>%
  select(MomAnimalID, year, Wt) %>%
  rename(animalID = MomAnimalID) %>%
  rename(season = year) %>%
  left_join(intrinsic_variables, by = c("animalID", "season")) %>%
  filter(!is.na(Wt)) %>%
  filter(!is.na(proportion)) %>%
  filter(Wt < 200) #filter superweaners

ks.test(weaner_data$Wt, 
        "pnorm", 
        mean = mean(weaner_data$Wt),
        sd = sd(weaner_data$Wt))

#histogram to test normality
ggplot(weaner_data, aes(x = Wt)) +
  geom_histogram(bins = 100, fill = "#2C7FB8", color = "navy") +
  theme_minimal() +
  labs(x = "Weaner weight (Wt)", y = "Frequency")

#Q-Q plot for raw data
ggplot(weaner_data, aes(sample = Wt)) +
  stat_qq() +
  stat_qq_line()

#1) Model with wean mass as response, age and proportion as predictors
mod_wean_age <- lmer(Wt ~ proportion + AgeYears + (1 | animalID_fct) + (1 | season_fct), 
                     data = weaner_data,
                     REML = FALSE); summary(mod_wean_age)
ranef(mod_wean_age)
simulateResiduals(mod_wean_age, plot = TRUE) #plot residuals
check_collinearity(mod_wean_age) #check collinearity
#normality checks
qqnorm(residuals(mod_wean_age)); qqline(residuals(mod_wean_age))
hist(residuals(mod_wean_age))
plot(fitted(mod_wean_age), residuals(mod_wean_age)); abline(h = 0, col = "red")

#2) Predicted effect of proportion
pred_wean_mass <- ggpredict(mod_wean_age, terms = "proportion")

#3) Wean weight mod prediction plot
plot_wean <- ggplot() +
  geom_jitter(data = weaner_data, 
              aes(x = proportion, y = Wt),
              size = 0.3) +
  # confidence band
  geom_ribbon(data = pred_wean_mass,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "#4C78A8",
              alpha = 0.18) +
   # fitted line
  geom_line(data = pred_wean_mass,
            aes(x = x, y = predicted),
            color = "#1F4E79",
            linewidth = 1.3) +
  labs(x = "Mother-offspring association",
       y = "Weaning mass") +
  theme_classic(base_size = 16) +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x  = element_text(size = 16),
        axis.text.y  = element_text(size = 16)); plot_wean

ggsave("./TablesFigures/wean_mass_figure.png", plot_wean, width = 10, height = 8, dpi = 600)
