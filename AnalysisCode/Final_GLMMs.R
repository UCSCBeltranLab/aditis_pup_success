##run the source code from Data Processing
source("./DataCurationCode/Data_processing_MPA.R")

#### These are the main models Maddie and I think make most sense:

##The first model includes all predictors with a smaller subset of data (2016-2023)
##The second model allows us to test for senescence/non-linearity of the age predictor using the full dataset (1996-2025)

##Marm in case you're wondering why I don't just use a single dataframe for all the models - it's just for figure making consistency!

########## 1a) 2016-2023 subset model with all predictors ###############

# 1) Subset the data for only seasons 2016 - 2023
intrinsic_2016_2023 <- intrinsic_variables %>%
  filter(season %in% 2016 : 2023) %>%
  filter(!is.na(avg_density)) %>%
  droplevels()

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

############ 1b) 2016-2023 model with experience instead of age ###################
mod_binom_exp_2016_2023 <- glmer(proportion ~ experience_prior + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                 weights = total_resights,
                                 family = binomial(link= "logit"),
                                 control = glmerControl(optimizer = "bobyqa"),
                                 data = intrinsic_2016_2023); summary(mod_binom_exp_2016_2023) #model summary
ranef(mod_binom_exp_2016_2023)
exp(fixef(mod_binom_exp_2016_2023)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_binom_exp_2016_2023, plot = TRUE) #plot residuals
check_collinearity(mod_binom_exp_2016_2023) #check predictor VIFs

#predicted effect of pupping experience
plot(ggpredict(mod_binom_exp_2016_2023, terms = "experience_prior"))

########## Linear age figure ##########

# 1) Define x-axis ages and seasons used for season-specific curves
ages <- sort(unique(intrinsic_2016_2023$AgeYears)) #unique ages present (x values)
seasons <- levels(intrinsic_2016_2023$season_fct) #season factor levels used in predictions

# 2) Compute pooled observed proportion by age (black points + binomial CIs)
observed_data_2016_2023 <- intrinsic_2016_2023 %>%
  group_by(AgeYears) %>% #group by rows within each age
  summarise(n_age = n(), #sample size at each age for labels
            n_success = sum(count_1_pup, na.rm = TRUE), #total successes
            n_trials = sum(total_resights, na.rm = TRUE), #total effort
            avg_prop = n_success/n_trials, #observed proportion
            lwr = binom.test(n_success, n_trials)$conf.int[1], #binomial CI lower
            upr = binom.test(n_success, n_trials)$conf.int[2], #binomial CI upper
            .groups = "drop")

# 3) Build standardized newdata for each age for post-stratification
nd_by_age <- lapply(ages, \(A) transform(intrinsic_2016_2023, #retain covariate distribution
                                         AgeYears = A)) #set all seals to age A

# 4) Post-stratified prediction curve (row predictions -> weighted mean)
post_curve <- function(fit, season = NULL){
  vapply(seq_along(nd_by_age), function(a){
    nd <- nd_by_age[[a]] #dataset for one age
    if(!is.null(season)) nd$season_fct <- factor(season, levels = seasons) #force season level
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL) #predict including all RE
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE) #total resight weighted mean
  }, numeric(1))
}

# 5) Thin season curves (one post-stratified curve per season)
season_lines_2016_2023 <- bind_rows(lapply(seasons, \(S)
                                           tibble(season_fct = S, #season grouping
                                                  AgeYears = ages, #x-axis
                                                  pred = post_curve(mod_binom_2016_2023, season = S)))) #season-specific predictions

# 6) Overall curve + parametric bootstrap CI
set.seed(123) #set for reproducibility
nsim <- 100 #bootstrap draws
post_curve_overall <- function(fit) post_curve(fit, season = NULL) #overall curve

overall_2016_2023 <- tibble(AgeYears = ages,
                            pred = post_curve_overall(mod_binom_2016_2023)) %>% #mean prediction
  left_join(observed_data_2016_2023 %>% select(AgeYears, n_age), by = "AgeYears") #add sample sizes

boot_2016_2023 <- bootMer(mod_binom_2016_2023, FUN = post_curve_overall,
                          nsim = nsim, type = "parametric", use.u = FALSE) #parametric bootstrap

# 7) Convert bootstrap draws into 95% CI ribbon
ci_2016_2023 <- overall_2016_2023 %>%
  mutate(lo = apply(boot_2016_2023$t, 2, quantile, 0.025, na.rm = TRUE), #lower CI
         hi = apply(boot_2016_2023$t, 2, quantile, 0.975, na.rm = TRUE)) #upper CI

# 8) Plot: thin season curves + bootstrap ribbon + mean curve + observed data
plot_age_2016_2023 <- ggplot() +
  geom_line(data = season_lines_2016_2023, #thin season-specific curves
            aes(AgeYears, pred, group = season_fct, color = season_fct),
            alpha = 0.8) +
  geom_ribbon(data = ci_2016_2023, #bootstrap CI band
              aes(AgeYears, ymin = lo, ymax = hi),
              fill = "#83C05A", alpha = 0.28) +
  geom_line(data = ci_2016_2023, #overall predicted curve
            aes(AgeYears, pred),
            color = "#83C05A", linewidth = 1.5) +
  geom_pointrange(data = observed_data_2016_2023, #observed pooled proportions
                  aes(AgeYears, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_text(data = observed_data_2016_2023, #sample size labels
            aes(AgeYears, 1.01, label = n_age),
            vjust = -0.5) +
  scale_x_continuous(breaks = ages) +
  coord_cartesian(ylim = c(0.75, 1.02), clip = "off") +
  theme_few() +
  labs(x = "Age (Years)", y = "Mother-offspring association"); plot_age_2016_2023

################### Density figure (raw) #####################

# Set plotting colors for seasons
pal <- c("#E57373","#E6A64C","#E5C76B","#7FBF7B","#2E7D32","#5C7EE5","#6DAEDB","#D77CC8")

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

# 6) Make density plot: bootstrap ribbon + weighted mean line + raw points colored by season
plot_density_jitter <- ggplot() +
  geom_ribbon(data = ci_density_2016_2023, #bootstrap CI band for overall curve
              aes(avg_density, ymin = lo, ymax = hi),
              fill = "#916EB4",
              alpha = 0.3) +
  geom_jitter(data = intrinsic_2016_2023, #raw observations colored by season
              aes(avg_density, proportion, color = season_fct),
              alpha = 0.7,
              size = 3) +
  geom_line(data = ci_density_2016_2023, #overall post-stratified prediction line
            aes(avg_density, pred),
            color = "#916EB4",
            linewidth = 1.3) +
  scale_color_manual(values = pal, name = "Year") +
  coord_cartesian(ylim = c(0, 1.01), clip = "off") +
  theme_few() +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x  = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)) +
  labs(x = "Within-colony location density",
       y = "Mother-offspring association"); plot_density_jitter

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
              fill = "#916EB4",
              alpha = 0.3) +
  geom_line(data = ci_density_2016_2023, #overall post-stratified prediction line
            aes(avg_density, pred),
            color = "#916EB4",
            linewidth = 1.3) +
  geom_pointrange(data = observed_density_2016_2023, #observed binned proportions by season
                  aes(avg_density, avg_prop, ymin = lwr, ymax = upr),
                  color = "black",
                  size = 1.2,
                  linewidth = 1.2) +
  geom_text(data = observed_density_2016_2023, #sample size labels
            aes(avg_density, 1, label = n_bin),
            vjust = -0.5) +
  scale_color_manual(values = pal, name = "Year") +
  coord_cartesian(ylim = c(0.8, 1), clip = "off") + 
  scale_y_continuous(n.breaks = 6) +
  theme_few() +
  labs(x = "Within-colony location density",
       y = "Mother-offspring association"); plot_density

ggsave("./TablesFigures/plot_density.png", plot_density, width = 14, height = 12, dpi = 600)

########### Extreme events figure (binned) ############

# 1) Define x-axis extreme event values and seasons used for plotting
extreme_vals <- sort(unique(intrinsic_2016_2023$n_extreme_both)) #unique extreme event counts present in data
seasons <- levels(intrinsic_2016_2023$season_fct) #factor levels for seasons (used for plotting colors)

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


# 8) Make extreme events plot: observed seasonal points + bootstrap ribbon + weighted mean line + raw points by season
plot_n_extreme <- ggplot() +
  geom_pointrange(data = observed_extreme_2016_2023, #observed points with exact binomial CI by season
                  aes(n_extreme_both, avg_prop, ymin = lwr, ymax = upr, color = season_fct),
                  size = 1.2,
                  linewidth = 1.3) +
  geom_ribbon(data = ci_extreme_2016_2023, #bootstrap CI band for overall curve
              aes(n_extreme_both, ymin = lo, ymax = hi),
              fill = "#916EB4",
              alpha = 0.4) +
  geom_line(data = ci_extreme_2016_2023, #overall post-stratified prediction line
            aes(n_extreme_both, pred),
            color = "#916EB4",
            linewidth = 1.2) +
  geom_text(data = observed_extreme_2016_2023, #sample size labels
            aes(n_extreme_both, 0.98, label = n_obs),
            vjust = -0.5, hjust = 0.4) +
  scale_color_manual(values = pal, name = "Year") +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 3) +
  coord_cartesian(ylim = c(0.8, 1), clip = "off") +
  theme_few() +
  labs(x = "Number of extreme wave and tide events",
       y = "Mother-offspring association"); plot_n_extreme

ggsave("./TablesFigures/plot_extreme_events.png", plot_n_extreme, width = 14, height = 12, dpi = 600)

############# Extreme events figure (raw) #######################

# 9) Make extreme events plot with jittered raw points
plot_n_extreme_jitter <- ggplot() +
  geom_pointrange(data = observed_extreme_2016_2023, #observed points with exact binomial CI by season
                  aes(n_extreme_both, avg_prop, ymin = lwr, ymax = upr, color = season_fct),
                  size = 1.2,
                  linewidth = 1.3) +
  geom_jitter(data = intrinsic_2016_2023, #raw observations colored by season
              aes(n_extreme_both, proportion, color = season_fct),
              alpha = 0.7,
              size = 4) +
  geom_ribbon(data = ci_extreme_2016_2023, #bootstrap CI band for overall curve
              aes(n_extreme_both, ymin = lo, ymax = hi),
              fill = "#916EB4",
              alpha = 0.4) +
  geom_line(data = ci_extreme_2016_2023, #overall post-stratified prediction line
            aes(n_extreme_both, pred),
            color = "#916EB4",
            linewidth = 1.2) +
  scale_color_manual(values = pal, name = "Year") +
  scale_x_continuous(n.breaks = 10) +
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

############# 2a) 1996-2025 full dataset model for age ###############

## model for full dataset to estimate age effects
## age with piecewise threshold (best supported, see 2a)
## age_cat = "Young", "Old" based on age_senesce
## age10 = (AgeYears - age_senesce) / 10) scaled numeric version of age centered at senescence threshold
## random effects of animalID and year

mod_binom_1996_2025 <- glmer(proportion ~ age_cat : age10 + (1 | animalID_fct) + (1 | season_fct),
                             weights = total_resights,
                             family = binomial(link= "logit"),
                             control = glmerControl(optimizer = "bobyqa"),
                             data = intrinsic_variables); summary(mod_binom_1996_2025)
ranef(mod_binom_1996_2025)
exp(fixef(mod_binom_1996_2025)) #converts fixed-effect log-odds to odds ratios
resid <- simulateResiduals(mod_binom_1996_2025, plot = TRUE) #plot residuals
plotResiduals(resid, form = intrinsic_variables$age_cat) #age residuals

################# 2b) 1996-2025 full dataset model for pupping experience ###############

### compare to model with experience threshold ###
intrinsic_variables <- intrinsic_variables %>%
  mutate(
    experience_cat = case_when(
      experience_prior < 5  ~ "early",
      experience_prior >= 5 ~ "late"
    ),
    experience_cat = factor(experience_cat, levels = c("early", "late")),
    
    # center at threshold (5 pups) and scale per 10 pups
    exp10 = (experience_prior - 5) / 10
  )

mod_exp_1996_2025 <- glmer(
  proportion ~ experience_cat : exp10 +
    (1 | animalID_fct) + (1 | season_fct),
  weights = total_resights,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa"),
  data = intrinsic_variables
); summary(mod_exp_1996_2025)

# test all possible experience thresholds
cutoff <- 1:14

# For each cutoff:
# 1) define Early/Late at a, and center experience at a
# 2) fit the same model

threshold_test_exp <- function(a) {
  
  d <- intrinsic_variables %>%
    mutate(
      experience_cat = factor(
        if_else(experience_prior >= a, "Late", "Early"),
        levels = c("Early", "Late")
      ),
      exp10 = (experience_prior - a) / 10
    )
  
  m <- glmer(
    proportion ~ experience_cat : exp10 + (1 | animalID_fct) + (1 | season_fct),
    weights = total_resights,
    family = binomial(link = "logit"),
    control = glmerControl(optimizer = "bobyqa"),
    data = d
  )
  
  # Return model fit statistic
  tibble(
    Threshold = a,
    AIC = AIC(m),
    logLik = as.numeric(logLik(m))
  )
}

# Fit all thresholds automatically
aic_experience_threshold_comparison <- map_dfr(cutoff, threshold_test_exp) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) %>%
  arrange(Threshold)

threshold_comparison_tbl <- aic_experience_threshold_comparison %>%
  transmute(
    Threshold = Threshold,
    AIC = round(AIC, 1),
    `ΔAIC` = round(delta_AIC, 2),
    `AIC weight` = round(AIC_weight, 3)
  ) %>%
  arrange(Threshold)

# identify best threshold
best_threshold <- aic_experience_threshold_comparison %>%
  slice_min(AIC, n = 1) %>%
  pull(Threshold)

# bold the lowest AIC threshold
threshold_comparison_ft <- flextable(threshold_comparison_tbl) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(threshold_comparison_tbl$Threshold == best_threshold), part = "body"); threshold_comparison_ft

############## Threshold age figure ################

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
            aes(AgeYears, pred, group = interaction(season_fct, seg), color = age_cat), #separate lines
            alpha = 0.28) +
  geom_ribbon(data = ci_1996_2025, #bootstrap CI band for overall curve
              aes(AgeYears, ymin = lo, ymax = hi, fill = age_cat, group = seg), #separate ribbons on each side of threshold
              alpha = 0.28) +
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
           x = age_senesce - 0.5,
           y = 0.8, label = "Age threshold (9 years)",
           angle = 90, vjust = -0.8, size = 3.5) +
  coord_cartesian(ylim = c(0.75, 1.02), clip = "off") +
  geom_text(data = observed_data_1996_2025, aes(AgeYears, 1.01, label = n_age), vjust = -0.5) + #sample size labels per age
  scale_color_manual(name = "Age class", values = COL) +
  scale_fill_manual(name = "Age class", values = COL) +
  scale_x_continuous(breaks = ages) +
  theme_few() +
  labs(x = "Age (Years)", y = "Mother-offspring association"); plot_age_1996_2025

# 10) Plot with colors for year random effect lines
plot_age_season_1996_2025 <- ggplot() +
  geom_line(data = season_lines_1996_2025, # thin season-specific curves
            aes(AgeYears, pred,
                group = interaction(season_fct, seg),
                color = season_fct),
            alpha = 0.45,
            linewidth = 0.6) +
  geom_ribbon(data = ci_1996_2025, # bootstrap CI band for overall curve
              aes(AgeYears, ymin = lo, ymax = hi,
                  fill = age_cat, group = seg),
              alpha = 0.28) +
  geom_line(data = subset(ci_1996_2025, AgeYears < age_senesce), # overall pre-threshold curve
            aes(AgeYears, pred, group = seg),
            linewidth = 1.3,
            color = "#92BAEE") +
  geom_line(data = subset(ci_1996_2025, AgeYears >= age_senesce), # overall post-threshold curve
            aes(AgeYears, pred, group = seg),
            linewidth = 1.3,
            linetype = "dashed",
            color = "#EB99D2") +
  geom_pointrange(data = observed_data_1996_2025, # pooled observed points
                  aes(AgeYears, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_vline(xintercept = age_senesce, linetype = "dashed") +
  annotate("text",
           x = age_senesce,
           y = 0.8,
           label = "Age threshold (9 years)",
           angle = 90,
           vjust = -0.8,
           size = 3.5) +
  geom_text(data = observed_data_1996_2025, # sample size labels
            aes(AgeYears, 1.01, label = n_age),
            vjust = -0.5) +
  coord_cartesian(ylim = c(0.75, 1.02), clip = "off") +
  scale_color_hue(name = "Season") +
  scale_fill_manual(name = "Age class", values = COL) +
  scale_x_continuous(breaks = ages) +
  theme_few() +
  labs(x = "Age (Years)", y = "Mother-offspring association"); plot_age_season_1996_2025

#################### 3) Interaction models with age, extreme events, density ##########################
intrinsic_variables$Age_c <- (intrinsic_variables$AgeYears - 
  mean(intrinsic_variables$AgeYears, na.rm = TRUE)) / 5

intrinsic_variables$extreme_c <- (intrinsic_variables$n_extreme_both - 
  mean(intrinsic_variables$n_extreme_both, na.rm = TRUE)) / 5

intrinsic_variables$density_c <- (intrinsic_variables$avg_density - 
                                    mean(intrinsic_variables$avg_density, na.rm = TRUE)) / 5

mod_int_2016_2025 <- glmer(proportion ~ AgeYears * avg_density + (1 | animalID_fct) + (1 | season_fct),
                           weights = total_resights,
                           family = binomial(link= "logit"),
                           control = glmerControl(optimizer = "bobyqa"),
                           data = intrinsic_variables); summary(mod_int_2016_2025)
ranef(mod_int_2016_2025)
exp(fixef(mod_int_2016_2025)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_int_2016_2025, plot = TRUE) #plot residuals
check_collinearity(mod_int_2016_2025) #check predictor VIFs

pred_int_density <- ggpredict(
  mod_int_2016_2025,
  terms = c("avg_density [all]", "AgeYears [all]"))

ggplot(pred_int_density, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  labs(
    x = "Average location density",
    y = "Predicted MOA",
    color = "AgeYears",
    fill = "AgeYears"
  ) +
  theme_minimal()

mod_int_1996_2023 <- glmer(proportion ~ age_cat : age10 * n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                           weights = total_resights,
                           family = binomial(link = "logit"),
                           control = glmerControl(optimizer = "bobyqa"),
                           data = intrinsic_variables); summary(mod_int_1996_2023)
ranef(mod_int_1996_2023)
exp(fixef(mod_int_1996_2023)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_int_1996_2023, plot = TRUE) #plot residuals
check_collinearity(mod_int_1996_2023) #check predictor VIFs

### model with interaction between age*n_extreme_both with all predictors ###
mod_int_weather_2016_2023 <- glmer(proportion ~ AgeYears * extreme_c + avg_density + (1 | animalID_fct) + (1 | season_fct),
                                   weights = total_resights,
                                   family = binomial(link = "logit"),
                                   control = glmerControl(optimizer = "bobyqa"),
                                   data = intrinsic_variables); summary(mod_int_weather_2016_2023)
ranef(mod_int_weather_2016_2023)
exp(fixef(mod_int_weather_2016_2023)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_int_weather_2016_2023, plot = TRUE) #plot residuals
check_collinearity(mod_int_weather_2016_2023) #check predictor VIFs

# 1) Predictions for the interaction
pred_mod_int <- ggpredict(mod_int_weather_2016_2023,
                          terms = c("AgeYears [all]", "extreme_c [all]"))

# make numeric for gradient
pred_mod_int$n_extreme_both_num <- parse_number(as.character(pred_mod_int$group))

# 2) Plot: ribbon + line, color gradient by n_extreme_both intensity
ggplot(pred_mod_int, aes(x = x, y = predicted, group = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = n_extreme_both_num),
              alpha = 0.12, color = NA) +
  geom_line(aes(color = n_extreme_both_num),
            linewidth = 0.9, alpha = 0.85) +
  scale_colour_gradientn(name = "n_extreme_both",
                         colors = c("pink", "#f768a1", "#ae017e")) +
  scale_fill_gradientn(colors = c("pink", "#f768a1", "#ae017e"),
                       guide = "none") +
  labs(x = "Maternal age", y = "Predicted MOA") +
  theme_minimal()

pred_mod_int$AgeYears <- parse_number(as.character(pred_mod_int$x))

ggplot(pred_mod_int, aes(x = n_extreme_both_num, y = predicted, group = x)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.08, color = NA) +
  geom_line(aes(color = Age_num),
            linewidth = 1, alpha = 0.9) +
  labs(x = "Number of extreme events",
       y = "Mother-offspring association") +
  theme_minimal()

### model with interaction between age*density with all predictors ###
mod_int_density_2016_2023 <- glmer(proportion ~ AgeYears * avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                   weights = total_resights,
                                   family = binomial(link = "logit"),
                                   control = glmerControl(optimizer = "bobyqa"),
                                   data = intrinsic_variables); summary(mod_int_density_2016_2023)
ranef(mod_int_density_2016_2023)
exp(fixef(mod_int_density_2016_2023)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_int_density_2016_2023, plot = TRUE) #plot residuals
check_collinearity(mod_int_density_2016_2023) #check predictor VIFs

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

#################### 4a) 1996-2025 linear vs. quadratic vs. piecewise model AIC comparison ########################

### threshold ###

mod_binom_1996_2025 <- glmer(proportion ~ age_cat : age10 + (1 | animalID_fct) + (1 | season_fct),
                             weights = total_resights,
                             family = binomial(link= "logit"),
                             control = glmerControl(optimizer = "bobyqa"),
                             data = intrinsic_variables); summary(mod_binom_1996_2025) #model summary

### linear ###

mod_linear <- glmmTMB(proportion ~ AgeYears + (1 | animalID_fct) + (1 | season_fct),
                      weights = total_resights,
                      family = binomial(link = "logit"),
                      data = intrinsic_variables); summary(mod_linear)

### quadratic ###

mod_quad <- glmmTMB(proportion ~ age10 + I(age10^2) + (1 | animalID_fct) + (1 | season_fct),
                    weights = total_resights,
                    family = binomial(link = "logit"),
                    data = intrinsic_variables); summary(mod_quad)

## make table with model type comparisons

mod_comparisons <- list(Linear = mod_linear,
                        Threshold = mod_binom_1996_2025,
                        Quadratic = mod_quad)

aic_table_1996_2025 <- tibble(Model = names(mod_comparisons),
                              AIC = sapply(mod_comparisons, AIC),
                              K = sapply(mod_comparisons, function(m) attr(logLik(m), "df"))) %>%
  mutate(delta_AIC = AIC - min(AIC),
         aic_weights = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(delta_AIC)

aic_table_1996_2025 <- aic_table_1996_2025 %>%
  transmute(Model = Model, 
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(aic_weights, 3)) %>%
  arrange("ΔAIC")

best_model <- aic_table_1996_2025$Model[1]

#bold lowest AIC model
aic_table_1996_2025 <- flextable(aic_table_1996_2025) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(aic_table_1996_2025$Model == best_model),
       part = "body"); aic_table_1996_2025

#save final table
save_as_docx(aic_table_1996_2025, path = "./TablesFigures/1996_2025_Age_Predictor_Comparison.docx")

########## 4b) 2016-2023 linear vs. quadratic vs. piecewise model AIC comparison ###################

### threshold ###

mod_binom_2016_2023 <- glmer(proportion ~ age_cat : age10 + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                             weights = total_resights,
                             family = binomial(link= "logit"),
                             control = glmerControl(optimizer = "bobyqa"),
                             data = intrinsic_variables); summary(mod_binom_2016_2023) #model summary

### linear ###

mod_linear_2016_2023 <- glmmTMB(proportion ~ AgeYears + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                weights = total_resights,
                                family = binomial(link = "logit"),
                                data = intrinsic_variables); summary(mod_linear_2016_2023)

### quadratic ###

mod_quad_2016_2023 <- glmmTMB(proportion ~ age10 + I(age10^2) + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                    weights = total_resights,
                    family = binomial(link = "logit"),
                    data = intrinsic_variables); summary(mod_quad_2016_2023)


## name models to compare
mod_comparisons <- list(Linear = mod_linear_2016_2023,
                        Quadratic = mod_quad_2016_2023,
                        Threshold = mod_binom_2016_2023)

# 2) make table with AIC comparisons
aic_table_2016_2023 <- tibble(Model = names(mod_comparisons),
                              AIC = sapply(mod_comparisons, AIC)) %>%
  mutate(delta_AIC = AIC - min(AIC),
         aic_weights = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(delta_AIC)

aic_table_2016_2023 <- aic_table_2016_2023 %>%
  transmute(Model = Model,
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(aic_weights, 3)) %>%
  arrange("ΔAIC")

best_model <- aic_table_2016_2023$Model[1]

#bold lowest AIC model
aic_table_2016_2023 <- flextable(aic_table_2016_2023) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(aic_table_2016_2023$Model == best_model),
       part = "body"); aic_table_2016_2023

#save final table
save_as_docx(aic_table_2016_2023, path = "./TablesFigures/2016_2023_Age_Predictor_Comparison.docx")

############ 4c) age threshold AIC comparison table ####################

#test all possible thresholds
cutoff <- 5:14

# For each cutoff:
# 1) define Young/Old at a, and center age at a
# 2) fit the same model

threshold_test <- function(a) {
  
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
aic_threshold_comparison <- map_dfr(cutoff, threshold_test) %>%
  mutate(delta_AIC = AIC - min(AIC),
         AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(Threshold)

threshold_comparison_tbl <- aic_threshold_comparison %>%
  transmute(Threshold = Threshold,
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(AIC_weight, 3)) %>%
  arrange("Threshold")

#bold the lowest AIC threshold
threshold_comparison_tbl <- flextable(threshold_comparison_tbl) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(threshold_comparison_tbl$Threshold == 9), part = "body"); threshold_comparison_tbl

save_as_docx(threshold_comparison_tbl, path = "./TablesFigures/AIC_Threshold_Comparison.docx")

############# 4d) age threshold significance comparison figure #################

#test all possible thresholds
cutoff <- 5:14

# For each cutoff:
# 1) define Young/Old at a, and center age at a
# 2) fit the same model
# 3) pull the post-threshold slope term (Old × age10), its SE, and p-value
thr_res <- map_dfr(cutoff, \(a){
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
threshold_comparison_figure <- ggplot(thr_res, aes(age_cutoff, coef)) +
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
  labs(x = "Senescence threshold",
       y = "Estimated old slope (log-odds)") +
  theme_classic(); threshold_comparison_figure

ggsave("./TablesFigures/threshold_comparison_figure.png", threshold_comparison_figure, width = 8, height = 10, dpi = 600)

############### Individual level consistency figures (individual ID effects) ###############

# 1) Build a dataframe pre-merge for modelling to calculate individual-level consistency
season_level <- Proportion_MPA %>% #has animalID, season, count_1_pup, total_resights
  mutate(proportion = count_1_pup / total_resights) %>% #Includes ALL data
  left_join(metadata %>% select(animalID, season, AgeYears),
            by = c("animalID","season")) %>%
  distinct() %>% #reduces to 1 animalID/season combo per observation (no duplicates)
  mutate(animalID_fct = factor(animalID),
         season_fct = factor(season))

#helper for weighted variance calculations
wvar <- function(x, w) {
  w <- w / sum(w) # (1) Normalize weights so they sum to 1
  mu <- sum(w * x) # (2) Compute the weighted mean of x
  sum(w * (x - mu)^2) # (3) Compute weighted average squared deviation from mean
}

# 2) Calculate cosnsitency score for MPA for each individual
##Summarize per individual:
#  - n_seasons: number of seasons observed
#  - n_obs_total: total number of resights
#  - mean_assoc: weighted mean proportion of association (0–1)
#  - sd_assoc: weighted SD across seasons, indicating behavioral variability
#  - frac_zero / frac_one: fraction of seasons where proportion = 0 or 1
#  - consistency: rescaled measure (1 − SD/0.5), bounded between 0–1 where 1 = perfectly consistent, 0 = highly variable
individual_consistency <- season_level %>%
  group_by(animalID_fct) %>%
  mutate(n_seasons   = n(),
         n_obs_total = sum(total_resights),
         mean_assoc  = weighted.mean(proportion, w = total_resights, na.rm = TRUE),     
         sd_assoc    = sqrt(wvar(proportion, w = total_resights)),     
         frac_zero   = mean(proportion == 0),
         frac_one    = mean(proportion == 1),
         .groups = "drop") %>%
  mutate(
    #Consistency: standardized so that SD = 0 (no variation) → 1, and SD = 0.5 (max) → 0
    consistency = pmax(pmin(1 - (sd_assoc / 0.5), 1), 0)   # 0–1, higher = more consistent
  ) %>%
  filter(n_seasons >= 3) #require at least 3 breeding seasons per individual

### Foundational consistency figures (allows us to look at specific animalIDs) ###

##make a data table to create a star label for animals with consistent prop == 1
star_data <- individual_consistency %>%
  group_by(animalID_fct) %>%
  summarize(all_one = all(sd_assoc == 0),
            max_y = max(proportion)) %>%
  filter(all_one) %>%
  mutate(y_position = max_y + 0.05)

## boxplot with color gradient for standard deviation + stars for "highly consistent" MOA indiciduals
ggplot(individual_consistency, aes(x = animalID_fct, y = proportion, fill = sd_assoc)) +
  geom_boxplot(width = 0.7, outlier.size = 0.6, outlier.color = "#04BBB2") +
  geom_jitter(size = 1, alpha = 0.6, color = "#04BBB2") +
  scale_fill_gradient(low = "#D8FFF7", high = "#073481", name = "SD of Proportion") +
  geom_text(data = star_data, aes(x = animalID_fct, y = y_position, label = "*"), color = "#C699E1", size = 7, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 1.03)) +
  scale_y_continuous(n.breaks = 10) +
  theme_few() +
  labs(x = "Individual ID", 
       y = "Mother-offspring association") +
  theme(axis.text.x = element_text(angle = 90, size = 6))

############ Age-corrected consistency figure (Maddie's version) ###########

## 1) Compute age-specific means and deviations (season level)
#   For each age (AgeYears), compute the mean maternal association across all seals,
#   then calculate each observation's deviation from that age-specific mean.
season_level_age_centered <- individual_consistency %>%
  group_by(AgeYears) %>%
  mutate(
    mean_prop_age = mean(proportion, na.rm = TRUE),          # population mean at this age
    dev_from_age  = proportion - mean_prop_age               # deviation from age-specific mean
  ) %>%
  ungroup()

# Visualize within-age variation in maternal association w/ boxplots for each age
ggplot(season_level_age_centered,
       aes(x = factor(AgeYears), y = dev_from_age)) +
  geom_boxplot(fill = "#A1D3B2", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +           # 0 = age-specific mean
  labs(
    x = "Age (years)",
    y = "Deviation from age-specific mean",
    title = "Within-age variation in maternal association"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))


## 1) Mean age per individual (weighted by resights)
#   Compute a weighted mean age per individual, weighting by total resights,
#   so individuals with more resight effort contribute more to their mean age.
age_by_ind <- season_level_age_centered %>%
  group_by(animalID_fct) %>%
  summarise(
    mean_age = weighted.mean(AgeYears, w = total_resights, na.rm = TRUE),
    .groups = "drop"
  )

## 2) Classify direction and consistency of deviation
#   dev_eps: threshold for how far from 0 a mean deviation must be to count as
#   "meaningfully above" or "meaningfully below" age peers.
dev_eps      <- 0.02   # deviation > 0.02 or < -0.02 is considered meaningful

# cons_thresh: threshold for high consistency in deviation across years.
# (e.g., consistency_dev >= 0.8)
cons_thresh  <- 0.7   

ind_consistency_age <- season_level_age_centered %>%
  group_by(animalID_fct) %>%
  summarize(
    n_seasons   = n_distinct(season),
    n_obs_total = sum(total_resights, na.rm = TRUE),
    mean_assoc  = mean(proportion, na.rm = TRUE),
    mean_dev    = mean(dev_from_age, na.rm = TRUE),
    sd_dev      = sd(dev_from_age, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Consistency metric (0..1): high when sd_dev is small relative to signal
    consistency_dev = 1 - (sd_dev / (sd_dev + abs(mean_dev) + 1e-6))
  )

# Using the individual-level object ind_consistency_age, which already contains:
# - mean_dev: mean of dev_from_age across seasons
# - consistency_dev: consistency in that deviation across seasons
# - mean_assoc: lifetime mean maternal association

ind_consistency_age <- ind_consistency_age %>%
  mutate(
    # Direction of mean deviation relative to age peers
    dev_dir = case_when(
      mean_dev >  dev_eps  ~ "Above peers",
      mean_dev < -dev_eps  ~ "Below peers",
      TRUE                 ~ "Near age mean"
    ),
    
    # Consistency of that deviation across years
    dev_consistency_class = case_when(
      dev_dir == "Above peers"    & consistency_dev >= cons_thresh ~ "Consistently above peers",
      dev_dir == "Below peers"    & consistency_dev >= cons_thresh ~ "Consistently below peers",
      dev_dir == "Near age mean"  & consistency_dev >= cons_thresh ~ "Consistently near mean",
      TRUE ~ "Inconsistent"  # fluctuates around peers, no stable position
    ),
    
    # Combine direction + consistency into a single age-normalized strategy label
    strategy_cross = case_when(
      dev_consistency_class == "Consistently above peers" ~ "Consistently above peers",
      dev_dir == "Above peers"    & dev_consistency_class == "Inconsistent" ~ "Inconsistent above peers",
      dev_dir == "Below peers"    & dev_consistency_class == "Inconsistent" ~ "Inconsistent below peers",
      dev_dir == "Near age mean"  & dev_consistency_class == "Inconsistent" ~ "Inconsistent near mean",
      TRUE ~ "Other"  # catch-all; should ideally be rare or empty
    ),
    
    # Order strategy levels for plotting (top: specialists, bottom: below peers)
    strategy_cross = factor(
      strategy_cross,
      levels = c(
        "Consistently above peers",
        "Inconsistent above peers",
        "Inconsistent near mean",
        "Inconsistent below peers"
      )
    )
  )


## 3) Summary statistics for strategy classes
# Overall summary: how many individuals in each consistency class, and their mean maternal association
summary_overall <- ind_consistency_age %>%
  group_by(dev_consistency_class) %>%
  summarise(
    n         = n(),                                 # number of individuals
    prop      = n / sum(n),                          # proportion of sample
    mean_assoc = mean(mean_assoc, na.rm = TRUE),     # average maternal association
    .groups   = "drop"
  ); summary_overall

# Direction-only summary: how many individuals are on average above, below, etc. and their mean MPA
summary_dir <- ind_consistency_age %>%
  group_by(dev_dir) %>%
  summarise(
    n         = n(),
    prop      = n / sum(n),
    mean_assoc = mean(mean_assoc, na.rm = TRUE),
    .groups   = "drop"
  ); summary_dir

# Cross-tab: direction x consistency class: how many individuals are in different strategies
summary_cross <- ind_consistency_age %>%
  group_by(dev_dir, dev_consistency_class) %>%
  summarise(
    n               = n(),
    prop_within_dir = n / sum(n),                    # proportion within each dev_dir
    mean_assoc      = mean(mean_assoc, na.rm = TRUE),
    .groups         = "drop"
  ); summary_cross


## 4) Join individual-level consistency back to season rows
individual_consistency_age_rows <- season_level_age_centered %>%
  left_join(
    ind_consistency_age %>%
      select(animalID_fct, consistency_dev),
    by = "animalID_fct"
  )

names(individual_consistency_age_rows)


## 5) Boxplot of within-individual variation colored by consistency
#   For each individual, show distribution of maternal association (proportion) across seasons
ggplot(individual_consistency_age_rows,
       aes(x = animalID_fct, y = proportion, fill = consistency_dev)) +
  geom_boxplot(width = 0.7, outlier.size = 1.5, outlier.color = "black") +
  geom_jitter(size = 1, alpha = 0.6, color = "#04BBB2") +
  scale_fill_gradient(
    low  = "#D8FFF7",
    high = "#073481",
    name = "Age-normalized\nconsistency"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(n.breaks = 10) +
  theme_few() +
  labs(
    x = "Animal ID",
    y = "Proportion maternal association"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = "right"
  )


## 6) Strategy scatterplot in "strategy space"
#   Create a compact dataframe for plotting strategy classes
ind_consistency_age_strat <- ind_consistency_age %>%
  select(animalID_fct, mean_assoc, mean_dev, consistency_dev, n_obs_total, strategy_cross)

# Plot individuals in "strategy space":
#   - x-axis: mean maternal association (0–1)
#   - y-axis: age-normalized consistency
#   - vertical dashed lines: thresholds for mean association
#   - horizontal dashed line: threshold for high consistency
consistency_figure <- ggplot(ind_consistency_age_strat,
                             aes(x = mean_assoc, y = consistency_dev, color = strategy_cross)) +
  geom_point(aes(size = n_obs_total, alpha = 0.8)) +
  ## geom_vline(xintercept = 0, size = 0.6, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = cons_thresh, size = 0.6, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c(
    "Consistently above peers"  = "#1f78b4",
    "Inconsistent above peers"  = "#6F73D2",
    "Inconsistent near mean"    = "#66C2A5",
    "Inconsistent below peers"  = "#B74F6F"
  )) +
  scale_size_continuous(name = "Total observations", range = c(2, 10)) +
  labs(x     = "Mean mother-offspring association",
       y     = "Consistency in deviation from age mean",
       color = "Age-normalized\nstrategy class",
       title = "Age-normalized maternal strategies",
       subtitle = "Deviation from age peers (x-axis) and consistency in that deviation (y-axis)") +
  guides(alpha = "none") +
  theme_few(base_size = 14) +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.title = element_text(face = "bold")); consistency_figure

ggsave("./TablesFigures/consistency_figure.png", consistency_figure, width = 12, height = 8, dpi = 600)

## check who the most consistent individuals are!
ind_consistency_age$animalID_fct[ind_consistency_age$strategy_cross == "Consistently above peers"]

ind_consistency_age_perfect <- ind_consistency_age %>%
  filter(strategy_cross == "Consistently above peers")

############ Age-normalized consistency figure (Aditi's pass) ###########

# 1) helper for weighted variance
wvar <- function(x, w) {
  w <- w / sum(w)
  mu <- sum(w * x)
  sum(w * (x - mu)^2)
}

# 2) Compute deviation from age-specific mean:
# For each age, calculate the weighted average MOA,
# Calculate how far above/below its age peers it is
season_level_age_centered <- intrinsic_variables %>%
  select(animalID_fct, season, AgeYears, proportion, total_resights) %>%
  group_by(AgeYears) %>%
  mutate(mean_prop_age = weighted.mean(proportion, w = total_resights, na.rm = TRUE),
         dev_from_age  = proportion - mean_prop_age) %>%
  ungroup()

# 3) Individual consistency in the deviation:
# For each animal, calculate lifetime weighted mean deviation from age-specific mean MOA
# Calculate the consistency or sd of that deviation using wvar
ind_consistency_age <- season_level_age_centered %>%
  group_by(animalID_fct) %>%
  summarise(n_seasons   = n_distinct(season),
            n_obs_total = sum(total_resights, na.rm = TRUE),
            mean_dev = weighted.mean(dev_from_age, w = total_resights, na.rm = TRUE),
            sd_dev = sqrt(wvar(dev_from_age, w = total_resights)),
            .groups = "drop") %>%
  filter(n_seasons >= 3) %>% #only animals observed at least 3 times in lifetime
  mutate(consistency_dev = exp(-sd_dev)) #rescale so consistency is 0-1 and 1 is highly consistent

##Arbitrary thresholds (can change)
dev_eps     <- 0.05 #change of 0.05 in dev from mean is different from mean
cons_thresh <- 0.90 #above 0.9 is highly consistent

# 4) Individual mean and consistency classifications based on the thresholds
ind_consistency_age <- ind_consistency_age %>%
  mutate(dev_dir = case_when(
    mean_dev >  dev_eps ~ "Above peers",
      mean_dev < -dev_eps ~ "Below peers",
      TRUE                ~ "Near age mean"
    ),
    dev_consistency_class = case_when(
      dev_dir == "Above peers"   & consistency_dev >= cons_thresh ~ "Consistently above peers",
      dev_dir == "Below peers"   & consistency_dev >= cons_thresh ~ "Consistently below peers",
      dev_dir == "Near age mean" & consistency_dev >= cons_thresh ~ "Consistently near mean",
      TRUE ~ "Inconsistent" #all others are inconsistent
    ),
    strategy_cross = case_when(
      dev_consistency_class == "Consistently above peers" ~ "Consistently above peers",
      dev_consistency_class == "Consistently below peers" ~ "Consistently below peers",
      dev_consistency_class == "Consistently near mean"   ~ "Consistently near mean",
      dev_dir == "Above peers"   & dev_consistency_class == "Inconsistent" ~ "Inconsistent above peers",
      dev_dir == "Below peers"   & dev_consistency_class == "Inconsistent" ~ "Inconsistent below peers",
      dev_dir == "Near age mean" & dev_consistency_class == "Inconsistent" ~ "Inconsistent near mean"
    ) 
  ) %>%
  mutate(
    strategy_cross = factor(
      strategy_cross,
      levels = c(
        "Consistently above peers",
        "Consistently near mean",
        "Consistently below peers",
        "Inconsistent above peers",
        "Inconsistent near mean",
        "Inconsistent below peers"
      )
    )
  )

max_abs <- max(abs(ind_consistency_age$mean_dev), na.rm = TRUE) #use to rescale x-axis so each side is even

# 6) Consistency figure: 
#mean deviation on the x, consistency in that deviation on y
#colors for maternal "strategy"
#dashed lines at 0 (no dev from the mean) and 0.9 (above 0.9 is consistent)
#points weighted by total lifetime resights per animal
consistency_figure <- ggplot(ind_consistency_age,
                             aes(x = mean_dev, y = consistency_dev, color = strategy_cross)) +
  geom_point(aes(size = n_obs_total), alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = cons_thresh, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c(
    "Consistently above peers" = "#1f78b4",
    "Consistently near mean"   = "#4DB6AC",
    "Consistently below peers" = "orange",
    "Inconsistent above peers" = "#6F73D2",
    "Inconsistent near mean"   = "pink",
    "Inconsistent below peers" = "#D17A92",
    "Other"                    = "grey60"
  )) +
  scale_size_continuous(name = "Total observations", range = c(2, 10)) +
  labs(
    x = "Deviation from age-specific mean MOA",
    y = "Consistency of deviation",
    color = "Strategy class",
    title = "Age-normalized maternal strategies",
    subtitle = "Mean deviation from age peers (x) and consistency in deviation (y)"
  ) +
  guides(alpha = "none") +
  theme_few(base_size = 14) +
  scale_x_continuous(limits = c(-max_abs, max_abs)) +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.title = element_text(face = "bold")
  ); consistency_figure

ggsave("./TablesFigures/consistency_figure.png",
       consistency_figure,
       width = 12, height = 8, dpi = 600)

