##run the source code from Data Processing
source("./DataCurationCode/Data_processing_MPA.R")

#### These are the main models Maddie and I think make most sense:

##The first model includes all predictors with a smaller subset of data (2016-2023)
##The second model allows us to test for senescence/non-linearity of the age predictor using the full dataset (1996-2025)

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

# 3) Compare base model to model with experience instead of age 
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

### Linear age figure ###

# 2) Define x-axis ages and seasons used for season-specific curves
ages <- sort(unique(intrinsic_2016_2023$AgeYears)) #unique ages present (x values)
seasons <- levels(intrinsic_2016_2023$season_fct) #season factor levels used in predictions

# 3) Compute pooled observed proportion by age (black points + binomial CIs)
observed_data_2016_2023 <- intrinsic_2016_2023 %>%
  group_by(AgeYears) %>% #group by rows within each age
  summarise(n_age = n(), #sample size at each age for labels
            n_success = sum(count_1_pup, na.rm = TRUE), #total successes
            n_trials = sum(total_resights, na.rm = TRUE), #total effort
            avg_prop = n_success/n_trials, #observed proportion
            lwr = binom.test(n_success, n_trials)$conf.int[1], #binomial CI lower
            upr = binom.test(n_success, n_trials)$conf.int[2], #binomial CI upper
            .groups = "drop")

# 4) Build standardized newdata for each age for post-stratification
nd_by_age <- lapply(ages, \(A) transform(intrinsic_2016_2023, #retain covariate distribution
                                         AgeYears = A)) #set all seals to age A

# 5) Post-stratified prediction curve (row predictions → weighted mean)
post_curve <- function(fit, season = NULL){
  vapply(seq_along(nd_by_age), function(i){
    nd <- nd_by_age[[i]] #dataset for one age
    if(!is.null(season)) nd$season_fct <- factor(season, levels = seasons) #force season level
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL) #predict including all RE
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE) #total resight weighted mean
  }, numeric(1))
}

# 6) Thin season curves (one post-stratified curve per season)
season_lines_2016_2023 <- bind_rows(lapply(seasons, \(S)
                                           tibble(season_fct = S, #season grouping
                                                  AgeYears = ages, #x-axis
                                                  pred = post_curve(mod_binom_2016_2023, season = S)))) #season-specific predictions

# 7) Overall curve + parametric bootstrap CI
set.seed(123) #set for reproducibility
nsim <- 100 #bootstrap draws
post_curve_overall <- function(fit) post_curve(fit, season = NULL) #overall curve

overall_2016_2023 <- tibble(AgeYears = ages,
                            pred = post_curve_overall(mod_binom_2016_2023)) %>% #mean prediction
  left_join(observed_data_2016_2023 %>% select(AgeYears, n_age), by = "AgeYears") #add sample sizes

boot_2016_2023 <- bootMer(mod_binom_2016_2023, FUN = post_curve_overall,
                          nsim = nsim, type = "parametric", use.u = FALSE) #parametric bootstrap

# 8) Convert bootstrap draws into 95% CI ribbon
ci_2016_2023 <- overall_2016_2023 %>%
  mutate(lo = apply(boot_2016_2023$t, 2, quantile, 0.025, na.rm = TRUE), #lower CI
         hi = apply(boot_2016_2023$t, 2, quantile, 0.975, na.rm = TRUE)) #upper CI

# 9) Plot: thin season curves + bootstrap ribbon + mean curve + observed data
ggplot() +
  geom_line(data = season_lines_2016_2023, #thin season-specific curves
            aes(AgeYears, pred, group = season_fct),
            color = "navy", alpha = 0.25) +
  geom_ribbon(data = ci_2016_2023, #bootstrap CI band
              aes(AgeYears, ymin = lo, ymax = hi),
              fill = "#92BAEE", alpha = 0.28) +
  geom_line(data = ci_2016_2023, #overall predicted curve
            aes(AgeYears, pred),
            color = "#7299D3", linewidth = 1.2) +
  geom_pointrange(data = observed_data_2016_2023, #observed pooled proportions
                  aes(AgeYears, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_text(data = observed_data_2016_2023, #sample size labels
            aes(AgeYears, 1.01, label = n_age),
            vjust = -0.5) +
  scale_x_continuous(breaks = ages) +
  coord_cartesian(ylim = c(0.75, 1.01), clip = "off") +
  theme_few() +
  labs(x = "Age (Years)", y = "Mother-offspring association")

### plot effects of density ###

# 1) predicted effect of density on flipped proportion
pred_density_2016_2023 <- ggpredict(mod_binom_2016_2023, terms = "avg_density [all]")

# 2) plot flipped proportion and average density
ggplot(data = pred_density_2016_2023, aes(x = x, y = predicted)) +
  
  #Raw points
  geom_point(data = intrinsic_2016_2023,
             aes(x = avg_density, y = proportion),
             color = "#583478",
             alpha = 0.4,
             size = 2,
             inherit.aes = FALSE) +
  
  #Confidence ribbon
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "#5E397F",
              alpha = 0.3) +
  
  #Prediction line
  geom_line(color = "#976ABF",
            linewidth = 1.2) +
  coord_cartesian(ylim = c(0, 1), expand = TRUE) +
  labs(x = "Wthin-colony location density",
       y = "Mother-offspring association") +
  theme_few()

# 1) predicted effect of extreme events by season
pred_extreme_2016_2023 <- ggpredict(mod_binom_2016_2023,
                                    terms = c("n_extreme_both [all]", "season_fct"))

# 2) force season order in the raw data
intrinsic_2016_2023 <- intrinsic_2016_2023 %>%
  mutate(season_fct = factor(season_fct,
                             levels = sort(unique(as.numeric(as.character(season_fct))))))

# 3) force the same order in ggpredict output (its season column is `group`)
pred_extreme_2016_2023$group <- factor(pred_extreme_2016_2023$group, levels = seasons)

# 4) pink→green palette (named + ordered)
pal <- setNames(colorRampPalette(c("#D1497A", "#E07A9A", "#7FC97F", "#2E8B57"))(length(seasons)),
                seasons) #season-specific colors

# 5) plot: points + ribbons + lines, with ordered legend
ggplot(pred_extreme_2016_2023, aes(x = x, y = predicted, group = group, color = group, fill = group)) +
  geom_point(data = intrinsic_2016_2023,
             aes(x = n_extreme_both, y = proportion, color = season_fct),
             alpha = 0.6, size = 2, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = pal, breaks = seasons, name = "Year") +
  scale_fill_manual(values = pal, breaks = seasons, name = "Year") +
  labs(x = "Number of extreme wave and tide events",
       y = "Mother-offspring association") +
  theme_few()

############# 1b) 1996-2025 full dataset model for age ###############

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
            lwr = binom.test(n_success, n_trials)$conf.int[1], #exact binomial CI lower bound
            upr = binom.test(n_success, n_trials)$conf.int[2], #exact binomial CI upper bound
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
  vapply(seq_along(nd_by_age), function(i){ #loop over ages
    nd <- nd_by_age[[i]] #newdata for one age
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

boot_1996_2025 <- bootMer(mod_binom_1996_2025, FUN = post_curve_overall, nsim = nsim, #bootstrap curves by refitting/simulating
                          type = "parametric", use.u = FALSE) #parametric bootstrap = resimulate random effects each time

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
ggplot() +
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
  geom_pointrange(data = observed_data_1996_2025, #pooled observed points with exact binomial CI
                  aes(AgeYears, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  annotate("text", 
           x = age_senesce - 0.5,
           y = 0.78, label = "Age threshold (9 years)",
           angle = 90, vjust = -0.8, size = 3.5) +
  coord_cartesian(ylim = c(0.75, 1.01), clip = "off") +
  geom_text(data = observed_data_1996_2025, aes(AgeYears, 1.01, label = n_age), vjust = -0.5) + #sample size labels per age
  scale_color_manual(name = "Age class", values = COL) +
  scale_fill_manual(name = "Age class", values = COL) +
  scale_x_continuous(breaks = ages) +
  theme_few() +
  labs(x = "Age (Years)", y = "Mother-offspring association")

#################### 1c) Interaction models with age, extreme events, density ##########################

mod_int_2016_2025 <- glmer(proportion ~ scale(AgeYears) * scale(avg_density) + (1 | animalID_fct) + (1 | season_fct),
                           weights = total_resights,
                           family = binomial(link= "logit"),
                           control = glmerControl(optimizer = "bobyqa"),
                           data = intrinsic_variables); summary(mod_int_2016_2025)
ranef(mod_int_2016_2025)
exp(fixef(mod_int_2016_2025)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_int_2016_2025, plot = TRUE) #plot residuals
check_collinearity(mod_int_2016_2025) #check predictor VIFs

mod_int_1996_2023 <- glmer(proportion ~ age10 : age_cat * scale(n_extreme_both) + (1 | animalID_fct) + (1 | season_fct),
                           weights = total_resights,
                           family = binomial(link = "logit"),
                           control = glmerControl(optimizer = "bobyqa"),
                           data = intrinsic_variables); summary(mod_int_1996_2023)
ranef(mod_int_1996_2023)
exp(fixef(mod_int_1996_2023)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_int_1996_2023, plot = TRUE) #plot residuals
check_collinearity(mod_int_1996_2023) #check predictor VIFs

### model with interaction between age*n_extreme_both with all predictors ###
mod_int_weather_2016_2023 <- glmer(proportion ~ scale(AgeYears) * scale(n_extreme_both) + scale(avg_density) + (1 | animalID_fct) + (1 | season_fct),
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
                          terms = c("AgeYears [all]", "n_extreme_both [all]"))

pred_df$n_extreme_both_num <- parse_number(as.character(pred_df$group))

# 2) Plot: ribbon + line, color gradient by n_extreme_both intensity
ggplot(pred_df, aes(x = x, y = predicted, colour = n_extreme_both_num, group = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = n_extreme_both_num),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  scale_colour_gradientn(name = "n_extreme_both", colors = c("pink", "#f768a1", "#ae017e")) +
  scale_fill_gradientn(colors = c("pink", "#f768a1", "#ae017e"), guide = "none") +
  labs(x = "Maternal age", y = "Mother-offspring association (predicted probability)") +
  theme_minimal()

### model with interaction between age*density with all predictors ###
mod_int_density_2016_2023 <- glmer(proportion ~ scale(AgeYears) * scale(avg_density) + scale(n_extreme_both) + (1 | animalID_fct) + (1 | season_fct),
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

#################### 2a) 1996-2025 linear vs. quadratic vs. piecewise model AIC comparison ########################

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
save_as_docx(aic_table_1996_2025, path = "1996_2025_Age_Predictor_Comparison.docx")

########## 2b) 2016-2023 linear vs. quadratic vs. piecewise model AIC comparison ###################

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
save_as_docx(aic_table_2016_2023, path = "2016_2023_Age_Predictor_Comparison.docx")

############ 3a) age threshold AIC comparison table ####################

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

save_as_docx(threshold_comparison_tbl, path = "AIC_Threshold_Comparison.docx")

############# 3b) age threshold significance comparison figure #################

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
  tab <- as.data.frame(summary(m)$coefficients$cond) %>%
    rownames_to_column("term")
  
  #interaction term = slope of age10 for Old
  term_name <- "age_catOld : age10"
  
  row <- filter(tab, term == term_name)
  
  tibble(age_cutoff = a,
         coef = row$Estimate,
         se = row$`Std. Error`,
         p = 2 * pnorm(abs(row$`z value`), lower.tail = FALSE))
}) %>%
  mutate(si = p <= 0.05, # flag significance
         logp_cap = pmin(-log10(p), 5)) # point size weighted by p

# Plot: coefficient ± 95% CI; bigger points = smaller p; color = significance
ggplot(thr_res, aes(age_cutoff, coef)) +
  geom_hline(yintercept = 0, linetype = "dashed") +  # no post-threshold slope
  geom_pointrange(aes(ymin = coef - 1.96 * se,
                      ymax = coef + 1.96 * se,
                      color = sig,
                      size = logp_cap)) +
  scale_color_manual(values = c(`TRUE` = "deepskyblue", `FALSE` = "tomato"),
                     labels = c(`TRUE` = "p ≤ 0.05", `FALSE` = "p > 0.05"),
                     name = "Significance") +
  scale_size_continuous(range = c(0.5, 2.2), guide = "none") +
  scale_x_continuous(breaks = cutoff) +
  labs(x = "Senescence threshold",
       y = "Estimated old slope (log-odds)") +
  theme_classic()

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

############ Age-corrected consistency figure ###########

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
  select(animalID_fct, mean_assoc, consistency_dev, n_obs_total, strategy_cross)

# Plot individuals in "strategy space":
#   - x-axis: mean maternal association (0–1)
#   - y-axis: age-normalized consistency
#   - vertical dashed lines: thresholds for mean association
#   - horizontal dashed line: threshold for high consistency
ggplot(ind_consistency_age_strat,
       aes(x = mean_assoc, y = consistency_dev, color = strategy_cross)) +
  geom_point(aes(size = n_obs_total, alpha = 0.8)) +
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = cons_thresh, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c(
    "Consistently above peers"  = "#1f78b4",
    "Inconsistent above peers"  = "#6F73D2",
    "Inconsistent near mean"    = "#66C2A5",
    "Inconsistent below peers"  = "#B74F6F"
  )) +
  scale_size_continuous(name = "Total observations", range = c(2, 10)) +
  labs(
    x     = "Mean maternal association (0–1)",
    y     = "Consistency in deviation from age mean",
    color = "Age-normalized\nstrategy class",
    title = "Age-normalized maternal strategies",
    subtitle = "Deviation from age peers (x-axis) and consistency in that deviation (y-axis)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.box      = "vertical",
    legend.title    = element_text(face = "bold")
  )

## check who the most consistent individuals are!
ind_consistency_age$animalID_fct[ind_consistency_age$strategy_cross == "Consistently above peers"]

ind_consistency_age_perfect <- ind_consistency_age %>%
  filter(strategy_cross == "Consistently above peers")
