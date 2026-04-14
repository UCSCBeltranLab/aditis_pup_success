############ 6) Age-normalized among-individual variation (Aditi's pass) ###############

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
  mutate(
    mean_prop_age = mean(proportion, na.rm = TRUE),
    dev_from_age  = proportion - mean_prop_age
  ) %>%
  ungroup()

ind_consistency_age <- season_level_age_centered %>%
  group_by(animalID_fct) %>%
  summarise(
    n_seasons   = n_distinct(season),
    n_obs_total = sum(total_resights, na.rm = TRUE),
    mean_dev    = mean(dev_from_age, na.rm = TRUE),
    sd_dev      = sd(dev_from_age, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_seasons >= 3)

ggplot(ind_consistency_age, aes(x = mean_dev, y = sd_dev)) +
  geom_point(aes(size = n_obs_total), alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    x = "Mean deviation from age-specific mean MOA",
    y = "Within-individual SD of deviation",
    title = "Age-normalized maternal strategies"
  )

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
dev_eps     <- 0.03 #change of 0.03 in dev from mean is different from mean
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
# mean deviation on the x, consistency in that deviation on y
# colors for maternal "strategy"
# dashed lines at 0 (no dev from the mean) and 0.9 (above 0.9 is consistent)
# points weighted by total lifetime resights per animal
consistency_figure <- ggplot(ind_consistency_age,
                             aes(x = mean_dev, y = consistency_dev, color = strategy_cross)) +
  geom_point(aes(size = n_obs_total), alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = cons_thresh, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("Consistently above peers" = "#1f78b4",
                                "Consistently near mean"   = "#4DB6AC",
                                "Consistently below peers" = "#D17A92",
                                "Inconsistent above peers" = "pink",
                                "Inconsistent near mean"   = "#6F73D2",
                                "Inconsistent below peers" = "orange",
                                "Other"                    = "grey60")) +
  scale_size_continuous(name = "Total observations", range = c(2, 10)) +
  labs(x = "Deviation from age-specific mean MOA",
       y = "Consistency of deviation",
       color = "Strategy class",
       title = "Age-normalized maternal strategies",
       subtitle = "Mean deviation from age peers (x) and consistency in deviation (y)") +
  
  annotate("text", x =  max_abs * 0.7, y = cons_thresh + 0.07,
           label = "More consistent",
           size = 4) +
  
  annotate("text", x =  max_abs * 0.7, y = cons_thresh - 0.07,
           label = "Less consistent",
           size = 4) +
  
  annotate("text", x =  max_abs * 0.6, y = 0.6,
           label = "Above age mean",
           size = 4) +
  
  annotate("text", x = -max_abs * 0.6, y = 0.6,
           label = "Below age mean",
           size = 4) +
  guides(alpha = "none") +
  theme_few(base_size = 14) +
  scale_x_continuous(limits = c(-max_abs, max_abs)) +
  coord_cartesian(ylim = c(0.6, 1), clip = "off") +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.title = element_text(face = "bold")); consistency_figure

ggsave("./TablesFigures/consistency_figure.png",
       consistency_figure,
       width = 13, height = 10, dpi = 600)

###################### 6b) Among-individual variation (Aditi's second pass!) ###################

## 1) Add model-based age-specific expectation and Pearson residuals
intrinsic_dev <- intrinsic_variables %>%
  mutate(expected_age_moa = predict(
    mod_binom_1996_2025,
    type = "response",
    re.form = NA),
    dev_from_age_expectation = proportion - expected_age_moa,
    pearson_resid = residuals(mod_binom_1996_2025,
                              type = "pearson"))

## 2) Keep only animals observed at least 3 times
intrinsic_dev_3 <- intrinsic_dev %>%
  group_by(animalID_fct) %>%
  filter(n_distinct(season_fct) >= 3) %>%
  ungroup()

## 3) Summarise individuals
## mean_dev = average above/below age-specific expectation
## sd_resid = within-individual consistency after accounting for mean-variance structure
ind_consistency_desc <- intrinsic_dev_3 %>%
  group_by(animalID_fct) %>%
  summarise(n_obs_total = n(),
            n_seasons = n_distinct(season_fct),
            total_resights_life = sum(total_resights, na.rm = TRUE),
            mean_dev = mean(dev_from_age_expectation, na.rm = TRUE),
            median_dev = median(dev_from_age_expectation, na.rm = TRUE),
            min_dev = min(dev_from_age_expectation, na.rm = TRUE),
            max_dev = max(dev_from_age_expectation, na.rm = TRUE),
            sd_dev = sd(dev_from_age_expectation, na.rm = TRUE),
            sd_resid = sd(pearson_resid, na.rm = TRUE),
            .groups = "drop")

## 4) descriptive labels for "strategies" :) this could be where it gets a little ambiguous
dev_eps <- 0.02

ind_consistency_desc <- ind_consistency_desc %>%
  mutate(direction = case_when(mean_dev >  dev_eps ~ "Above age-specific mean",
                               mean_dev < -dev_eps ~ "Below age-specific mean",
                               TRUE ~ "Near age-specific mean"),
         consistency_rank = percent_rank(-sd_resid),
         consistency_class = case_when(consistency_rank >= 0.75 ~ "More consistent",
                                       consistency_rank <= 0.25 ~ "Less consistent",
                                       TRUE ~ "Intermediate"),
         consistency_plot = -sd_resid)

plot_individual_variation <- ggplot(ind_consistency_desc,
                                    aes(x = mean_dev, y = sd_resid, color = direction)) +
  geom_point(aes(size = total_resights_life), alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    x = "Mean deviation from model-based age-specific expected MOA",
    y = "Within-individual SD of Pearson residuals",
    color = "Direction",
    size = "Lifetime observations") +
  theme_few(base_size = 14) +
  theme(legend.position = "top"); plot_individual_variation

### strategy plot figure ###

cons_thresh <- median(ind_consistency_desc$sd_resid, na.rm = TRUE)

# create individual strategies based on whether the animal is above or below the median?
ind_consistency_desc <- ind_consistency_desc %>%
  mutate(consistency_level = if_else(sd_resid <= cons_thresh,
                                     "More consistent",
                                     "Less consistent"),
         
         strategy = case_when(
           direction == "Above age-specific mean" & consistency_level == "More consistent" ~ "Low variation, above",
           direction == "Below age-specific mean" & consistency_level == "More consistent" ~ "Low variation, below",
           direction == "Near age-specific mean"  & consistency_level == "More consistent" ~ "Low variation, average",
           
           direction == "Above age-specific mean" & consistency_level == "Less consistent" ~ "High variation, above",
           direction == "Below age-specific mean" & consistency_level == "Less consistent" ~ "High variation, below",
           direction == "Near age-specific mean"  & consistency_level == "Less consistent" ~ "High variation, average"
         )
  )

plot_maternal_strategy <- ggplot(ind_consistency_desc,
                                 aes(x = mean_dev, y = sd_resid, color = strategy)) +
  geom_point(aes(size = total_resights_life), alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = cons_thresh, linetype = "dashed", color = "grey50") +
  labs(
    x = "Mean deviation from age-specific expected MOA",
    y = "Within-individual variability (SD of Pearson residuals)",
    color = "Strategy",
    size = "Lifetime observations",
    title = "Individual maternal strategies relative to age-specific expectation") +
  theme_few(base_size = 14) +
  theme(legend.position = "top"); plot_maternal_strategy

########## extract population-level repeatability from the model ###############

vc <- as.data.frame(VarCorr(mod_binom_1996_2025))

var_ind  <- vc[vc$grp == "animalID_fct", "vcov"]
var_seas <- vc[vc$grp == "season_fct", "vcov"]

# For binomial models, residual variance is fixed at:
var_resid <- (pi^2) / 3

R <- var_ind / (var_ind + var_seas + var_resid); R

get_repeatability <- function(mod, id_var = "animalID_fct", season_var = "season_fct") {
  vc <- as.data.frame(VarCorr(mod))
  
  var_id <- vc$vcov[vc$grp == id_var]
  var_season <- vc$vcov[vc$grp == season_var]
  
  if (length(var_id) == 0) var_id <- 0
  if (length(var_season) == 0) var_season <- 0
  
  var_resid <- (pi^2) / 3
  
  R <- var_id / (var_id + var_season + var_resid)
  return(R)
}

R_obs <- get_repeatability(mod_binom_1996_2025)
R_obs

permute_repeatability <- function(dat, n_perm = 1000, seed = 123) {
  set.seed(seed)
  
  R_null <- rep(NA_real_, n_perm)
  failed <- 0
  
  for (i in seq_len(n_perm)) {
    dat_perm <- dat %>%
      mutate(animalID_fct = sample(animalID_fct, replace = FALSE))
    
    mod_perm <- try(
      glmer(
        proportion ~ age_cat:age10 + (1 | animalID_fct) + (1 | season_fct),
        weights = total_resights,
        family = binomial(link = "logit"),
        control = glmerControl(optimizer = "bobyqa"),
        data = dat_perm
      ),
      silent = TRUE
    )
    
    if (inherits(mod_perm, "try-error")) {
      failed <- failed + 1
      next
    }
    
    R_null[i] <- get_repeatability(mod_perm)
  }
  
  message("Failed fits: ", failed)
  R_null
}

R_null <- permute_repeatability(intrinsic_variables, n_perm = 100, seed = 1)

summary(R_null)
quantile(R_null, c(0.025, 0.5, 0.975), na.rm = TRUE)

p_emp <- mean(R_null >= R_obs, na.rm = TRUE)
p_emp

hist(R_null,
     breaks = 30,
     col = "grey80",
     border = "white",
     main = "Null distribution of repeatability",
     xlab = "Repeatability (R)")

abline(v = R_obs, col = "red", lwd = 2)
abline(v = null_ci, col = "blue", lty = 2, lwd = 2)

mod_predict <- glmmTMB(
  proportion ~ age_cat:age10 + (1 | animalID_fct) + (1 | season_fct),
  weights = total_resights,
  dispformula = ~ (1 | animalID_fct),
  family = binomial(link = "logit"),
  data = intrinsic_dev_3
)
summary(mod_predict)
VarCorr(mod_predict)
summary(mod_predict)$varcor
sigma(mod_predict)

mod_predict_null <- glmmTMB(
  proportion ~ age_cat:age10 + (1 | animalID_fct) + (1 | season_fct),
  weights = total_resights,
  dispformula = ~ 1,
  family = binomial(link = "logit"),
  data = intrinsic_dev_3
); mod_predict_null

anova(mod_predict_null, mod_predict)