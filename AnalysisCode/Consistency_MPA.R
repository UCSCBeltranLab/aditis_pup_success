##Consistency models and figures
library(performance)

###################### Playing around with consistency ##########################

#plot distribution of standard deviations
ggplot(data = consistency_filtered, aes(x = sd_proportion)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  labs(title = "Histogram of std dev",
       x = "SD", y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  theme_minimal()

##calculating consistency metrics
consistency_by_individual <- intrinsic_variables %>%
  group_by(animalID) %>%
  summarise(mean_proportion = mean(proportion, na.rm = TRUE),
            sd_proportion = sd(proportion, na.rm = TRUE),
            n_obs = n()) %>%
  arrange(sd_proportion)

#filter by greater than 3 calculations across lifetime
consistency_filtered <- consistency_by_individual %>%
  filter(n_obs >= 3)

##plot for consistency- standard dev per individual
ggplot(data = consistency_filtered, aes(x = reorder(animalID, sd_proportion), y = sd_proportion)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Within-Individual Variability in Proportion", x = "Animal ID (sorted)", y = "SD of Proportion")

##plot avg prop by std deviation per individual
ggplot(data = consistency_filtered, aes(x = mean_proportion, y = sd_proportion)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Average Proportion", y = "Standard Deviation",
       title = "Individual Consistency in Mom-Pup Association")

##filter intrinsic variables for >= 3 observations
intrinsic_variables_filt <- intrinsic_variables %>%
  filter(animalID %in% consistency_filtered$animalID)

n_distinct(intrinsic_variables_filt$animalID) #check that the resulting animals match previous table

##create a combined consistency data table
consistency_data <- intrinsic_variables_filt %>%
  left_join(consistency_filtered, by = "animalID")

##make a data table to create a star label for animals with consistent prop == 1
star_data <- consistency_data %>%
  group_by(animalID_fct) %>%
  summarise(all_one = all(proportion == 1),
            max_y = max(proportion)) %>%
  filter(all_one) %>%
  mutate(y_position = max_y + 0.05)

##boxplot with 2 color gradient for consistency
ggplot(consistency_data, aes(x = animalID_fct, y = proportion, fill = sd_proportion)) +
  geom_boxplot(width = 0.7, outlier.size = 1.5, outlier.color = "black") +
  geom_jitter(size = 1, alpha = 0.6, color = "#04BBB2") +
  scale_fill_gradient(low = "#D8FFF7", high = "#073481", name = "SD of Proportion") +
  geom_text(data = star_data, aes(x = animalID_fct, y = y_position, label = "*"), color = "#C699E1", size = 7, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 1.03)) +
  scale_y_continuous(n.breaks = 10) +
  theme_few() +
  labs(x = "Animal ID", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, size = 6))

##violin plot with 2 color gradient for consistency
ggplot(consistency_data, aes(x = animalID_fct, y = proportion, fill = sd_proportion)) +
  geom_violin(scale = "width", trim = TRUE, color = "black") +  # violin instead of boxplot
  geom_jitter(size = 1, alpha = 0.6, color = "#04BBB2", width = 0.2) +
  scale_fill_gradient(low = "#D8FFF7", high = "#073481", name = "SD of Proportion") +
  geom_text(data = star_data, aes(x = animalID_fct, y = y_position, label = "*"), 
            color = "#C699E1", size = 7, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 1.03)) +
  scale_y_continuous(n.breaks = 10) +
  theme_classic() +
  labs(x = "Animal ID", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, size = 6))

###################### Consistency Model using SD  ################################
##creating a consistency model using animal SD as the response
#consistency as a function of covariates
consistency_data <- consistency_data %>%
  mutate(age_last_seen = as.numeric(age_last_seen)) %>%
  mutate(n_obs = as.numeric(n_obs))

#now let's find out how many pups each mom had across her lifetime
total_pups_lifetime <- metadata %>%
  select(animalID, season) %>%
  distinct() %>%
  group_by(animalID) %>%
  summarise(total_pups = n(), .groups = "drop")

consistency_data <- consistency_data %>%
  left_join(total_pups_lifetime, by = "animalID")

mod_gamma_sd <- glmmTMB(sd_proportion ~ age_last_seen + n_obs,
                 data = consistency_data[consistency_data$sd_proportion > 0, ],
                 family = Gamma(link = "log"))
summary(mod_gamma_sd)
exp(fixef(mod_gamma_sd)$cond) 
DHARMa::simulateResiduals(mod_gamma_sd, plot = TRUE)
check_collinearity(mod_gamma_sd) #checking for collinearity

################################ Consistency MPA DHGLM ######################################

library(brms)
library("bayestestR")
install.packages("tidybayes")
library(tidybayes)

#trying brms version with 0-1 inflated beta
mod_consistency <- brm(
  bf(proportion ~ AgeYears + total_resights + (1 | animalID_fct) + (1 | season_fct), #mean structure
     phi ~ AgeYears + (1 | animalID_fct), #residual variation = consistency
     zoi ~ AgeYears + total_resights),  #zero-one inflation
  family = zero_one_inflated_beta(),
  data = consistency_data,
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.999))
summary(mod_consistency)
pp_check(mod_consistency)

check_collinearity(mod_consistency)

############### Calculate population repeatability

# Calculate mean proportion
mu <- mean(consistency_data$proportion, na.rm = TRUE)

# Extract posterior draws from your brms model
draws <- as_draws_df(mod_consistency)

# Calculate mean phi on original scale (precision parameter)
mean_phi_log <- mean(draws$b_phi_Intercept)
phi <- exp(mean_phi_log)

# Variance components
var_animal_logit <- 0.84^2  # extracted SD^2 from model output

# Transform between-animal variance to response scale using delta method
# Derivative of inverse logit at mean proportion
d_invlogit <- mu * (1 - mu)

# Transformed between-animal variance
var_animal_resp <- (d_invlogit^2) * var_animal_logit

# Residual variance on response scale
var_resid <- (mu * (1 - mu)) / (1 + phi)

# Repeatability calculation
R <- var_animal_resp / (var_animal_resp + var_resid)

# Round and print
round(R, 3)

######### Calculate within-individual Repeatability

# Step 1: Convert draws to data frame
draws <- as_draws_df(mod_consistency)

# Step 2: Identify phi-related random effect columns
phi_cols <- grep("^r_animalID_fct__phi\\[.*?,Intercept\\]$", names(draws), value = TRUE)

# Step 3: Compute posterior means for each animal’s phi RE
re_phi_df <- draws[, phi_cols] |> 
  summarise(across(everything(), mean, na.rm = TRUE)) |> 
  t() |> 
  as.data.frame()

# Step 4: Clean up animal IDs from rownames
re_phi_df$animalID <- gsub("r_animalID_fct__phi\\[|,Intercept\\]", "", rownames(re_phi_df))
colnames(re_phi_df)[1] <- "re_phi"

# Step 5: Get fixed effect for phi (log scale)
fix_phi <- mean(draws$Intercept_phi, na.rm = TRUE)

# Step 6: Compute phi_i for each animal
re_phi_df$phi_i <- exp(fix_phi + re_phi_df$re_phi)

# Step 7: Compute residual and repeatability
mu_hat <- mean(consistency_data$proportion, na.rm = TRUE)
re_phi_df$var_resid_i <- (mu_hat * (1 - mu_hat)) / (1 + re_phi_df$phi_i)
re_phi_df$R_i <- var_animal_resp / (var_animal_resp + re_phi_df$var_resid_i)

hist(re_phi_df$R_i, breaks = 30, main = "Individual Repeatability (R_i)", xlab = "R_i")
  



