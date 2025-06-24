##Consistency models and figures

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

##plot with 2 color gradient for consistency
ggplot(consistency_data, aes(x = animalID_fct, y = proportion, fill = sd_proportion)) +
  geom_boxplot(coef = 1, width = 0.7, outlier.size = 1.5, outlier.color = "black") +
  geom_jitter(size = 1, alpha = 0.6, color = "#04BBB2") +
  scale_fill_gradient(low = "#D8FFF7", high = "#073481", name = "SD of Proportion") +
  geom_text(data = star_data, aes(x = animalID_fct, y = y_position, label = "*"), color = "#C699E1", size = 7, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 1.03)) +
  scale_y_continuous(n.breaks = 10) +
  theme_classic() +
  labs(x = "Animal ID", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, size = 6))

##creating a consistency model using animal SD as the response
#consistency as a function of other variables
mod_gamma <- glmmTMB(sd_proportion ~ mean_proportion + age_last_seen + n_obs,
                 data = consistency_data[consistency_data$sd_proportion > 0, ],
                 family = Gamma(link = "log"))
summary(mod_gamma)
DHARMa::simulateResiduals(mod_gamma, plot = TRUE)

##Modeling consistency using a DHGLM?
#This only works for values between 0 and 1 though
consistency_data_sub <- consistency_data %>%
  filter(proportion < 1)

mod_dhglm <- glmmTMB(proportion ~ AgeYears + (1 | animalID_fct),
  dispformula = ~ 1 + (1 | animalID_fct),
  family = beta_family(link = "logit"),
  data = consistency_data_sub)
summary(mod_dhglm)
DHARMa::simulateResiduals(mod_dhglm, plot = TRUE)


##Trying a gamma with flipped proportion?
flipped_prop <- max(consistency_data$proportion) - consistency_data$proportion
flipped_prop <- flipped_prop - min(flipped_prop) + 0.001
consistency_data$flipped_prop <- flipped_prop

mod_dhglm_gamma <- glmmTMB(
  flipped_prop ~ AgeYears + (1 | animalID),
  dispformula = ~ 1 + (1 | animalID),
  family = Gamma(link = "log"),
  data = consistency_data
)
summary(mod_dhglm_gamma)
DHARMa::simulateResiduals(mod_dhglm_gamma, plot = TRUE)




