##Modeling intrinsic variables!

##General distribution of proportions
ggplot(data = Proportion_MPA, aes(x = proportion)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Proportion MPA",
       x = "Proportion MPA", y = "Frequency") +
  scale_y_continuous(n.breaks = 10)
  theme_few()

##Graph age against proportion
ggplot(data = intrinsic_variables, aes(x = AgeYears, y = proportion)) +
  geom_point() +
  labs(title = "Proportion MPA by Age", x = "Age", y = "MPA") +
  theme_few()

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

##Zero-One Inflated Beta Regression model for age and proportion
install.packages("gamlss")
library(gamlss)
install.packages("gamlss.dist")
library(gamlss.dist)
packageVersion("gamlss")

##Model with age and year born as independent effects
mod_01_BR_intrinsic <- gamlss(proportion ~ AgeYears + year_born + random(as.factor(animalID)) + random(as.factor(season)), family = "BEINF", data = intrinsic_variables) #the model
summary(mod_01_BR_Age_YearBorn) #summary statistics
exp(coef(mod_01_BR_Age_YearBorn)) #test for covariance
wp(mod_01_BR_intrinsic, ylim.all = 4) #worm plot
wp(mod_01_BR_intrinsic, xvar = intrinsic_variables$AgeYears)

#Another model with age and year born, but made AgeYears smooth/non-linear
mod_02_BR_intrinsic_smooth <- gamlss(proportion ~ pb(AgeYears) + year_born + random(as.factor(animalID)) + random(as.factor(season)),
  sigma.formula = ~ pb(AgeYears),
  nu.formula    = ~ pb(AgeYears),
  tau.formula   = ~ pb(AgeYears),
  family = BEINF,
  data = intrinsic_variables)
summary(mod_02_BR_intrinsic_smooth) #summary stats
exp(coef(mod_02_BR_intrinsic_smooth)) #test for covariance
wp(mod_02_BR_intrinsic_smooth, ylim.all = 4) #worm plot
wp(mod_02_BR_intrinsic_smooth, xvar = intrinsic_variables$AgeYears, ylim.worm = 4) #worm plot for Age

##try this for assessing collinearity
cor(intrinsic_variables$AgeYears, intrinsic_variables$year_born)

#try cohorts to see if it fixes adding cohort to model
intrinsic_variables$cohort <- cut(
  intrinsic_variables$year_born,
  breaks = c(1980, 1990, 2000, 2010, 2015, 2022),
  labels = c("1981–89", "1990–99", "2000–09", "2010–14", "2015–21"),
  right = FALSE)

##Cohorts work yay
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

#remove the random 0 proportion in data
intrinsic_variables <- intrinsic_variables %>%
  filter(proportion > 0)

##GLMM for age, year_born and proportion
install.packages("glmmTMB")
library(glmmTMB)
library(splines)

##new type of model, use GLMM for Gamma distribution with poly for Age
mod_gamma <- glmmTMB(
  proportion ~ poly(AgeYears, 3) + as.factor(year_born) + (1 | animalID) + (1 | season),
  data = intrinsic_variables,
  family = Gamma(link = "log")
)
summary(mod_gamma)

##calculate age last observed and add it to intrinsic_variables
intrinsic_variables <- intrinsic_variables %>%
  group_by(animalID) %>%
  mutate(age_last_seen = max(AgeYears, na.rm = TRUE)) %>%
  ungroup()

##Now try GLMM with splines and add age last seen
mod_spline <- glmmTMB(
  proportion ~ ns(AgeYears, df = 3) + as.factor(year_born) + age_last_seen + (1 | animalID) + (1 | season),
  data = intrinsic_variables,
  family = Gamma(link = "log")
)
summary(mod_spline)

#have year_born as a factor in the dataframe
intrinsic_variables$year_born <- as.factor(intrinsic_variables$year_born)

##using simulated data for predicted values
median_age_last <- median(intrinsic_variables$age_last_seen, na.rm = TRUE)

new_intrinsic_data <- data.frame(
  AgeYears = seq(min(intrinsic_variables$AgeYears, na.rm = TRUE),
                 max(intrinsic_variables$AgeYears, na.rm = TRUE),
                 length.out = 100),
  year_born = factor("2012", levels = levels(intrinsic_variables$year_born)),
  age_last_seen = median_age_last,
  animalID = NA,
  season = NA
)

new_intrinsic_data$pred <- predict(mod_spline, newdata = new_intrinsic_data, type = "response")

ggplot(new_intrinsic_data, aes(x = AgeYears, y = pred)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Predicted Proportion vs. Age (Year Born = 2012)",
    x = "Age (Years)",
    y = "Predicted Proportion"
  ) +
  theme_minimal()

#ggpredict version
install.packages("ggeffects") 
library(ggeffects)

age_effect <- ggpredict(mod_spline, terms = "AgeYears [all]")
plot(age_effect) +
  ggtitle("Effect of AgeYears on Predicted Proportion") +
  ylab("Predicted Proportion") +
  xlab("Age (Years)")


##Graph year against proportion
##violin plot for distribution of MPA across all seasons
ggplot(data = Proportion_MPA, aes(x = factor(season), y = proportion)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.5) + 
  labs(title = "MPA distribution by season", x = "Season", y = "Proportion")

##model using a binomial to 





 

