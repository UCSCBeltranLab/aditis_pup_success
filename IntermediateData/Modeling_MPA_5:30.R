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

##Poisson model for age and proportion - a little better
mod_Poisson_Age_by_Proportion <- glm(proportion ~ AgeYears, data = metadata, family = "poisson")
summary(mod_Poisson_Age_by_Proportion)

library(MASS)
##Negative Binomial for age and proportion
mod_NB_Age_by_Proportion <- glm.nb(proportion ~ AgeYears, data = metadata)
summary(mod_NB_Age_by_Proportion)

##Zero-One Inflated Beta Regression model for age and proportion
install.packages("gamlss")
library(gamlss)

##Model with age and year born as independent effects
mod_01_BR_Age_YearBorn <- gamlss(proportion ~ AgeYears + year_born, family = "BEINF", data = na.omit(metadata))
summary(mod_01_BR_Age_YearBorn)
exp(coef(mod_01_BR_Age_YearBorn))

##Model with only age
mod_01_BR_YearBorn <- gamlss(proportion ~ year_born, family = "BEINF", data = na.omit(metadata))
summary(mod_01_BR_YearBorn)
##find out how to control for animalID

##Graph year against proportion
##violin plot for distribution of MPA across all seasons
ggplot(data = Proportion_MPA, aes(x = factor(season), y = proportion)) +
  geom_violin(fill = "lightblue", color = "darkblue") +
  geom_jitter(width = 0.1, alpha = 0.5) + 
  labs(title = "MPA distribution by season", x = "Season", y = "Proportion")



