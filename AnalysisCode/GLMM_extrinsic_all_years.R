##Let's make the extrinsic model for all years!

##First we need the wave and tide datasets together with the proportion data!
extrinsic_weather <- metadata %>%
  select(animalID, season, proportion, total_resights, AgeYears) %>%
  distinct() %>%
  filter(!is.na(proportion), proportion > 0) %>%
  left_join(harem_assignment, by = c("animalID", "season")) %>%
  left_join(tide_wave_flagged, by = "season") %>%
  mutate(is_one = as.numeric(proportion == 1)) %>%
  mutate(season_fct = as.factor(season),
         animalID_fct = as.factor(animalID),
         area_fct = as.factor(area))

extrinsic_weather <- extrinsic_weather %>%
  mutate(group = case_when(
    grepl("^N", area_fct) | area_fct %in% c("BBNN", "BBNS", "BBN") ~ "NP",
    TRUE ~ "SP"))

##Make the subset for extrinsic weather gamma model
##create a version of intrinsic_variables with proportion < 1
extrinsic_weather_sub <- extrinsic_weather %>%
  filter(proportion < 1)

##transform proportion for this filtered data to make it suitable for gamma distribution 
flipped_prop <- max(extrinsic_weather_sub$proportion) - extrinsic_weather_sub$proportion #subtract all proportions from the maximum
flipped_prop <- flipped_prop - min(flipped_prop) + 0.0000001 #ensure no 0s by adding small constant
extrinsic_weather_sub$flipped_prop <- flipped_prop

##create a version of transformed flipped prop for Gaussian model
extrinsic_weather_sub <- extrinsic_weather_sub %>%
  mutate(prop_trans = car::logit(proportion, adjust = 0.001))


#################also try including area as a fixed effect by NP or SP areas#################

##binomial
mod_abiotic_bin <- glmmTMB(is_one ~ n_extreme_both + total_resights + (1 | season_fct) + (1 | animalID_fct) + (1 | area_fct),
                            data = extrinsic_weather,
                            family = binomial(link = "logit"))
summary(mod_abiotic_bin)
exp(fixef(mod_abiotic_bin)$cond) 
DHARMa::simulateResiduals(mod_abiotic_bin, plot = TRUE)
DHARMa::testDispersion(DHARMa::simulateResiduals(mod_abiotic_bin))
check_collinearity(mod_abiotic_bin)

##gamma
mod_abiotic_gam <- glmmTMB(flipped_prop ~ n_extreme_both + total_resights + (1 | season_fct) + (1 | animalID_fct) + (1 | area_fct),
                            data = extrinsic_weather_sub,
                            family = Gamma(link = "log"))
summary(mod_abiotic_gam)
exp(fixef(mod_abiotic_gam)$cond) 
DHARMa::simulateResiduals(mod_abiotic_gam, plot = TRUE)
DHARMa::testDispersion(DHARMa::simulateResiduals(mod_abiotic_gam))
check_collinearity(mod_abiotic_gam)

##gaussian
mod_gauss_abiotic <- glmmTMB(prop_trans ~ n_extreme_both + total_resights + (1 | animalID_fct) + (1 | season_fct) + (1 | area_fct),
                               family = gaussian(),
                               data = extrinsic_weather_sub)
summary(mod_gauss_abiotic)
exp(fixef(mod_gauss_abiotic$cond))
DHARMa::simulateResiduals(mod_gauss_abiotic, plot = TRUE)
DHARMa::testDispersion(DHARMa::simulateResiduals(mod_gauss_abiotic))
check_collinearity(mod_gauss_abiotic)
