######## Overall distribution of proportions #########

ggplot(data = intrinsic_variables, aes(x = proportion)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "darkblue") +
  labs(x = "Mother-offspring association", 
       y = "Frequency") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme_few()


############ facet the two age and experience figures together ###########

plot_age <- cowplot::plot_grid(plot_age_1996_2025 + 
                                 theme(axis.title.x = element_blank(),
                                       axis.text.x  = element_blank(),
                                       axis.ticks.x = element_blank()),
                               plot_age_2016_2023,
                               ncol = 1,
                               align = "v",
                               axis = "lr",
                               rel_heights = c(1, 1),
                               labels = c("(a)", "(b)"),
                               label_size = 14,
                               label_fontface = "bold",
                               label_colour = "black",
                               label_x = 0.02,
                               label_y = 0.98,
                               hjust = 0); plot_age

ggsave("./TablesFigures/age_figure.png", plot = plot_age, width = 8, height = 10, dpi = 600)

plot_experience <- cowplot::plot_grid(plot_exp_1996_2025 +
                                        theme(axis.title.x = element_blank(),
                                              axis.text.x = element_blank(),
                                              axis.ticks.x = element_blank()),
                                      plot_exp_2016_2023,
                                      ncol = 1,
                                      align = "v",
                                      axis = "lr",
                                      rel_heights = c(1, 1),
                                      labels = c("(a)", "(b)"),
                                      label_size = 14,
                                      label_fontface = "bold",
                                      label_colour = "black",
                                      label_x = 0.02,
                                      label_y = 0.98,
                                      hjust = 0); plot_experience

ggsave("./TablesFigures/experience_figure.png", plot = plot_experience, width = 8, height = 10, dpi = 600)

############ conceptual MOA figure #########

example_ids <- c("36358", "48257", "51872")

plot_df <- metadata %>%
  mutate(
    date = as.Date(date),
    BirthDate = as.Date(BirthDate),
    animalID = as.character(animalID),
    withpup = as.character(withpup)
  ) %>%
  filter(season == 2024, animalID %in% example_ids) %>%
  transmute(
    animalID, season, AgeYears, proportion,
    date,
    pup_status = case_when(
      withpup == "0" ~ "0",
      withpup == "1" ~ "1",
      TRUE ~ NA_character_
    ),
    label = withpup
  ) %>%
  filter(!is.na(date), !is.na(pup_status)) %>%
  bind_rows(
    metadata %>%
      mutate(
        BirthDate = as.Date(BirthDate),
        animalID = as.character(animalID)
      ) %>%
      filter(season == 2024, animalID %in% example_ids, !is.na(BirthDate)) %>%
      transmute(
        animalID, season, AgeYears, proportion,
        date = BirthDate,
        pup_status = "birth",
        label = NA_character_)
  ) %>%
  arrange(animalID, date, desc(pup_status == "birth")) %>%
  distinct(animalID, season, date, .keep_all = TRUE) %>%
  group_by(animalID, season) %>%
  mutate(
    panel_lab = paste0(
      animalID, " (Age ", first(AgeYears),
      ", MOA=", sprintf("%.2f", first(proportion)), ")")) %>%
  ungroup()

plot_df <- plot_df %>%
  group_by(animalID, season) %>%
  mutate(panel_lab = paste0(
    "<b>", animalID, "</b>",
    "<br>Age ", first(AgeYears),
    "<br>MOA = ", sprintf("%.2f", first(proportion))
  )) %>%
  ungroup()

all_dates <- seq(min(plot_df$date), max(plot_df$date), by = "day")

plot_df_full <- plot_df %>%
  mutate(label = if_else(pup_status == "birth", "1", label)) %>%
  select(panel_lab, date, pup_status, label) %>%
  group_by(panel_lab) %>%
  complete(date = all_dates) %>%
  ungroup()

conceptual_MOA_figure <- ggplot(plot_df_full, aes(date, 1, fill = pup_status)) +
  geom_tile(
    width = 0.95, height = 0.95,
    color = "white", linewidth = 0.3
  ) +
  geom_text(
    data = subset(plot_df_full, !is.na(label)),
    aes(label = label),
    size = 4.5,
    fontface = "bold",
    color = "black"
  ) +
  facet_grid(panel_lab ~ ., switch = "y") +
  scale_fill_manual(
    name = "Pup status",
    values = c(
      "0" = "#86c5e3",
      "1" = "#7bc96f",
      "birth" = "#f2b6c6"
    ),
    breaks = c("0", "1", "birth"),
    labels = c("0" = "0", "1" = "1", "birth" = "Approximate birth"),
    na.value = "white"
  ) +
  scale_x_date(breaks = seq(min(plot_df_full$date), max(plot_df_full$date), by = "3 days"),
               date_labels = "%b %d",
               expand = c(0, 0)) +
  scale_y_continuous(breaks = NULL) +
  labs(x = "Date", y = NULL) +
  theme_classic(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    strip.text.y.left = element_markdown(
      angle = 0,
      hjust = 0,
      size = 15,
      lineheight = 1.3,
    ),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(color = "black",),
    axis.text.x = element_text(
      angle = 70,
      hjust = 1,
      size = 13
    ),
    panel.spacing.y = unit(1.2, "lines"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 14)
  ); conceptual_MOA_figure

ggsave("./TablesFigures/conceptual_MOA_figure.png", conceptual_MOA_figure, width = 14, height = 8)

################ density conceptual plot ###################
library(sf)
library(dplyr)
library(lubridate)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(viridis)
library(maptiles)
library(tidyterra)

# 1) Read beaches shapefile
beaches <- st_read("./RawData/ANM map/Ano Nuevo Map Final.shp", quiet = TRUE) %>%
  select(-id)

# 2) Read density data
seal.density <- read.csv("./RawData/seal.density.csv")

# 3) Mean density per beach across years
area_density_mean <- seal.density %>%
  rename(dominant_area = Beach) %>%
  group_by(dominant_area) %>%
  summarize(mean_density = mean(density, na.rm = TRUE), .groups = "drop")

# 4) Check the beach-name field in your shapefile
names(beaches)

# Replace "Beach" below if your shapefile uses a different name field
beaches_density <- beaches %>%
  left_join(area_density_mean, by = c("Beach" = "dominant_area"))

# 5) Make labels: beach name + mean density
label_points <- st_point_on_surface(beaches_density) %>%
  mutate(
    label_text = paste0(
      Beach, "\n",
      round(mean_density, 2)
    )
  )

# 6) Get extent for zooming and tile download
bbox <- st_bbox(beaches_density)

# Convert bbox to sf object for maptiles
bbox_sf <- st_as_sfc(bbox) %>%
  st_as_sf()

# expand bbox a little before downloading tiles
buffer_x <- 0.0007
buffer_y <- 0.0007

bbox_expanded <- bbox
bbox_expanded["xmin"] <- bbox["xmin"] - buffer_x
bbox_expanded["xmax"] <- bbox["xmax"] + buffer_x
bbox_expanded["ymin"] <- bbox["ymin"] - buffer_y
bbox_expanded["ymax"] <- bbox["ymax"] + buffer_y

bbox_sf_expanded <- st_as_sfc(bbox_expanded) %>%
  st_as_sf()

# re-download tiles using expanded extent
sat_map <- get_tiles(
  x = bbox_sf_expanded,
  provider = "Esri.WorldImagery",
  crop = TRUE,
  zoom = 16
)

seal_density_map <- ggplot() +
  geom_spatraster_rgb(data = sat_map) +
  geom_sf(
    data = beaches_density,
    aes(fill = mean_density),
    color = "white",
    linewidth = 0.4,
    alpha = 0.8
  ) +
  scale_fill_viridis_c(
    option = "mako",
    direction = -1,
    name = "Conspecific density
(seals per 10m radius)"
  ) +
  coord_sf(
    xlim = c(bbox_expanded["xmin"], bbox_expanded["xmax"]),
    ylim = c(bbox_expanded["ymin"], bbox_expanded["ymax"]),
    expand = FALSE
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ); seal_density_map

ggsave("./TablesFigures/seal_density_map.png", plot = seal_density_map, width = 12, height = 10, dpi = 300, bg = "white")

################### model output tables ###################

library(broom.mixed)
library(dplyr)
library(flextable)
library(officer)

make_model_table <- function(model, fixed_labels = NULL, random_labels = NULL) {
  
  fixef_tbl <- broom.mixed::tidy(model, effects = "fixed") %>%
    mutate(Predictor = if (is.null(fixed_labels)) term else ifelse(term %in% names(fixed_labels), fixed_labels[term], term)) %>%
    transmute(
      Predictor,
      Estimate = round(estimate, 2),
      SE = round(std.error, 2),
      Z = round(statistic, 2),
      `P-value` = case_when(
        is.na(p.value) ~ NA_character_,
        p.value < 0.001 ~ paste0(format(p.value, digits = 2, scientific = TRUE), " ***"),
        p.value < 0.01 ~ paste0(format(round(p.value, 4), nsmall = 4), " **"),
        p.value < 0.05 ~ paste0(format(round(p.value, 4), nsmall = 4), " *"),
        TRUE ~ format(round(p.value, 4), nsmall = 4)
      )
    )
  
  ranef_tbl <- as.data.frame(VarCorr(model)) %>%
    mutate(Predictor = if (is.null(random_labels)) grp else ifelse(grp %in% names(random_labels), random_labels[grp], grp)) %>%
    transmute(
      Predictor = paste0("Random effect: ", Predictor),
      Estimate = round(vcov, 3),
      SE = NA,
      Z = round(sdcor, 3),
      `P-value` = NA
    )
  
  bind_rows(fixef_tbl, ranef_tbl) %>%
    flextable() %>%
    align(align = "center", part = "all") %>%
    autofit()
}

mod_table_age_2016_2023 <- make_model_table(
  mod_binom_2016_2023,
  fixed_labels = c(
    "(Intercept)" = "Intercept",
    "AgeYears" = "Maternal age",
    "avg_density" = "Average conspecific density",
    "n_extreme_both" = "Number of extreme wave and tide events"
  ),
  random_labels = c(
    "animalID_fct" = "individualID",
    "season_fct" = "year"
  )
); mod_table_age_2016_2023

#save final table
save_as_docx(mod_table_age_2016_2023, path = "./TablesFigures/Model_Output_Age_2016_2023.docx")

mod_table_exp_2016_2023 <- make_model_table(
  mod_binom_exp_2016_2023,
  fixed_labels = c(
    "(Intercept)" = "Intercept",
    "experience_prior" = "Prior pupping experience",
    "avg_density" = "Average conspecific density",
    "n_extreme_both" = "Number of extreme wave and tide events"
  ),
  random_labels = c(
    "animalID_fct" = "individualID",
    "season_fct" = "year"
  )
); mod_table_exp_2016_2023

#save final table
save_as_docx(mod_table_exp_2016_2023, path = "./TablesFigures/Model_Output_Experience_2016_2023.docx")

mod_table_age_1996_2025 <- make_model_table(
  mod_binom_1996_2025,
  fixed_labels = c(
    "(Intercept)" = "Intercept",
    "age_catYoung:age10" = "Age : Pre-threshold ( < 9 years)",
    "age_catOld:age10" = "Age : Post-threshold ( ≥ 9 years)"
  ),
  random_labels = c(
    "animalID_fct" = "individual",
    "season_fct" = "year"
  )
); mod_table_age_1996_2025

#save final table
save_as_docx(mod_table_age_1996_2025, path = "./TablesFigures/Model_Output_Age_1996_2025.docx")


mod_table_exp_1996_2025 <- make_model_table(
  mod_exp_1996_2025,
  fixed_labels = c(
    "(Intercept)" = "Intercept",
    "experience_catinexperienced:exp10" = "Prior pupping experience : Pre-threshold ( < 5 previous pups)",
    "experience_catexperienced:exp10" = "Prior pupping experience : Post-threshold ( ≥ 5 previous pups)"
  ),
  random_labels = c(
    "animalID_fct" = "individual",
    "season_fct" = "year"
  )
); mod_table_exp_1996_2025

#save final table
save_as_docx(mod_table_exp_1996_2025, path = "./TablesFigures/Model_Output_Experience_1996_2025.docx")


mod_wean_mass_table <- make_model_table(
  mod_wean_age,
  fixed_labels = c(
    "(Intercept)" = "Intercept",
    "AgeYears" = "Maternal age",
    "proportion" = "Mother-offspring association"
  ),
  random_labels = c(
    "animalID_fct" = "individual",
    "season_fct" = "year"
  )
); mod_wean_mass_table

save_as_docx(mod_wean_mass_table, path = "./TablesFigures/Model_Output_Wean_Mass.docx")




