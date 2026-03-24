############ facet the two age figures together ###########

plot_age <- plot_grid(
  plot_age_1996_2025 + 
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    ),
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
  hjust = 0,
  vjust = 1
); plot_age

ggsave("./TablesFigures/experience_figure.png", plot = plot_exp_1996_2025, width = 10, height = 6, dpi = 600)

############ conceptual MOA figure #########

# 1) make sure dates are as Date and animalIDs are characters
metadata <- metadata %>%
  mutate(date = as.Date(date),
         BirthDate = as.Date(BirthDate))

# 2) Choose 3 animals from 2025
example_ids <- c("44397","48257", "51872") ##select animals

meta_subset <- metadata %>%
  filter(season == 2024, animalID %in% example_ids)

# 3) Resight tiles
resight_tiles <- meta_subset %>%
  filter(!is.na(date), !is.na(withpup)) %>%
  distinct(animalID, season, date, .keep_all = TRUE) %>%
  mutate(withpup = as.character(withpup),
         pup_status = case_when(withpup == "0" ~ "0",
                                withpup == "1" ~ "1",
                                TRUE ~ NA_character_)) %>%
  select(animalID, season, date, AgeYears, proportion, pup_status, withpup, area)

# 4) Approximate birth tiles
birth_tiles <- meta_subset %>%
  distinct(animalID, season, BirthDate, AgeYears, proportion, area, .keep_all = TRUE) %>%
  filter(!is.na(BirthDate)) %>%
  transmute(
    animalID,
    season,
    date = BirthDate,
    AgeYears,
    proportion,
    area,
    pup_status = "approximate birth",
    withpup = "1"
  )

# 5) Combine tiles
# keep approximate birth if same day as a resight
plot_df <- bind_rows(resight_tiles, birth_tiles) %>%
  mutate(
    pup_status = factor(pup_status, levels = c("approximate birth", "0", "1"))
  ) %>%
  arrange(animalID, season, date, pup_status) %>%
  distinct(animalID, season, date, .keep_all = TRUE)

# 6) Panel labels
label_df <- plot_df %>%
  group_by(animalID, season) %>%
  summarise(
    AgeYears = first(na.omit(AgeYears)),
    proportion = first(na.omit(proportion)),
    .groups = "drop"
  ) %>%
  mutate(
    panel_lab = paste0(animalID, " (Age ", AgeYears, ", MOA=", sprintf("%.2f", proportion), ")")
  )

plot_df <- plot_df %>%
  left_join(label_df, by = c("animalID", "season"))

panel_order <- label_df %>%
  arrange(desc(AgeYears), animalID) %>%
  pull(panel_lab)

# 7) Create one shared daily grid
all_dates <- seq(min(plot_df$date, na.rm = TRUE),
                 max(plot_df$date, na.rm = TRUE),
                 by = "day")

date_lookup <- tibble(
  date = all_dates,
  day_index = seq_along(all_dates)
)

plot_df_full <- plot_df %>%
  select(panel_lab, date, pup_status, withpup, area) %>%
  group_by(panel_lab) %>%
  complete(date = all_dates) %>%
  ungroup() %>%
  left_join(date_lookup, by = "date") %>%
  mutate(
    panel_lab = factor(panel_lab, levels = panel_order),
    y = 1
  )

# choose actual dates to label, but keep every day as a tile
label_every <- 3   # try 2, 3, or 5
x_breaks <- date_lookup$day_index[seq(1, nrow(date_lookup), by = label_every)]
x_labels <- format(date_lookup$date[seq(1, nrow(date_lookup), by = label_every)], "%b %d")

conceptual_MOA_plot <- ggplot(
  plot_df_full,
  aes(x = day_index, y = y, fill = pup_status)
) +
  geom_tile(
    data = subset(plot_df_full, !is.na(pup_status)),
    width = 0.98,
    height = 7,
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(
    data = subset(plot_df_full, withpup %in% c("0", "1")),
    aes(label = withpup),
    size = 10,
    color = "black",
    fontface = "bold",
    vjust = -0.3
  ) +
  geom_text(
    data = subset(plot_df_full, !is.na(pup_status) & !is.na(area)),
    aes(label = area),
    size = 5,
    color = "black",
    vjust = 1.8
  ) +
  facet_grid(
    panel_lab ~ .,
    switch = "y",
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_manual(
    name = "Pup status",
    values = c(
      "approximate birth" = "#f2b6c6",
      "0" = "#8fd0eb",
      "1" = "#86d37c"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  labs(x = "Date", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 20) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0,
      hjust = 0,
      size = 35,
      face = "bold",
      margin = margin(r = 12)
    ),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 30
    ),
    axis.title.x = element_text(size = 40, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.spacing.y = unit(1.5, "lines"),
    legend.position = "right",
    legend.title = element_text(size = 35, face = "bold"),
    legend.text = element_text(size = 35),
  ); conceptual_MOA_plot

ggsave(
  "./TablesFigures/conceptual_MOA_plot.pdf",
  conceptual_MOA_plot,
  width = 49,
  height = 20
)

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
    name = "Mean location density
(seals per 10m radius)"
  ) +
  coord_sf(
    xlim = c(bbox_expanded["xmin"], bbox_expanded["xmax"]),
    ylim = c(bbox_expanded["ymin"], bbox_expanded["ymax"]),
    expand = FALSE
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ); seal_density_map

ggsave("./TablesFigures/seal_density_map.png", plot = seal_density_map, width = 30, height = 30, dpi = 300, bg = "white")

