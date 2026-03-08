############ facet the two age figures together ###########

plot_age <- plot_grid(
  plot_age_1996_2025,
  plot_age_2016_2023,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1))

plot_age

ggsave("./TablesFigures/age_plot.png",
       plot = plot_age,
       width = 16,
       height = 20,
       dpi = 1000)

############ conceptual MOA figure #########

# 1) make sure dates are as Date and animalIDs are characters
metadata <- metadata %>%
  mutate(date = as.Date(date),
         BirthDate = as.Date(BirthDate),
         animalID = as.character(animalID))

# 2) Choose 3 animals from 2025
example_ids <- c("36681", "46960", "51527")

meta_subset <- metadata %>%
  filter(season == 2025, animalID %in% example_ids)

# 3) Resight tiles
resight_tiles <- meta_subset %>%
  filter(!is.na(date), !is.na(withpup)) %>%
  distinct(animalID, season, date, .keep_all = TRUE) %>%
  mutate(withpup = as.character(withpup),
         pup_status = case_when(withpup == "0" ~ "0",
                                withpup == "1" ~ "1",
                                TRUE ~ NA_character_)) %>%
  select(animalID, season, date, AgeYears, proportion, pup_status, withpup)

# 4) Approximate birth tiles
birth_tiles <- meta_subset %>%
  distinct(animalID, season, BirthDate, AgeYears, proportion, .keep_all = TRUE) %>%
  filter(!is.na(BirthDate)) %>%
  transmute(animalID,
            season,
            date = BirthDate,
            AgeYears,
            proportion,
            pup_status = "approximate birth",
            withpup = "1")

# 5) Combine tiles
# keep approximate birth if same day as a resight
plot_df <- bind_rows(resight_tiles, birth_tiles) %>%
  mutate(pup_status = factor(pup_status,
                             levels = c("approximate birth", "0", "1"))) %>%
  arrange(animalID, season, date, pup_status) %>%
  distinct(animalID, season, date, .keep_all = TRUE)

# 6) Panel labels
label_df <- plot_df %>%
  group_by(animalID, season) %>%
  summarise(AgeYears = first(na.omit(AgeYears)),
            proportion = first(na.omit(proportion)),
            .groups = "drop") %>%
  mutate(panel_lab = paste0(animalID,
                            " (Age ", AgeYears,
                            ", MOA=", sprintf("%.2f", proportion), ")"))

plot_df <- plot_df %>%
  left_join(label_df, by = c("animalID", "season"))

panel_order <- label_df %>%
  arrange(desc(AgeYears), animalID) %>%
  pull(panel_lab)

# 7) Create one shared daily grid for 2025
all_dates <- seq(min(plot_df$date, na.rm = TRUE),
                 max(plot_df$date, na.rm = TRUE),
                 by = "day")

date_lookup <- tibble(date = all_dates,
                      day_index = seq_along(all_dates))

plot_df_full <- plot_df %>%
  select(panel_lab, date, pup_status, withpup) %>%
  group_by(panel_lab) %>%
  complete(date = all_dates) %>%
  ungroup() %>%
  left_join(date_lookup, by = "date") %>%
  mutate(panel_lab = factor(panel_lab, levels = panel_order),
         y = 1)

# 8) Plot with actual dates on x-axis
conceptual_MOA_plot <- ggplot(plot_df_full, aes(x = day_index, y = y, fill = pup_status)) +
  geom_tile(data = subset(plot_df_full, !is.na(pup_status)),
            width = 0.9,
            height = 4,
            color = NA) +
  geom_text(data = subset(plot_df_full, withpup %in% c("0","1")),
            aes(label = withpup),
            size = 3,
            color = "black") +
  facet_grid(panel_lab ~ ., switch = "y") +
  scale_fill_manual(name = "pup status",
                    values = c("approximate birth" = "#f2b6c6",
                               "0" = "#8fd0eb",
                               "1" = "#86d37c"),
                    drop = FALSE) +
  scale_x_continuous(breaks = date_lookup$day_index,
                     labels = format(date_lookup$date, "%b %d"),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = NULL,
                     expand = c(0, 0)) +
  labs(x = "Date", y = NULL) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 0, size = 15, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.spacing.y = unit(0.8, "lines"),
        legend.position = "right")

conceptual_MOA_plot

# optional save
ggsave("./TablesFigures/conceptual_MOA_plot.png",
       plot = conceptual_MOA_plot,
       width = 16,
       height = 10,
       dpi = 1000)

