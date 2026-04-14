
# Load libraries and data -------------------------------------------------
library(dplyr)
library(lubridate)
library(sf)
library(purrr)

options(scipen = 999)
sf_use_s2(FALSE)

# read in beach shapefile
list.files("..", pattern = "Ano Nuevo.*\\.shp$", recursive = TRUE, full.names = TRUE)

beaches <- st_read("./RawData/ANM map/Ano Nuevo Map Final.shp", quiet = TRUE) %>%
  select(-id)

# read in seal detections
base_dir <- "./RawData/Picterra_outputs_updated"

shps <- list.files(base_dir, recursive = TRUE, pattern = "polygons\\.shp$", full.names = FALSE)

dates <- dirname(shps)
years <- substr(dates, 1, 4)

read_one <- function(rel_shp, d, y) {
  
  shp_path <- file.path(base_dir, rel_shp)
  
  x <- st_read(shp_path, quiet = TRUE) %>%
    mutate(date = d, year = y)
  
  if (is.na(st_crs(x))) st_crs(x) <- 4326
  x
}

picterra.list <- pmap(list(rel_shp = shps, d = dates, y = years), read_one) |> compact()

picterra.output <- do.call(rbind, picterra.list)

# Convert to NAD83 / California zone 3 ------------------------------------
st_crs(picterra.output)

picterra.output <- st_transform(picterra.output, "EPSG:26943")

st_crs(picterra.output)

beaches <- st_transform(beaches, "EPSG:26943")

# Assign age/sex class to seals -------------------------------------------
uas.data <- picterra.output %>%
  mutate(age_sex = tolower(as.character(age_sex))) %>%
  filter(age_sex %in% c("female", "male", "pup"))

# Assign seals to beach ---------------------------------------------------
centroids <- st_centroid(uas.data)

seals_by_beach <- st_intersection(beaches, centroids)

# Calculate female density ------------------------------------------------
density.df <- data.frame()

for (i in 1:length(dates)) {
  
  survey.subset <- seals_by_beach %>%
    filter(date == dates[i],
           age_sex == "female")
  
  if (nrow(survey.subset) == 0) next
  
  seal.centroids <- st_centroid(survey.subset)
  
  ## set to 10m, but you can modify that to whatever you would like!
  seal.buffer <- st_buffer(seal.centroids, 10)
  
  int <- st_intersects(seal.buffer, seal.centroids)
  
  survey.subset$density <- lengths(int) - 1
  
  density.df <- rbind(density.df, survey.subset)
}

# Save data ---------------------------------------------------------------

seal.density <- density.df %>%
  mutate(
    centroid = st_centroid(geometry),
    lon = st_coordinates(centroid)[,1],
    lat = st_coordinates(centroid)[,2]
  ) %>%
  select(date, age_sex, lat, lon, Beach, density) %>%
  st_drop_geometry()

write.csv(seal.density, "./RawData/seal.density.csv",
          row.names = FALSE)

# Plot an example ---------------------------------------------------------

###low density version 

# grab a LIGHT color from mako
light_mako <- viridis(190, option = "mako", direction = -1)[20]

# find a seal with 5 neighbors (or closest to 5 if none exist)
target_n <- 5
best_diff <- Inf
low_example <- NULL

for (i in 1:nrow(uas.data)) {
  
  focal <- uas.data[i, ]
  
  if (focal$age_sex[1] != "female") next
  
  same <- uas.data %>%
    filter(date == focal$date[1],
           age_sex == "female")
  
  hits <- st_intersects(
    st_buffer(st_centroid(focal), 10),
    st_centroid(same)
  )[[1]]
  
  n_neighbors <- length(hits) - 1
  
  if (n_neighbors == target_n) {
    low_example <- focal
    break
  }
  
  if (abs(n_neighbors - target_n) < best_diff) {
    best_diff <- abs(n_neighbors - target_n)
    low_example <- focal
  }
}

# plot
seal.buffer <- st_buffer(st_centroid(low_example), 10)

same <- uas.data %>%
  filter(date == low_example$date[1],
         age_sex == "female")

hits <- st_intersects(seal.buffer, st_centroid(same))[[1]]

int.poly <- same[hits[-1], ]   # drop focal seal itself

png("./TablesFigures/density_example_low.png",
    width = 2000,
    height = 2000,
    res = 300,
    bg = "transparent")   # <-- makes outside transparent

par(bg = NA, mar = c(0, 0, 0, 0))  # remove outer margins

plot(st_geometry(seal.buffer), border = "grey40")

# surrounding seals
plot(st_geometry(int.poly),
     add = TRUE,
     col = light_mako,
     border = "grey40")

# focal seal (draw last so it's on top)
plot(st_geometry(low_example),
     add = TRUE,
     col = light_mako,
     border = "black",
     lwd = 3)

dev.off()

###high density version

dark_mako <- viridis(50, option = "mako", direction = 1)[20]

# find a seal with MANY neighbors (target ~18, or closest)
target_n <- 18
best_diff <- Inf
high_example <- NULL

for (i in 1:nrow(uas.data)) {
  
  focal <- uas.data[i, ]
  
  if (focal$age_sex[1] != "female") next
  
  same <- uas.data %>%
    filter(date == focal$date[1],
           age_sex == "female")
  
  hits <- st_intersects(
    st_buffer(st_centroid(focal), 10),
    st_centroid(same)
  )[[1]]
  
  n_neighbors <- length(hits) - 1
  
  if (n_neighbors == target_n) {
    high_example <- focal
    break
  }
  
  if (abs(n_neighbors - target_n) < best_diff) {
    best_diff <- abs(n_neighbors - target_n)
    high_example <- focal
  }
}

# plot
seal.buffer <- st_buffer(st_centroid(high_example), 10)

same <- uas.data %>%
  filter(date == high_example$date[1],
         age_sex == "female")

hits <- st_intersects(seal.buffer, st_centroid(same))[[1]]
int.poly <- same[hits[-1], ]   # drop focal seal

png("./TablesFigures/density_example_high.png",
    width = 2000,
    height = 2000,
    res = 300,
    bg = "transparent")   # outside = transparent

par(bg = NA, mar = c(0, 0, 0, 0))  # remove margins

plot(st_geometry(seal.buffer), border = "grey40")

# surrounding seals
plot(st_geometry(int.poly),
     add = TRUE,
     col = dark_mako,
     border = "grey40")

# focal seal (on top)
plot(st_geometry(high_example),
     add = TRUE,
     col = dark_mako,
     border = "black",
     lwd = 5)

dev.off()

