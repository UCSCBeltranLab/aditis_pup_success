library(dplyr)
library(sf)
library(purrr)
library(stringr)
install.packages("rnaturalearth")
library(rnaturalearth)

options(scipen = 999)
sf_use_s2(FALSE)

# Read beaches 
beaches <- st_read("./RawData/ANM map/Ano Nuevo Map Final.shp", quiet = TRUE) %>%
  select(-id)

plot(st_geometry(beaches))

# Picterra outputs (already include "polygons.shp")
base_dir <- "./RawData/Picterra_outputs_updated"

# returns paths like "2016/polygons.shp", "20210129/polygons.shp", etc
shps  <- list.files(base_dir, recursive = TRUE, pattern = "polygons\\.shp$", full.names = FALSE)

# Extract "date" as the folder name right above polygons.shp
# e.g., "2016/polygons.shp" -> "2016"
dates <- dirname(shps)
years <- substr(dates, 1, 4)

read_one <- function(rel_shp, d, y) {
  
  shp_path <- file.path(base_dir, rel_shp)
  if (!file.exists(shp_path)) return(NULL)
  
  x <- st_read(shp_path, quiet = TRUE) %>%
    mutate(date = d, year = y)
  
  if (is.na(st_crs(x))) st_crs(x) <- 4326
  x
}

picterra.list <- pmap(list(rel_shp = shps, d = dates, y = years), read_one) |> compact()

stopifnot(length(picterra.list) > 0)

# rbind keeps sf class reliably
picterra.output <- do.call(rbind, picterra.list)
stopifnot(inherits(picterra.output, "sf"))

# Convert to projected CRS ------------------------------------------------
picterra.output <- st_transform(picterra.output, 26943)
beaches         <- st_transform(beaches, 26943)

# Use existing classifications in age_sex (as in your second script) ------
uas.data <- picterra.output %>%
  mutate(age_sex = tolower(as.character(age_sex))) %>%
  filter(age_sex %in% c("female", "male", "pup"))

# Assign seals to beach ---------------------------------------------------
centroids      <- st_centroid(uas.data)
seals_by_beach <- st_intersection(beaches, centroids)

# Calculate female density (FEMALES ONLY) --------------------------------
density.df <- list()

for (d in unique(dates)) {
  
  survey.subset <- seals_by_beach %>%
    filter(date == d, age_sex == "female")  # <-- females only
  
  if (nrow(survey.subset) == 0) next
  
  seal.centroids <- st_centroid(survey.subset)
  seal.buffer    <- st_buffer(seal.centroids, 10)
  
  int <- st_intersects(seal.buffer, seal.centroids)
  survey.subset$density <- lengths(int) - 1
  
  density.df[[d]] <- survey.subset
}

density.df <- bind_rows(density.df)

seal.density <- density.df %>%
  mutate(
    centroid = st_centroid(geometry),
    lon = st_coordinates(centroid)[,1],
    lat = st_coordinates(centroid)[,2]
  ) %>%
  select(date, age_sex, lat, lon, Beach, density) %>%
  st_drop_geometry()

write.csv(seal.density, "./IntermediateData/seal.density.csv", row.names = FALSE)

