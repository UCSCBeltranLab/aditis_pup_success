
  library(dplyr)
  library(sf)
  library(purrr)
  
  read_one_date <- function(d, y) {
    
    filepath_adults <- file.path("./RawData/Picterra_outputs_og", d, "adults", "adults.shp")
    filepath_pups   <- file.path("./RawData/Picterra_outputs_og", d, "pups",   "pups.shp")
    
    out <- list()
    
    if (file.exists(filepath_adults)) {
      out[[length(out) + 1]] <- st_read(filepath_adults, quiet = TRUE) %>%
        mutate(date = d, year = y, class = "adults")
    }
    
    if (file.exists(filepath_pups)) {
      out[[length(out) + 1]] <- st_read(filepath_pups, quiet = TRUE) %>%
        mutate(date = d, year = y, class = "pups")
    }
    
    if (length(out) == 0) return(NULL)
    
    x <- bind_rows(out)   # <- key change (sf-safe)
    
    # set CRS if missing (adjust if not lon/lat)
    if (is.na(st_crs(x))) st_crs(x) <- 4326
    
    x
  }
  
  picterra.list <- map2(dates, years, read_one_date)
  picterra.list <- compact(picterra.list)
  
  picterra.output <- bind_rows(picterra.list)  # <- key change (sf-safe)
  
  # confirm it stayed sf
  stopifnot(inherits(picterra.output, "sf"))
  

# Convert to projected CRS -------------------------------------------------
# NOTE: 26943 is NOT California Zone 3; see note below.
picterra.output <- st_transform(picterra.output, 26943)
beaches <- st_transform(beaches, 26943)

# Assign age/sex class to seals --------------------------------------------
uas.data <- picterra.output %>%
  filter(
    (class == "adults" & area_m2 < 3.8 & area_m2 > 1.3) |
      (class == "pups"  & area_m2 < 0.6 & area_m2 > 0.14)
  ) %>%
  mutate(class = as.factor(case_when(
    class == "adults" & area_m2 < 2.1  ~ "female",
    class == "adults" & area_m2 >= 2.1 ~ "male",
    class == "pups"                   ~ "pup"
  )))

# Assign seals to beach ----------------------------------------------------
centroids <- st_centroid(uas.data)
seals_by_beach <- st_intersection(beaches, centroids)

# Calculate female density -------------------------------------------------
density.df <- NULL

for (i in seq_along(dates)) {
  
  survey.subset <- seals_by_beach %>%
    filter(date == dates[i], class == "female")
  
  if (nrow(survey.subset) == 0) next
  
  seal.centroids <- st_centroid(survey.subset)
  seal.buffer <- st_buffer(seal.centroids, 10)
  
  int <- st_intersects(seal.buffer, seal.centroids)
  survey.subset$density <- lengths(int) - 1
  
  density.df <- if (is.null(density.df)) survey.subset else dplyr::bind_rows(density.df, survey.subset)
  

# Save data ---------------------------------------------------------------
seal.density <- density.df %>%
  select(date, class, lat, lon, Beach, density) %>%
  st_drop_geometry()

write.csv(seal.density, "./IntermediateData/seal.density.csv", row.names = FALSE)



