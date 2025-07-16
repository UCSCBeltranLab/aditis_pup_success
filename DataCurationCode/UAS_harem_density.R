# Load libraries and data -------------------------------------------------

library(dplyr)
library(lubridate)
library(sf)

options(scipen = 999)

#read in beach shapefile

beaches <- st_read("../RawData/ANM map/Ano Nuevo Map Final.shp") %>%
  select(-id)

#read in seal detections
dates <- list.files("../RawData/Picterra_outputs/")

years <- substr(dates, 0, 4)

picterra.output <- data.frame()

for(i in 1:length(dates)) {
  filepathadults <- paste0("../RawData/Picterra_outputs/", dates[i], "/adults/adults.shp")
  
  filepathpups <- paste0("../RawData/Picterra_outputs/", dates[i], "/pups/pups.shp")
  
  picterra.output <- rbind(picterra.output, 
                           cbind(st_read(filepathadults), date = dates[i], year = years[i], class = "adults"),
                           cbind(st_read(filepathpups), date = dates[i], year = years[i], class = "pups"))
  
}



# Convert to NAD83 / California zone 3 -------------------------------------------------------------


st_crs(picterra.output)

picterra.output = st_transform(picterra.output,"EPSG:26943") #transform shape file 

st_crs(picterra.output)

beaches <- st_transform(beaches,"EPSG:26943") #transform shape file 


# Assign age/sex class to seals -----------------------------------------------


uas.data <- picterra.output %>%
  filter(class == "adults" & area_m2 < 3.8 &
           class == "adults" &area_m2 > 1.3 |
           class == "pups" & area_m2 <.6 &
           class == "pups" & area_m2 > .14) %>%
  mutate(class = as.factor(case_when(class == "adults" & area_m2 < 2.1 ~ "female", 
                                     class == "adults" & area_m2 >= 2.1 ~ "male",
                                     class == "pups"~ "pup")))

# Assign seals to beach ------------------------------------------------

centroids <- st_centroid(uas.data)

##turn off spherical geometry
sf_use_s2(F)

seals_by_beach <- st_intersection(beaches, centroids)




# Calculate female density ------------------------------------------------
density.df <- data.frame()

for (i in 1:length(dates)) {
  
  survey.subset <- seals_by_beach %>%
    filter(date == dates[i],
           class == "female")
  
  seal.centroids <- st_centroid(survey.subset)
  
  ##set to 10m, but you can modify that to whatever you would like!
  seal.buffer <- st_buffer(seal.centroids, 10)
  
  
  int <- st_intersects(seal.buffer, seal.centroids) 
  
  survey.subset$density <- lengths(int) - 1
  
  density.df <- rbind(density.df,survey.subset )
}



# Plot an example ---------------------------------------------------------

seal.buffer <- st_buffer(st_centroid(uas.data[1,]), 10)

plot(st_geometry(seal.buffer))

plot(st_geometry(uas.data[1,]), add = TRUE, col = "red")

int <- uas.data %>%
  filter(date == uas.data[[1, 6]]) %>%
  st_centroid() %>%
  st_intersection(seal.buffer, 10)


int.poly <- uas.data[row.names(uas.data) %in% row.names(int),]


plot(st_geometry(int.poly), add = TRUE)



# Save data ---------------------------------------------------------------


seal.density <- density.df %>%
 select(date, class, lat, lon, Beach, density) %>%
 st_drop_geometry()


write.csv(seal.density, "../IntermediateData/seal.density.csv",
          row.names = FALSE)

