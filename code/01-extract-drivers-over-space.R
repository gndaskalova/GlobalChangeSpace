# Collate population and biodiversity data with driver data
# Gergana Daskalova
# 12th Jan 2022
# gndaskalova@gmail.com

# The goal of this script is to extract the intensities of different global change drivers
# for the locations represented by the BioTIME, Living Planet and PREDICTS databases

# The driver data come from Bowler et al. 2020 People and Nature.

library(raster)
library(tidyverse)
library(gridExtra)
library(rgdal)
library(sp)
library(ggalt)
library(ggthemes)
library(dggridR)
library(rasterVis)
# might need to install dggridR from GitHub if you don't have it installed already
# library(devtools) # Use `install.packages('devtools')` if need be
# the dggridR package could take a while to install depending on your internet
# install_github('r-barnes/dggridR', vignette=TRUE)

# Load data ----
# Population data from Living Planet Database
mus <- read.csv("data/input/LPR2020data_public.csv")

# Assemblage data from BioTIME
load("data/input/bt_grid_coord.RData")

# Space for time community data from PREDICTS
# Note this file is not on GitHub because of its size
# The file can be downloaded here 
# https://data.nhm.ac.uk/dataset/the-2016-release-of-the-predicts-database
# Choose Site level summaries in zipped CSV format
# The file is originally called sites.csv, I renamed it to predicts.csv
# The file path below needs to be updated if the file is put in another location
predicts <- read.csv("data/input/predicts_sites.csv")

# Drivers
# Note the .gri file has to be in the same folder as the grd one
cumulative <- stack("data/input/Cumulative_4_GD_imputeZero.grd")
atc_terr <- raster("data/input/clusterRaster_Trank06.tif")
atc_mar <- raster("data/input/clusterRaster_Mrank06.tif")

# Data filtering ----
# Include only time-series with 5 or more survey points
# Include only marine and terrestrial studies
# Include only taxa with enough representation, e.g. >40 time-series
# Excluding freshwater data because there isn't enough driver data available for them

# The two objects bt and bt2 are to combine the BioTIME database with a few additional
# new studies that were added later

bt <- bt_grid_coord %>% 
  filter(REALM != "Freshwater" & duration > 4) %>%
  dplyr::select(REALM, BIOME_MAP, TAXA, duration, startYear, endYear, rarefyID_x,
         rarefyID_y, rarefyID) %>%
  filter(!TAXA %in% c("Fungi", "Reptiles"))

bt$type <- "Biodiversity"
bt <- bt %>% dplyr::select(type, rarefyID, REALM, BIOME_MAP,
                           TAXA, duration, startYear, endYear,
                           rarefyID_x, rarefyID_y)

colnames(bt) <- c("type", "timeseries_id",
                  "realm", "biome", "taxa", "duration", "start_year",
                  "end_year", "long", "lat")

# commented out to keep all time series regardless of Class
# mus <- mus %>%
#  filter(Class %in% c("Aves", "Amphibia", "Mammalia", "Reptilia",
#                             "Actinopterygii", "Elasmobranchii"))

mus$type <- "Population"

# Turn data into long form
mus <- mus %>% gather(year, pop, 30:98)
mus$year <- parse_number(as.character(mus$year))
mus$pop <- as.factor(mus$pop)
levels(mus$pop)[levels(mus$pop) == "NULL"] <- NA
mus <- mus %>% drop_na(pop)
mus$pop <- parse_number(as.character(mus$pop))

# Calculate duration per time series
mus <- mus %>% group_by(ID) %>% 
  mutate(duration = max(year) - min(year),
         startYear = min(year),
         endYear = max(year)) %>%
  filter(System != "Freshwater")

mus <- mus %>% gather(realm_type, biome, c(22, 25))

mus <- mus %>%
  dplyr::select(type, ID, System, biome, Class, duration, startYear,
                endYear, Longitude, Latitude)

colnames(mus) <- c("type", "timeseries_id",
                  "realm", "biome", "taxa", "duration", "start_year",
                  "end_year", "long", "lat")

# Check columns before combining LPD and BioTIME
colnames(bt)
bt$timeseries_id <- as.character(bt$timeseries_id)
mus$timeseries_id <- as.character(mus$timeseries_id)

popbio <- rbind(mus, bt)
#save(popbio, file ="data/input/popbio2020.RData")

# Extract driver data -----
coords_sp <- SpatialPoints(cbind(popbio$long, popbio$lat), proj4string = CRS("+proj=longlat"))

# The driver data are extracted over equal cells of around 96km2
# Create the spatial polygons for the extraction

# Transforming coordinate system to that of the driver data
coords_sp <- spTransform(coords_sp, CRS("+proj=eck4 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
coords_df <- as.data.frame(coords_sp@coords)
colnames(coords_df) <- c("long", "lat")

# Creating corners at the appropriate locations around the locations of the time-series
coords_df$atop <- coords_df$lat + 0.04413495
coords_df$bottom <- coords_df$lat - 0.04413495
coords_df$leftb <- coords_df$long - 0.04413495
coords_df$left <- coords_df$long + 0.04413495

# Creating spatial polygons
coords_df$timeseries_id <- popbio$timeseries_id
coords <- coords_df %>% dplyr::select(left, leftb, atop, bottom, timeseries_id)
coords3 <- coords %>% gather(type, lon, 1:2)
coords3 <- coords3 %>% gather(direction, lat, 1:2)
coords3 <- coords3 %>% arrange(timeseries_id)

coords4 <- coords3 %>% group_by(timeseries_id) %>% filter(row_number() == 1)
coords4$order <- "2"
coords3$order <- "1"
coords5 <- full_join(coords3, coords4)
coords5$sort <- paste0(coords5$type, coords5$direction)
coords5 <- coords5 %>% arrange(timeseries_id, order, sort)

df_to_spp <- coords5 %>%
  group_by(timeseries_id) %>%
  do(poly = dplyr::select(., lon, lat) %>% Polygon()) %>%
  rowwise() %>%
  do(polys = Polygons(list(.$poly),.$timeseries_id)) %>%
  {SpatialPolygons(.$polys)}

# plot(df_to_spp) # Check distribution, takes a while to plot!

# Climate change
cc <- subset(cumulative, 1)
# Human use
hu <- subset(cumulative, 2)
# Human population
hp <- subset(cumulative, 3)
# Pollution
po <- subset(cumulative, 4)
# Invasions
inv <- subset(cumulative, 5)
# Cumulative
all <- subset(cumulative, 6)

# Plot drivers
levelplot(all)

# Note that this stage can take a while
df13 <- as.data.frame(raster::extract(cc, df_to_spp, fun = mean))
colnames(df13) <- "climate_change"
df14 <- as.data.frame(raster::extract(hu, df_to_spp, fun = mean))
colnames(df14) <- "human_use"
df15 <- as.data.frame(raster::extract(hp, df_to_spp, fun = mean))
colnames(df15) <- "human_population"
df16 <- as.data.frame(raster::extract(po, df_to_spp, fun = mean))
colnames(df16) <- "pollution"
df17 <- as.data.frame(raster::extract(inv, df_to_spp, fun = mean))
colnames(df17) <- "invasions"
df18 <- as.data.frame(raster::extract(all, df_to_spp, fun = mean))
colnames(df18) <- "cumulative"

popbio <- distinct(popbio)

ids <- coords4 %>% dplyr::select(timeseries_id)
drivers <- cbind(df13, df14, df15, df16, df17, df18)
drivers$timeseries_id <- ids$timeseries_id

popbio <- left_join(popbio, drivers, by = "timeseries_id")

# ATC extraction (combination of drivers, not used in the representation manuscript)
popbio_terr <- popbio %>% filter(realm == "Terrestrial")

coords_sp_terr <- SpatialPoints(cbind(popbio_terr$long, popbio_terr$lat), proj4string = CRS("+proj=longlat"))

# Transforming coordinate system to that of the driver data
coords_sp_terr <- spTransform(coords_sp_terr, CRS("+proj=eck4 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
coords_df_terr <- as.data.frame(coords_sp_terr@coords)
colnames(coords_df_terr) <- c("long", "lat")

popbio_terr_atc <- as.data.frame(raster::extract(atc_terr, coords_df_terr))
colnames(popbio_terr_atc) <- "ATC"
popbio_terr$ATC <- popbio_terr_atc$ATC

popbio_mar <- popbio %>% filter(realm == "Marine")
coords_sp_mar <- SpatialPoints(cbind(popbio_mar$long, popbio_mar$lat), proj4string = CRS("+proj=longlat"))
coords_sp_mar <- spTransform(coords_sp_mar, CRS("+proj=eck4 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
coords_df_mar <- as.data.frame(coords_sp_mar@coords)
colnames(coords_df_mar) <- c("long", "lat")

popbio_mar_atc <- as.data.frame(raster::extract(atc_mar, coords_df_mar))
colnames(popbio_mar_atc) <- "ATC"
popbio_mar$ATC <- popbio_mar_atc$ATC
popbio_mar_atc <- as.data.frame(raster::extract(atc_mar, coords_df_mar))

# Unite the marine and terrestrial again
popbio <- rbind(popbio_terr, popbio_mar)
save(popbio, file = "data/output/popbio2022.RData")

# Adding PREDICTS locations
predicts <- predicts %>% drop_na(Latitude)
coords_sp <- SpatialPoints(cbind(predicts$Longitude, predicts$Latitude), proj4string = CRS("+proj=longlat"))

# Transforming coordinate system to that of the driver data
coords_sp <- spTransform(coords_sp, CRS("+proj=eck4 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
coords_df <- as.data.frame(coords_sp@coords)

save(coords_sp, file = "data/input/predicts_coords_simple2022.RData")
save(coords_df, file = "data/input/predicts_coords_df_simple2022.RData")

predicts_ids <- predicts %>% dplyr::select(SS, Longitude, Latitude)
save(predicts_ids, file = "data/input/predicts_ids2022.RData")

colnames(coords_df) <- c("long", "lat")

# Creating corners at the appropriate locations around the locations of the time-series
coords_df$atop <- coords_df$lat + 0.04413495
coords_df$bottom <- coords_df$lat - 0.04413495
coords_df$leftb <- coords_df$long - 0.04413495
coords_df$left <- coords_df$long + 0.04413495

# Creating spatial polygons
coords_df$study_id <- predicts$Study_number
coords <- coords_df %>% dplyr::select(left, leftb, atop, bottom, study_id)
coords3 <- coords %>% gather(type, lon, 1:2)
coords3 <- coords3 %>% gather(direction, lat, 1:2)
coords3 <- coords3 %>% arrange(study_id)

coords4 <- coords3 %>% group_by(study_id) %>% filter(row_number() == 1)
coords4$order <- "2"
coords3$order <- "1"
coords5 <- full_join(coords3, coords4)
coords5$sort <- paste0(coords5$type, coords5$direction)
coords5 <- coords5 %>% arrange(study_id, order, sort)

df_to_spp <- coords5 %>%
  group_by(study_id) %>%
  do(poly = dplyr::select(., lon, lat) %>% Polygon()) %>%
  rowwise() %>%
  do(polys = Polygons(list(.$poly),.$study_id)) %>%
  {SpatialPolygons(.$polys)}

save(df_to_spp, file = "data/input/predicts_spp2022.RData")
save(coords4, file = "data/input/predicts_coords2022.RData")

# plot(df_to_spp) # Check distribution, takes a while to plot!

# Climate change
cc <- subset(cumulative, 1)
# Human use
hu <- subset(cumulative, 2)
# Human population
hp <- subset(cumulative, 3)
# Pollution
po <- subset(cumulative, 4)
# Invasions
inv <- subset(cumulative, 5)
# Cumulative
all <- subset(cumulative, 6)

df13 <- as.data.frame(raster::extract(cc, df_to_spp, fun = mean))
colnames(df13) <- "climate_change"
df14 <- as.data.frame(raster::extract(hu, df_to_spp, fun = mean))
colnames(df14) <- "human_use"
df15 <- as.data.frame(raster::extract(hp, df_to_spp, fun = mean))
colnames(df15) <- "human_population"
df16 <- as.data.frame(raster::extract(po, df_to_spp, fun = mean))
colnames(df16) <- "pollution"
df17 <- as.data.frame(raster::extract(inv, df_to_spp, fun = mean))
colnames(df17) <- "invasions"
df18 <- as.data.frame(raster::extract(all, df_to_spp, fun = mean))
colnames(df18) <- "cumulative"

ids <- coords4 %>% dplyr::select(study_id)
drivers <- cbind(df13, df14, df15, df16, df17, df18)
drivers$study_id <- ids$study_id

predicts_drivers <- drivers
save(predicts_drivers, file = "data/output/predicts_drivers2022.RData")

# Random sampling driver extraction ----

N <- 1000000  # How many cells to sample

dggs <- dgconstruct(res = 12, metric = FALSE, resround = 'nearest')

maxcell <- dgmaxcell(dggs)                     #Get maximum cell id
cells <- sample(1:maxcell, N, replace = FALSE) #Choose random cells
grid <- dgcellstogrid(dggs, cells, frame = TRUE, wrapcells = TRUE) #Get grid

#Get polygons for each country of the world
countries <- map_data("world")

#Plot everything on a flat map
(p <- ggplot() + 
    geom_polygon(data = countries, aes(x = long, y = lat, group = group), 
                 fill = "grey60", color = "grey60") +
    coord_proj("+proj=eck4") +
    theme_map() +
    geom_point(data = grid, aes(x = long, y = lat, group = group), 
               colour = "deepskyblue4", alpha = 0.05))

ggsave(p, filename = "figures/map_ideal_sampling.pdf", dpi = 300, device = "pdf",
       height = 5, width = 5)

# The world is covered!

world_coords <- grid %>% dplyr::select(cell, long, lat) %>% distinct() %>%
  group_by(cell) %>% slice(1)

colnames(world_coords) <- c("cell", "long", "lat")

coords_again <- SpatialPoints(world_coords[, c("long", "lat")])
world_spatial_df <- SpatialPointsDataFrame(coords_again, world_coords)

random_coords <- SpatialPoints(cbind(world_spatial_df$long, world_spatial_df$lat),
                               proj4string = CRS("+proj=longlat"))

random_coords <- spTransform(random_coords, CRS("+proj=eck4 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
random <- as.data.frame(random_coords@coords)
colnames(random) <- c("long", "lat")

df19 <- as.data.frame(raster::extract(cc, random))
colnames(df19) <- "climate_change"
df20 <- as.data.frame(raster::extract(hu, random))
colnames(df20) <- "human_use"
df21 <- as.data.frame(raster::extract(hp, random))
colnames(df21) <- "human_population"
df22 <- as.data.frame(raster::extract(po, random))
colnames(df22) <- "pollution"
df23 <- as.data.frame(raster::extract(inv, random))
colnames(df23) <- "invasions"
df24 <- as.data.frame(raster::extract(all, random))
colnames(df24) <- "cumulative"

random_drivers <- cbind(df19, df20, df21, df22, df23, df24)
random_drivers$lat <- world_coords$lat
random_drivers$long <- world_coords$long

random_drivers <- na.omit(random_drivers)

save(random_drivers, file = "data/output/random_samplingNov2020.RData")

# Plot random drivers
random_drivers_long <- random_drivers %>% gather(driver, value, 1:6)

(driver_densities <- ggplot(random_drivers_long) +
  geom_density(aes(x = lat, y = value)) +
  facet_grid(~driver))
