# Collate population and biodiversity data with driver data
# Gergana Daskalova
# 1st March 2019
# gndaskalova@gmail.com

# The goal of this script is to extract the intensities of different global change drivers
# for the locations represented by the BioTIME, Living Planet and PREDICTS databases

# The driver data come from Bowler et al. 2020 People and Nature.

library(raster)
library(tidyverse)
library(gridExtra)
library(rgdal)
library(dggridR)
# might need to install dggridR from GitHub if you don't have it installed already
# library(devtools) # Use `install.packages('devtools')` if need be
# install_github('r-barnes/dggridR', vignette=TRUE)
library(sp)
library(ggalt)
library(ggthemes)

# Load data ----
# Population data from Living Planet Database
mus <- read.csv("data/input/mus_1stMarch.csv")

# Assemblage data from BioTIME)
load("data/input/rarefied_medians2018_meta.RData")
load("data/input/rarefied_mediansOct2017.Rdata")
load("data/input/rarefied_medians2018.Rdata")

# Space for time community data from PREDICTS
predicts <- read.csv("data/input/predicts.csv")

# Drivers ----
# Note the .gri file has to be in the same folder as the grd one
cumulative <- stack("data/input/Cumulative_4_GD_imputeZero.grd")
atc_terr <- raster("data/input/clusterRaster_Trank06.tif")
atc_mar <- raster("data/input/clusterRaster_Mrank06.tif")

# Data filtering ----
# Include only time-series with 5 or more survey points
# Include only marine and terrestrial studies
# Include only taxa with enough representation, e.g. >40 time-series

bt <- rarefied_medians %>% 
  filter(REALM != "Freshwater" & duration > 4) %>%
  dplyr::select(REALM, BIOME_MAP, TAXA, duration, startYear, endYear, rarefyID_x,
         rarefyID_y, Jtu_base, YEAR, rarefyID, STUDY_ID) %>%
  mutate(year.test = YEAR == endYear) %>%
  filter(year.test == TRUE) %>%
  dplyr::select(-year.test, -YEAR) %>%
  filter(!TAXA %in% c("Marine plants", "Reptiles"))

bt$metric <- "Turnover"
bt$type <- "Biodiversity"
bt$species <- "Community" # Because these data are not species-based

bt <- bt %>% dplyr::select(type, rarefyID, STUDY_ID,
                           REALM, BIOME_MAP, TAXA, species, duration, 
                           startYear, endYear, rarefyID_x,
                           rarefyID_y, metric, Jtu_base)

colnames(bt) <- c("type", "timeseries_id", "study_id",
                  "realm", "biome", "taxa", "species", "duration", "start_year",
                  "end_year", "long", "lat", "metric", "value")

bt2 <- rarefied_medians2018 %>% 
  filter(REALM != "Freshwater" & duration > 4) %>%
  # Excluding freshwater data because there isn't enough driver data available for them
  dplyr::select(REALM, TAXA, BIOME_MAP, duration, startYear, endYear,
                Jtu_base, YEAR, rarefyID, STUDY_ID) %>%
  mutate(year.test = YEAR == endYear) %>%
  filter(year.test == TRUE) %>%
  dplyr::select(-year.test, -YEAR) %>%
  filter(!TAXA %in% c("Marine plants", "Reptiles"))

bt2$metric <- "Turnover"
bt2$type <- "Biodiversity"
bt2$species <- "Community"

# Standardise colunmn names and create one combined file
coords <- new_rarefied_medians_sept2_meta %>% 
  dplyr::select(STUDY_ID, CENT_LAT, CENT_LONG) %>% distinct()

bt2 <- left_join(bt2, coords, by = "STUDY_ID")

bt2 <- bt2 %>% dplyr::select(type, rarefyID, STUDY_ID,
                           REALM, BIOME_MAP, TAXA, species, duration, 
                           startYear, endYear, CENT_LONG,
                           CENT_LAT, metric, Jtu_base) %>%
  drop_na(CENT_LAT) %>% distinct()

colnames(bt2) <- c("type", "timeseries_id", "study_id",
                  "realm", "biome", "taxa", "species", "duration", "start_year",
                  "end_year", "long", "lat", "metric", "value")

bt <- rbind(bt, bt2)

mus <- mus %>%
  filter(Class %in% c("Aves", "Amphibia", "Mammalia", "Reptilia",
                             "Actinopterygii", "Elasmobranchii"))
mus$metric <- "Population change"
mus$type <- "Population"
mus$study_id <- mus$id
mus$species <- paste(mus$Genus, mus$Species, sep = " ")

mus <- mus %>%
  dplyr::select(type, id, study_id, system, biome, Class, species, duration.x, startYear,
                endYear, Decimal.Longitude, Decimal.Latitude,
                metric, mu)

colnames(mus) <- c("type", "timeseries_id", "study_id",
                  "realm", "biome", "taxa", "species", "duration", "start_year",
                  "end_year", "long", "lat", "metric", "value")

popbio <- rbind(mus, bt)

# Extract driver data -----
coords_sp <- SpatialPoints(cbind(popbio$long, popbio$lat), proj4string = CRS("+proj=longlat"))

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
coords3 <- coords %>% gather(type, lon, select = c(1:2))
coords3 <- coords3 %>% gather(direction, lat, select = c(1:2))
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
save(popbio, file = "data/output/popbio.RData")

# Adding PREDICTS locations
predicts <- predicts %>% drop_na(Latitude)
coords_sp <- SpatialPoints(cbind(predicts$Longitude, predicts$Latitude), proj4string = CRS("+proj=longlat"))

# Transforming coordinate system to that of the driver data
coords_sp <- spTransform(coords_sp, CRS("+proj=eck4 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
coords_df <- as.data.frame(coords_sp@coords)

save(coords_sp, file = "data/input/predicts_coords_simple.RData")
save(coords_df, file = "data/input/predicts_coords_df_simple.RData")

predicts_ids <- predicts %>% dplyr::select(SS, Longitude, Latitude)
save(predicts_ids, file = "data/input/predicts_ids.RData")

colnames(coords_df) <- c("long", "lat")

# Creating corners at the appropriate locations around the locations of the time-series
coords_df$atop <- coords_df$lat + 0.04413495
coords_df$bottom <- coords_df$lat - 0.04413495
coords_df$leftb <- coords_df$long - 0.04413495
coords_df$left <- coords_df$long + 0.04413495

# Creating spatial polygons
coords_df$study_id <- predicts$Study_number
coords <- coords_df %>% dplyr::select(left, leftb, atop, bottom, study_id)
coords3 <- coords %>% gather(type, lon, select = c(1:2))
coords3 <- coords3 %>% gather(direction, lat, select = c(1:2))
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

save(df_to_spp, file = "data/input/predicts_spp.RData")
save(coords4, file = "data/input/predicts_coords.RData")

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
save(predicts_drivers, file = "data/output/predicts_drivers.RData")

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
