# Collate population and biodiversity data with driver data
# Gergana Daskalova
# 1st March 2019

library(raster)
library(tidyverse)
library(gridExtra)
library(rgdal)

# Load data ----
# Population change
mus <- read.csv("data/input/mus_1stMarch.csv")

# Richness change and turnover
load("data/input/rarefied_medians2018_meta.RData")
load("data/input/rarefied_mediansOct2017.Rdata")
load("data/input/rarefied_medians2018.Rdata")

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
bt$species <- "Community"

bt <- bt %>% dplyr::select(type, rarefyID, STUDY_ID,
                           REALM, BIOME_MAP, TAXA, species, duration, 
                           startYear, endYear, rarefyID_x,
                           rarefyID_y, metric, Jtu_base)

colnames(bt) <- c("type", "timeseries_id", "study_id",
                  "realm", "biome", "taxa", "species", "duration", "start_year",
                  "end_year", "long", "lat", "metric", "value")

bt2 <- rarefied_medians2018 %>% 
  filter(REALM != "Freshwater" & duration > 4) %>%
  dplyr::select(REALM, TAXA, BIOME_MAP, duration, startYear, endYear,
                Jtu_base, YEAR, rarefyID, STUDY_ID) %>%
  mutate(year.test = YEAR == endYear) %>%
  filter(year.test == TRUE) %>%
  dplyr::select(-year.test, -YEAR) %>%
  filter(!TAXA %in% c("Marine plants", "Reptiles"))

bt2$metric <- "Turnover"
bt2$type <- "Biodiversity"
bt2$species <- "Community"

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

# ATC extraction
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
