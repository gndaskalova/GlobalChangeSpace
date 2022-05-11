# Gergana Daskalova
# April 2021
# gndaskalova@gmail.com

# The goal of this script is to extract the name of the marine ecoregion
# in which the marine data points from the Living Planet database are located

# For terrestrial ecoregion extraction for the LPD, see the Geographic representation
# section of the script 07-plot-figures.R

library(rgdal)
library(raster)
library(tidyverse)

# Change working directory temporarily to be set to where the marine ecoregion
# shape files are
# Note that if you are changing the working directory on a Windows computer
# the formatting of the file paths differ, e.g. no "~"
# You can also click on Session/Working directory/Change directory and navigate
# to the folder where the ecoregion data are saved
setwd("data/input/MEOW-TNC")
ogrInfo(".", "meow_ecos")

regions <- readOGR(".", "meow_ecos")

#plot(regions, axes=TRUE, border="gray")

# Function to extract the ecoregions
getRegionalInfo  <- function(lat1, long1){

  #first, extract the co-ordinates (x,y - i.e., Longitude, Latitude)
  coords <- cbind(long1, lat1)
  
  FB.sp <- SpatialPointsDataFrame(coords,data.frame(value = c(4)))
  
  proj4string(FB.sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  dsdat <- over(regions, FB.sp, add=T, fn = mean) 
  
  ret <- data.frame(ECOREGION = regions$ECOREGION[which(dsdat$value==4)],
                    PROVICE = regions$PROVINCE[which(dsdat$value==4)],
                    REALM = regions$REALM[which(dsdat$value==4)])
  
  if(nrow(ret)==0) ret <- data.frame(ECOREGION = NA,
                                     PROVICE = NA,
                                     REALM = NA)
  return(ret)
  
}

# Change back the working directory to the project
setwd("~/GlobalChangeSpace")
# Note that if you have saved the repository for this project somewhere else
# the file path would need to be updated accordingly
# You can also click on Session/Working directory/Change directory and navigate
# to the folder with the repo

# For just one point
# getRegionalInfo(50.00806, -127.4342)

# For all the points in the marine LPD

# Load data
# Import sample site information
# For the Living Planet Database
# Import sample site information
lpd <- read.csv("data/input/LPR2020data_public.csv")
lpd$type <- "Population"

# Turn data into long form
lpd <- lpd %>% gather(year, pop, 30:98)
lpd$year <- parse_number(as.character(lpd$year))
lpd$pop <- as.factor(lpd$pop)
levels(lpd$pop)[levels(lpd$pop) == "NULL"] <- NA
lpd <- lpd %>% drop_na(pop)
lpd$pop <- parse_number(as.character(lpd$pop))

# Calculate duration per time series
lpd <- lpd %>% group_by(ID) %>% 
  mutate(duration = max(year) - min(year),
         startYear = min(year),
         endYear = max(year)) %>%
  filter(System != "Freshwater")

lpd <- lpd %>% gather(realm_type, biome, c(22, 25))

lpd <- lpd %>%
  dplyr::select(type, ID, System, biome, Class, duration, startYear,
                endYear, Longitude, Latitude)

colnames(lpd) <- c("type", "timeseries_id",
                   "realm", "biome", "taxa", "duration", "start_year",
                   "end_year", "long", "lat")

samples_lpd <- lpd %>%
  filter(realm == "Marine") %>%
  dplyr::select(timeseries_id, long, lat) %>%
  distinct() %>% ungroup()

lpd.coords.marine.simple <- samples_lpd %>%
  dplyr::select(timeseries_id, lat, long) %>% distinct()

colnames(lpd.coords.marine.simple)[c(2,3)] <- c("lat1", "long1")

marine_ecoregions_lpd <- lpd.coords.marine.simple %>% 
  group_by(timeseries_id) %>% 
  do(getRegionalInfo(.$lat1, .$long1))
length(unique(marine_ecoregions_lpd$ECOREGION))

# 188 ecoregions
# Just the number is needed to test % representation
# There are data from 188 ecoregions for the marine Living Planet Database

# Extract the number of ecoregions for the BioTIME marine time series
load("data/output/popbio2022.RData") # Living Planet and BioTIME databases
popbio <- popbio %>% filter(biome != "NULL") # removing blank duplicates

bt_mar <- popbio %>%
  filter(type == "Biodiversity" & realm == "Marine") %>%
  dplyr::select(timeseries_id, lat, long) %>% distinct() %>% ungroup()

colnames(bt_mar)[c(2,3)] <- c("lat1", "long1")

marine_ecoregions_bt <- bt_mar %>% 
  group_by(timeseries_id) %>% 
  do(getRegionalInfo(.$lat1, .$long1))
length(unique(marine_ecoregions_bt$ECOREGION))

# 102
