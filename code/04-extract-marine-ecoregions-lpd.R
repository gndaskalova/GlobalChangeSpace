# Gergana Daskalova
# April 2021
# gndaskalova@gmail.com

# The goal of this script is to extract the name of the marine ecoregion
# in which the marine data points from the Living Planet database are located

# There is no need to do this for the BioTIME data because the data are already included
# in the metadata

# For terrestrial ecoregion extraction for the LPD, see the Geographic representation
# section of the script 07-plot-figures.R

library(rgdal)
library(raster) 

# Change working directory temporarily to be set to where the marine ecoregion
# shape files are
setwd("~/GlobalChangeSpace/data/input/MEOW-TNC")
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

# For just one point
# getRegionalInfo(50.00806, -127.4342)

# For all the points in the marine LPD

# Load data
# Import sample site information
load("data/output/popbio.RData") # Living Planet and BioTIME databases, taking just the 
# LPD marinecoordinates

mus.coords.marine <- popbio %>%
  filter(type == "Population" & realm == "Marine") %>%
  dplyr::select(timeseries_id, long, lat) %>%
  distinct()

mus.coords.marine.simple <- mus.coords.marine %>%
  dplyr::select(timeseries_id, lat, long) %>% distinct()

colnames(mus.coords.marine.simple)[c(2,3)] <- c("lat1", "long1")

marine_ecoregions_lpd <- mus.coords.marine.simple %>% group_by(timeseries_id) %>% do(getRegionalInfo(.$lat1, .$long1))
length(unique(marine_ecoregions_lpd$ECOREGION))

# 160 ecoregions
# Just the number is needed to test % representation
# There are data from 160 ecoregions for the marine Living Planet Database