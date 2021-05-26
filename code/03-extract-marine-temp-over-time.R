# Gergana Daskalova
# April 2021
# gndaskalova@gmail.com

# This script aims to extract temperature data over time for marine locations
# part of the BioTIME and Living Planet Databases

# PREDICTS locations not included because those data are not temporal and are only terrestrial


# Load packages
library(raster)
library(ncdf4)

# The file here is too big for GitHub, you can download it from:
# https://psl.noaa.gov/repository/entry/show/PSD+Climate+Data+Repository/Public/PSD+Datasets/NOAA+OI+SST/Weekly+and+Monthly/sst.mnmean.nc?entryid=cac1c2a6-a864-4409-bb77-1fdead8eeb6e&output=default.html
# Please look for the file called sst.mnmean.nc

# Changing the working directory temporarily because for me the file is in Downloads
setwd("~/Downloads")

# Load the NOAA into R 
tmp <- brick("sst.mnmean.nc", varname="sst", package = "raster") # Mean monthly temperature

# Change back directory to the project
setwd("~/GlobalChangeSpace")

# Import sample site information
load("data/input/cell_coords_newsept.Rdata")

samples <- rarefyID_cell_centre %>%
  left_join(realm_meta) %>%
  filter(realm == "Marine") %>%
  dplyr::select(rarefyID, rarefyID_x, rarefyID_y) %>%
  distinct()

samples_simple <- samples %>% dplyr::select(rarefyID_x, rarefyID_y)
colnames(samples_simple) <- c("lon", "lat")

# Extract climate data from the RasterBrick as a data.frame
tmp.sites <- data.frame(raster::extract(tmp, samples_simple, ncol = 2)) # Mean monthly temperature

# Add sample site names to the data.frame
tmp.sites$rarefyID <- samples$rarefyID

# Save the extracted climate data to a .RData file
save(tmp.sites, file = "data/output/NOAA_BioTIME.RData")

# Turn into long format
tmp.sites.long <- tmp.sites %>% gather(year, temperature, c(1:2007))

# Calculate yearly averages
tmp.sites.long2 <- tmp.sites.long %>%
  separate(year, ".")

tmp.sites.long2$year <- parse_number(tmp.sites.long2$year)

tmp.sites.long2 <- tmp.sites.long2 %>% dplyr::group_by(rarefyID, year) %>%
  dplyr::summarise(mean_temp = mean(temperature))
sst_sites_long <- tmp.sites.long2

save(sst_sites_long, file = "data/output/NOAA_BioTIME_mean.RData")

# For the Living Planet Database
# Import sample site information
load("data/output/popbio.RData") # Living Planet and BioTIME databases

samples_lpd <- popbio %>%
  filter(type == "Population" & realm == "Marine") %>%
  dplyr::select(timeseries_id, long, lat) %>%
  distinct()

samples_lpd_simple <- samples_lpd %>% dplyr::select(long, lat)
colnames(samples_lpd_simple) <- c("lon", "lat")

# Extract climate data from the RasterBrick as a data.frame
tmp.sites.lpd <- data.frame(raster::extract(tmp, samples_lpd_simple, ncol = 2)) # Mean monthly temperature

# Add sample site names to the data.frame
tmp.sites.lpd$timeseries_id <- samples_lpd$timeseries_id

# Save the extracted climate data to a .RData file
save(tmp.sites.lpd, file = "data/output/NOAA_LPD.RData")

# Turn into long format
tmp.sites.long.lpd <- tmp.sites.lpd %>% gather(year, temperature, c(1:2007))

# Calculate yearly averages
tmp.sites.long.lpd2 <- tmp.sites.long.lpd %>%
  separate(year, ".")
colnames(tmp.sites.long.lpd2)[2] <- "year"

tmp.sites.long.lpd2$year <- parse_number(tmp.sites.long.lpd2$year)

tmp.sites.long.lpd2 <- tmp.sites.long.lpd2 %>% dplyr::group_by(timeseries_id, year) %>%
  dplyr::summarise(mean_temp = mean(temperature))
sst_sites_long_lpd <- tmp.sites.long.lpd2

save(sst_sites_long_lpd, file = "data/output/NOAA_LPD_mean.RData")

