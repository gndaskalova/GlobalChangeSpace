# Gergana Daskalova
# April 2021
# gndaskalova@gmail.com

# This script aims to extract temperature data over time for marine locations
# part of the BioTIME and Living Planet Databases

# PREDICTS locations not included because those data are not temporal and are only terrestrial

# Load packages
library(raster)
library(ncdf4)
library(tidyverse)

# The file here is too big for GitHub, you can download it from:
# https://psl.noaa.gov/repository/entry/show/PSD+Climate+Data+Repository/Public/PSD+Datasets/NOAA+OI+SST/Weekly+and+Monthly/sst.mnmean.nc?entryid=cac1c2a6-a864-4409-bb77-1fdead8eeb6e&output=default.html
# Please look for the file called sst.mnmean.nc

# Load the NOAA into R 
tmp <- brick("data/input/sst.mnmean.nc", varname="sst", package = "raster") # Mean monthly temperature

# Import sample site information
load("data/input/biotime.RData")

samples <- bt %>%
  filter(realm == "Marine") %>%
  dplyr::select(timeseries_id, long, lat) %>%
  distinct()

samples_simple <- samples %>% dplyr::select(long, lat)
colnames(samples_simple) <- c("lon", "lat")

# Extract climate data from the RasterBrick as a data.frame
tmp.sites <- data.frame(raster::extract(tmp, samples_simple, ncol = 2)) # Mean monthly temperature

# Add sample site names to the data.frame
tmp.sites$timeseries_id <- samples$timeseries_id

# Save the extracted climate data to a .RData file
save(tmp.sites, file = "data/output/NOAA_BioTIME2022.RData")

# Turn into long format
tmp.sites.long <- tmp.sites %>% gather(year, temperature, c(1:481))

# Calculate yearly averages
tmp.sites.long2 <- tmp.sites.long %>%
  separate(year, ".")

tmp.sites.long2$year <- parse_number(tmp.sites.long2$year)

tmp.sites.long2 <- tmp.sites.long2 %>% dplyr::group_by(rarefyID, year) %>%
  dplyr::summarise(mean_temp = mean(temperature))
sst_sites_long <- tmp.sites.long2

save(sst_sites_long, file = "data/output/NOAA_BioTIME_mean2022.RData")

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

samples_lpd_simple <- samples_lpd %>% dplyr::select(long, lat)
colnames(samples_lpd_simple) <- c("lon", "lat")
# Extract climate data from the RasterBrick as a data.frame
tmp.sites.lpd <- data.frame(raster::extract(tmp, samples_lpd_simple, ncol = 2)) # Mean monthly temperature

# Add sample site names to the data.frame
tmp.sites.lpd$timeseries_id <- samples_lpd$timeseries_id

# Save the extracted climate data to a .RData file
save(tmp.sites.lpd, file = "data/output/NOAA_LPD2022.RData")

# Turn into long format
tmp.sites.long.lpd <- tmp.sites.lpd %>% gather(year, temperature, c(1:481))

# Calculate yearly averages
tmp.sites.long.lpd2 <- tmp.sites.long.lpd %>%
  separate(year, ".")

tmp.sites.long.lpd2 <- tmp.sites.long.lpd %>%
  separate(year, ".")

# Warning messages are fine, the code discarded the day and month data

colnames(tmp.sites.long.lpd2)[2] <- "year"

tmp.sites.long.lpd2$year <- parse_number(tmp.sites.long.lpd2$year)

tmp.sites.long.lpd2 <- tmp.sites.long.lpd2 %>% dplyr::group_by(timeseries_id, year) %>%
  dplyr::summarise(mean_temp = mean(temperature))
sst_sites_long_lpd <- tmp.sites.long.lpd2

save(sst_sites_long_lpd, file = "data/output/NOAA_LPD_mean2022.RData")

