# Gergana Daskalova
# April 2021
# gndaskalova@gmail.com

# This script aims to extract temperature data over time for terrestrial locations
# part of the BioTIME and Living Planet Databases

# PREDICTS locations not included because those data are not temporal

# Load packages
library(raster)
library(ncdf4)

# The file here is too big for GitHub, you can download it from:
#https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.01/cruts.1709081022.v4.01/tmp/

# Please look for the file called cru_ts4.01.1901.2016.tmp.dat.nc

# You can also click on Session/Working directory/Change directory and navigate
# to the folder where the temperature data are saved

# Load the CRU TS datasets into R 
tmp <- brick("data/input/cru_ts4.01.1901.2016.tmp.dat.nc", varname="tmp", package = "raster") # Mean monthly temperature

# Import sample site information
load("data/input/biotime.RData")

samples <- bt %>%
  filter(realm == "Terrestrial") %>%
  dplyr::select(timeseries_id, long, lat) %>%
  distinct()

samples_simple <- samples %>% dplyr::select(long, lat)
colnames(samples_simple) <- c("lon", "lat")

# Extract climate data from the RasterBrick as a data.frame
tmp.sites <- data.frame(raster::extract(tmp, samples_simple, ncol = 2)) # Mean monthly temperature

# Add sample site names to the data.frame
tmp.sites$timeseries_id <- samples$timeseries_id

# Change column names
years <- 1901:2016
month <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
names(tmp.sites) <- paste(rep(years, each=12), rep(month, times=116), sep="_")
names(tmp.sites)[1393] <- "rarefyID"

# Save the extracted climate data to a .csv file
save(tmp.sites, file = "data/output/CRU_BioTIME2022.RData")

# Turn into long format
tmp.sites.long <- tmp.sites %>% gather(year, temperature, c(1:1392))

# Calculate yearly averages
tmp.sites.long$year <- parse_number(tmp.sites.long$year)
tmp.sites.long <- tmp.sites.long %>% dplyr::group_by(rarefyID, year) %>%
  dplyr::summarise(mean_temp = mean(temperature))
tmp_sites_long <- tmp.sites.long

save(tmp_sites_long, file = "data/output/CRU_BioTIME_mean2022.RData")

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
  filter(realm == "Terrestrial") %>%
  dplyr::select(timeseries_id, long, lat) %>%
  distinct() %>% ungroup()

samples_lpd_simple <- samples_lpd %>% dplyr::select(long, lat)
colnames(samples_lpd_simple) <- c("lon", "lat")

# Extract climate data from the RasterBrick as a data.frame
tmp.sites.lpd <- data.frame(raster::extract(tmp, samples_lpd_simple, ncol = 2)) # Mean monthly temperature

# Change column names
names(tmp.sites.lpd) <- paste(rep(years, each = 12), rep(month, times = 116), sep = "_")

# Add sample site names to the data.frame
tmp.sites.lpd$timeseries_id <- samples_lpd$timeseries_id
names(tmp.sites.lpd)[1393] <- "timeseries_id"

# Save the extracted climate data to a .csv file
save(tmp.sites.lpd, file = "data/output/CRU_LPD2022.RData")

# Turn into long format
tmp.sites.long.lpd <- tmp.sites.lpd %>% gather(year, temperature, c(1:1392))

# Calculate yearly averages
tmp.sites.long.lpd$year <- parse_number(tmp.sites.long.lpd$year)
tmp.sites.long.lpd <- tmp.sites.long.lpd %>% dplyr::group_by(timeseries_id, year) %>%
  dplyr::summarise(mean_temp = mean(temperature))
tmp_sites_long_lpd <- tmp.sites.long.lpd

save(tmp_sites_long_lpd, file = "data/output/CRU_LPD_mean2022.RData")

