# Gergana Daskalova
# May 2021
# gndaskalova@gmail.com

# The goal of this script is to run statistical models
# comparing global change intensity across the sites
# represented in the Living Planet, BioTIME and PREDICTS databases
# as well as comparing them with global change across randomly sampled 
# sites around the world

# Libraries
library(tidyverse)
library(CoordinateCleaner)
library(rstan)
library(brms)
library(modelr)
library(bayesplot)
library(tidybayes)
#install_github("kassambara/factoextra")
library(factoextra)
library(sjstats)
library(parameters)
library(sjPlot)
library(ggeffects)
library(rgdal)
library(sp)

# Load population, biodiversity and driver data
load("data/output/popbio2022.RData") # Living Planet and BioTIME databases

# Space for time community data from PREDICTS
# Note this file is not on GitHub because of its size
# The file can be downloaded here 
# https://data.nhm.ac.uk/dataset/the-2016-release-of-the-predicts-database
# Choose Database in zipped CSV format
# The file is originally called database.csv, I renamed it to predicts.csv
# The file path below needs to be updated if the file is put in another location
predicts <- read.csv("data/input/predicts_sites.csv")
load("data/output/predicts_drivers.RData")
predicts_drivers <- distinct(predicts_drivers)
predicts_drivers$sampling <- "PREDICTS"
predicts_drivers$sampling2 <- "Terrestrial PREDICTS"
predicts_drivers$realm <- "Terrestrial"
load("data/output/random_samplingNov2020.RData")

# Get coordinates for database locations and their driver intensities
# LPD
lpd.coords <- popbio %>%
  filter(type == "Population")

lpd.coords.marine <- filter(lpd.coords, realm == "Marine")
lpd.coords.terr <- filter(lpd.coords, realm == "Terrestrial")

# BioTIME
bt.coords <- popbio %>%
  filter(type == "Biodiversity")

bt.coords.marine <- filter(bt.coords, realm == "Marine")
bt.coords.terr <- filter(bt.coords, realm == "Terrestrial")

# PREDICTS
predicts.coords <- predicts %>% 
  dplyr::select(Source_ID, Latitude, Longitude) %>%
  group_by(Source_ID) %>%
  summarise_all(mean) %>% dplyr::select(-Source_ID) %>%
  distinct()

# Create column with marine versus terrestrial
names(random_drivers)[7:8] <- c("decimallongitude", "decimallatitude")
shape <- readOGR(dsn = "data/input/ne_110m_land", layer = "ne_110m_land")
random_drivers_marine <- random_drivers %>% cc_sea(ref = shape)
random_drivers_marine$realm <- "Marine"
random_drivers_terr <- anti_join(random_drivers, random_drivers_marine)
random_drivers_terr$realm <- "Terrestrial"
random_drivers <- bind_rows(random_drivers_terr, random_drivers_marine)

lpd.coords$sampling <- "Living Planet Database"
lpd.coords$sampling2 <- paste0(lpd.coords$realm, lpd.coords$sampling)
random_drivers$sampling <- "Random sampling"
random_drivers$sampling2 <- as.factor(paste0(random_drivers$realm, random_drivers$sampling))

colnames(random_drivers)
colnames(lpd.coords)

lpd.coords <- lpd.coords %>% 
  ungroup() %>%
  dplyr::select(timeseries_id, climate_change, human_use,
                human_population, pollution, 
                invasions, cumulative,
                realm, sampling2)

combined_pop <- rbind(random_drivers[, c(1:6, 9, 11)], lpd.coords[, c(2:9)])

combined_pop$sampling2 <- factor(combined_pop$sampling2, 
                                 levels = c("MarineRandom sampling", "TerrestrialRandom sampling",
                                            "MarineLiving Planet Database", "TerrestrialLiving Planet Database"),
                                 labels = c("Marine Random sampling", "Terrestrial Random sampling",
                                            "Marine Living Planet Database", "Terrestrial Living Planet Database"))

combined_pop <- na.omit(combined_pop)
combined_pop_marine <- combined_pop %>% filter(realm == "Marine")
combined_pop_terr <- combined_pop %>% filter(realm == "Terrestrial")

# Standardise drivers between 0 and 1
predicts_drivers <- na.omit(predicts_drivers)
predicts_drivers$climate_change <- (predicts_drivers$climate_change - min(predicts_drivers$climate_change))/(max(predicts_drivers$climate_change)-min(predicts_drivers$climate_change))
predicts_drivers$cumulative <- (predicts_drivers$cumulative - min(predicts_drivers$cumulative))/(max(predicts_drivers$cumulative)-min(predicts_drivers$cumulative))
predicts_drivers$human_use <- (predicts_drivers$human_use - min(predicts_drivers$human_use))/(max(predicts_drivers$human_use)-min(predicts_drivers$human_use))
predicts_drivers$pollution <- (predicts_drivers$pollution - min(predicts_drivers$pollution))/(max(predicts_drivers$pollution)-min(predicts_drivers$pollution))

lpd.coords <- na.omit(lpd.coords)
lpd.coords$climate_change <- (lpd.coords$climate_change - min(lpd.coords$climate_change))/(max(lpd.coords$climate_change)-min(lpd.coords$climate_change))
lpd.coords$cumulative <- (lpd.coords$cumulative - min(lpd.coords$cumulative))/(max(lpd.coords$cumulative)-min(lpd.coords$cumulative))
lpd.coords$human_use <- (lpd.coords$human_use - min(lpd.coords$human_use))/(max(lpd.coords$human_use)-min(lpd.coords$human_use))
lpd.coords$pollution <- (lpd.coords$pollution - min(lpd.coords$pollution))/(max(lpd.coords$pollution)-min(lpd.coords$pollution))

bt.coords <- na.omit(bt.coords)
bt.coords$climate_change <- (bt.coords$climate_change - min(bt.coords$climate_change))/(max(bt.coords$climate_change)-min(bt.coords$climate_change))
bt.coords$cumulative <- (bt.coords$cumulative - min(bt.coords$cumulative))/(max(bt.coords$cumulative)-min(bt.coords$cumulative))
bt.coords$human_use <- (bt.coords$human_use - min(bt.coords$human_use))/(max(bt.coords$human_use)-min(bt.coords$human_use))
bt.coords$pollution <- (bt.coords$pollution - min(bt.coords$pollution))/(max(bt.coords$pollution)-min(bt.coords$pollution))

predicts_drivers$database <- "PREDICTS"
predicts_drivers <- predicts_drivers %>% gather(driver, intensity, 1:6)

lpd.coords$database <- "LPD"
lpd.coords <- lpd.coords %>% gather(driver, intensity, 2:7)
lpd.coords$sampling <- "LPD"
lpd.coords$sampling2 <- as.factor(paste0(lpd.coords$realm, lpd.coords$sampling))

bt.coords$database <- "BioTIME"
bt.coords <- bt.coords %>% gather(driver, intensity, 11:16)
bt.coords$sampling <- "BioTIME"
bt.coords$sampling2 <- as.factor(paste0(bt.coords$realm, bt.coords$sampling))

colnames(predicts_drivers)
colnames(lpd.coords)
colnames(bt.coords)
colnames(lpd.coords)[1] <- "study_id"
colnames(bt.coords)[2] <- "study_id"

lpd.coords.simple <- lpd.coords %>% dplyr::select(study_id, sampling, 
                                                  sampling2, realm, database, 
                                                  driver, intensity)
bt.coords.simple <- bt.coords %>% dplyr::select(study_id, sampling, 
                                                sampling2, realm, database, 
                                                driver, intensity)

str(predicts_drivers)
str(bt.coords.simple)
str(lpd.coords.simple)
bt.coords.simple$study_id <- as.factor(as.character(bt.coords.simple$study_id))
lpd.coords.simple$study_id <- as.factor(as.character(lpd.coords.simple$study_id))

drivers_combined <- bind_rows(predicts_drivers, lpd.coords.simple, bt.coords.simple)
summary(as.factor(drivers_combined$sampling2))

# Combine with random drivers
colnames(random_drivers)
colnames(drivers_combined)

# Scale random sampling too
random_drivers <- na.omit(random_drivers)
random_drivers$climate_change <- (random_drivers$climate_change - min(random_drivers$climate_change))/(max(random_drivers$climate_change)-min(random_drivers$climate_change))
random_drivers$cumulative <- (random_drivers$cumulative - min(random_drivers$cumulative))/(max(random_drivers$cumulative)-min(random_drivers$cumulative))
random_drivers$human_use <- (random_drivers$human_use - min(random_drivers$human_use))/(max(random_drivers$human_use)-min(random_drivers$human_use))
random_drivers$pollution <- (random_drivers$pollution - min(random_drivers$pollution))/(max(random_drivers$pollution)-min(random_drivers$pollution))

random_drivers_simple <- random_drivers %>% gather(driver, intensity, 1:6)
random_drivers_simple$study_id <- "NA"
random_drivers_simple$database <- "Random sampling"
colnames(random_drivers_simple)
random_drivers_simple <- random_drivers_simple %>% dplyr::select(study_id, sampling, 
                                                                 sampling2, realm, database, 
                                                                 driver, intensity)

drivers_combined_r <- bind_rows(drivers_combined, random_drivers_simple)

drivers_combined_r <- drivers_combined_r %>%
  mutate(linetype = as.factor(case_when(sampling2 == "Random sampling" ~ "random",
                                        sampling2 == "BioTIME" ~ "data",
                                        sampling2 == "LPD" ~ "data",
                                        sampling2 == "PREDICTS" ~ "data")))

random_sampling_only <- drivers_combined_r %>% filter(database == "Random sampling")
random_sampling_only$database <- "BioTIME"
random_sampling_only2 <- drivers_combined_r %>% filter(database == "Random sampling")
random_sampling_only2$database <- "LPD"

drivers_combined_r <- drivers_combined_r %>% filter(database != "Random sampling")

drivers_combined_r <- bind_rows(drivers_combined_r, random_sampling_only,
                                random_sampling_only2)

drivers_combined_r$driver <- factor(drivers_combined_r$driver,
                                    levels = c("human_use", "climate_change",
                                               "human_population", "pollution",
                                               "invasions", "cumulative"),
                                    labels = c("Human use", "Climate change",
                                               "Human population", "Pollution",
                                               "Invasions", "Cumulative"))

random_drivers_terr_scaled <- random_drivers_simple %>% filter(realm == "Terrestrial" &
                                                                 driver == "human_use")
min(random_drivers_terr_scaled$intensity)
max(random_drivers_terr_scaled$intensity)

predicts_drivers_scaled <- predicts_drivers %>% filter(driver == "human_use")
min(predicts_drivers_scaled$intensity)
max(predicts_drivers_scaled$intensity)

drivers_combined_r_pr_r <- drivers_combined_r %>% filter(sampling2 == "TerrestrialRandom sampling") 
drivers_combined_r_pr_r$database <- "PREDICTS"

drivers_combined_r <- bind_rows(drivers_combined_r, drivers_combined_r_pr_r)
#save(drivers_combined_r, file = "data/output/drivers_combined_r2022.RData")

# Models ----
#load("data/output/drivers_combined_r2022.RData")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Randomly select only a portion of the random sampling points to make the sample sizes more equal
set.seed(123)
random_drivers_simple_select <- sample_n(random_drivers_simple, 10000)

drivers_combined_r_r <- bind_rows(drivers_combined, random_drivers_simple_select)

# Set priors
prior <- c(set_prior(prior = 'normal(0,6)', class='Intercept', coef=''))	

summary(as.factor(drivers_combined_r_r$sampling2))
drivers_combined_r_r_marine <- drivers_combined_r_r %>% filter(realm == "Marine")
drivers_combined_r_r_terr <- drivers_combined_r_r %>% filter(realm == "Terrestrial")

drivers_combined_r_r_terr$sampling2 <- factor(drivers_combined_r_r_terr$sampling2,
                                              levels = c("TerrestrialRandom sampling",
                                                         "TerrestrialLPD",
                                                         "TerrestrialBioTIME",
                                                         "Terrestrial PREDICTS"),
                                              labels = c("Random sampling",
                                                         "Living Planet",
                                                         "BioTIME",
                                                         "PREDICTS"))

terr_driver_m <- brm(bf(intensity ~ sampling2*driver), 
                     data = drivers_combined_r_r_terr, 
                     prior = prior, iter = 6000,
                     warmup = 2000,
                     inits = '0',
                     control = list(adapt_delta = 0.80),
                     cores = 2, chains = 2)

# Check model and save output
summary(terr_driver_m)
plot(terr_driver_m)
save(terr_driver_m, file = "data/output/terr_driver_m2022.RData")

drivers_combined_r_r_marine$sampling2 <- factor(drivers_combined_r_r_marine$sampling2,
                                                levels = c("MarineRandom sampling",
                                                           "MarineLPD",
                                                           "MarineBioTIME"),
                                                labels = c("Random sampling",
                                                           "Living Planet",
                                                           "BioTIME"))

mar_driver_m <- brm(bf(intensity ~ sampling2*driver), 
                    data = drivers_combined_r_r_marine, 
                    prior = prior, iter = 6000,
                    warmup = 2000,
                    inits = '0',
                    control = list(adapt_delta = 0.80),
                    cores = 2, chains = 2)

# Check model and save output
summary(mar_driver_m)
plot(mar_driver_m)
save(mar_driver_m, file = "data/output/mar_driver_m2022.RData")
