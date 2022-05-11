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
                realm, sampling2) %>% distinct()

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

# Sampled versus global global change space ----

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

# Check sample sizes
samples_terr <- drivers_combined_r_r_terr %>%
  group_by(sampling2) %>%
  summarise(n = n())

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

# Temperature over time before and during time series ----
load("data/output/CRU_BioTIME_mean2022.RData")
load("data/output/CRU_LPD_mean2022.RData")
load("data/output/NOAA_BioTIME_mean2022.RData")
load("data/output/NOAA_LPD_mean2022.RData")

# Filter out the temperature records for periods
# equal to the duration of each time-series (before they started)

# e.g. if a time-series is from 1980 to 1990, the matching before period
# is 1970 to 1980
# Choose time-series with at least 10 years duration so that a temperature
# trend can be estimated
popbio_long <- popbio %>% filter(duration > 9)

bt <- popbio_long %>% filter(type == "Biodiversity")

bt$start_comparison <- bt$start_year - bt$duration
bt$end_comparison <- bt$start_year

bt_simpler <- bt %>% dplyr::select(timeseries_id, start_year, end_year, start_comparison,
                                   end_comparison, duration)
colnames(bt_simpler)[1] <- "rarefyID"

temp_bt <- left_join(tmp_sites_long, bt_simpler, by = "rarefyID")

str(temp_bt)

temp_bt2 <- temp_bt %>% dplyr::filter(year > start_comparison -1)
temp_bt3 <- temp_bt2 %>% dplyr::filter(year < end_year + 1)
temp_bt3 <- temp_bt3 %>% mutate(timing = case_when(year < end_comparison + 1 ~ "Before monitoring",
                                                   year > start_year - 1 ~ "During monitoring"))

temp_bt3 <- na.omit(temp_bt3)

# Set priors
hier_prior <- c(set_prior(prior = 'normal(0,6)', class='b', coef='year.scaled'), 	# global slope
                set_prior(prior = 'normal(0,6)', class='Intercept', coef=''), 		# global intercept
                set_prior(prior = 'cauchy(0,2)', class='sd'))							# group-level intercepts and slopes

# Temperature over time in the period before monitoring started
temp_bt3_bf <- temp_bt3 %>% filter(timing == "Before monitoring")

summary(temp_bt3_bf$year)  # starts in 1901, ends in 2016
temp_bt3_bf$year2 <- temp_bt3_bf$year - 1900  # make years go 1, 2, 3, 4...
temp_bt3_bf$year.scaled <- scale(temp_bt3_bf$year2, center = T)
summary(temp_bt3_bf$year.scaled)  # centered on zero

# One model for all time series - then I extract the slopes of the random effect
temp_m_before <- brm(bf(mean_temp ~ year.scaled + (year.scaled|rarefyID)), 
                      data = temp_bt3_bf, iter = 12000,
                   warmup = 2000, prior = hier_prior,
                   init = '0',
                   control = list(adapt_delta = 0.89),
                   cores = 4, chains = 4)

summary(temp_m_before)
plot(temp_m_before)
save(temp_m_before, file = "data/output/richness_model_all.RData")

# Extract slopes for each cell
slopes3 <- as.data.frame(coef(temp_m_before))
save(slopes3, file = "data/output/slopes_richness.RData")

hist(slopes3$rarefyID.Estimate.year.scaled)


temp_bt3 <- temp_bt3 %>% mutate(timing = case_when(year < end_comparison + 1 ~ "Before monitoring",
                                                   year > start_year - 1 ~ "During monitoring"))

temp_bt3b <- na.omit(temp_bt3)

temp_bt_before <- temp_bt3b %>% dplyr::filter(timing == "Before monitoring")
temp_bt_during <- temp_bt3b %>% dplyr::filter(timing == "During monitoring")



temp_models <- temp_bt3b %>% group_by(rarefyID, timing) %>%
  do(broom::tidy(lm(mean_temp ~ year, data = .)))

temp_models <- temp_models %>%
  spread(term, estimate)

temp_models_spread <- temp_models %>% 
  drop_na(year) %>%
  dplyr::select(rarefyID, timing, year) %>%
  spread(timing, year)

# Add labels for the squares, add degrees per year
(bt_climate1 <- ggplot(temp_models_spread, aes(y = `Before monitoring`, x = `During monitoring`)) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon",
                    alpha = 0.8) +
    labs(x = "\nTemperature change\nduring monitoring (slope)",
         y = "Temperature change\nbefore monitoring (slope)\n", fill = "Density") +
    scale_fill_gradient(low = "#eaa9bc", high = "#441540") +
    ggtitle("\nd     Climate warming over time (BioTIME terrestrial)") +
    coord_cartesian(y = c(-0.12, 0.12), xlim = c(-0.12, 0.12)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme_classic() +
    theme(legend.position = c(0.12, 0.8),
          legend.title = element_text(color = "black", size = 10),
          text = element_text(color = "#22211d"),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.line.x = element_line(size = 0.3),
          axis.line.y = element_line(size = 0.3),
          axis.ticks.length = unit( .1, "cm"),
          #plot.margin = unit(c(2, 2, 2, 2),"cm"),
          plot.title = element_text(size = 14, hjust = 0.05, 
                                    color = "grey20", face = "bold")))

ggsave(bt_climate1, filename = "figures/biotime_cru_terr2022.png", height = 5, width = 5)
ggsave(bt_climate1, filename = "figures/biotime_cru_terr2022.pdf", height = 5, width = 5)

temp_models_spread$test <- temp_models_spread$`During monitoring` >  temp_models_spread$`Before monitoring`
summary(temp_models_spread$test)
# Mode   FALSE    TRUE    NA's 
# logical     667    1269       2 

1269/(1269+667)

# 65% of BioTIME terr time series experienced more warming during time series than before

# LPD
lpd <- popbio %>% filter(type == "Population")

lpd$start_comparison <- lpd$start_year - lpd$duration
lpd$end_comparison <- lpd$start_year - 1

lpd_simpler <- lpd %>% dplyr::select(timeseries_id, start_year, end_year, start_comparison,
                                     end_comparison, duration)

str(tmp_sites_long_lpd)
str(lpd_simpler)
# Make timeseries_id columns the same type
tmp_sites_long_lpd$timeseries_id <- as.factor(as.character(tmp_sites_long_lpd$timeseries_id))
lpd_simpler$timeseries_id <- as.factor(as.character(lpd_simpler$timeseries_id))

temp_lpd <- left_join(tmp_sites_long_lpd, lpd_simpler, by = "timeseries_id")

str(temp_lpd)

temp_lpd2 <- temp_lpd %>% dplyr::filter(year > start_comparison -1)
temp_lpd3 <- temp_lpd2 %>% dplyr::filter(year < end_year + 1)

temp_lpd3 <- temp_lpd3 %>% mutate(timing = case_when(year < end_comparison + 1 ~ "Before monitoring",
                                                     year > start_year - 1 ~ "During monitoring"))

temp_lpd3b <- na.omit(temp_lpd3)

# Keep only time series with duration of more than 4 years
temp_lpd3b <- temp_lpd3b %>%
  filter(duration > 4)

temp_models_lpd <- temp_lpd3b %>% group_by(timeseries_id, timing) %>%
  do(broom::tidy(lm(mean_temp ~ year, data = .)))

temp_models_lpd <- temp_models_lpd %>%
  spread(term, estimate)

lpd_temp_models_spread <- temp_models_lpd %>% drop_na(year) %>%
  dplyr::select(timeseries_id, timing, year) %>%
  spread(timing, year)

(lpd_climate1 <- ggplot(lpd_temp_models_spread, aes(y = `Before monitoring`, x = `During monitoring`)) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon",
                    alpha = 0.8) +
    labs(x = "\nTemperature change\nduring monitoring (slope)",
         y = "Temperature change\nbefore monitoring (slope)\n", fill = "Density") +
    scale_fill_gradient(low = "#be925a", high = "#36230d") +
    ggtitle("\nc     Climate warming over time (Living Planet)") +
    coord_cartesian(y = c(-0.12, 0.12), xlim = c(-0.12, 0.12)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme_classic() +
    theme(legend.position = c(0.12, 0.8),
          legend.title = element_text(color = "black", size = 10),
          text = element_text(color = "#22211d"),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.line.x = element_line(size = 0.3),
          axis.line.y = element_line(size = 0.3),
          axis.ticks.length = unit( .1, "cm"),
          #plot.margin = unit(c(2, 2, 2, 2),"cm"),
          plot.title = element_text(size = 14, hjust = 0.05, 
                                    color = "grey20", face = "bold")))

warming_terr_panel <- grid.arrange(lpd_climate1, bt_climate1, ncol = 2)

ggsave(warming_terr_panel, filename = "figures/warming_panel2022.png", height = 5, width = 12)
ggsave(warming_terr_panel, filename = "figures/warming_panel2022.pdf", height = 5, width = 12)

lpd_temp_models_spread$test <- lpd_temp_models_spread$`During monitoring` >  lpd_temp_models_spread$`Before monitoring`
summary(lpd_temp_models_spread$test)
# Mode   FALSE    TRUE 
# logical    1298    2246 

2246/(2246+1298)

# 63% of LPD terr time series experienced more warming during time series than before

# Marine LPD
# Make timeseries_id columns the same type
sst_sites_long_lpd$timeseries_id <- as.factor(as.character(sst_sites_long_lpd$timeseries_id))

temp_lpd_mar <- left_join(sst_sites_long_lpd, lpd_simpler, by = "timeseries_id")

str(temp_lpd_mar)

temp_lpd_mar2 <- temp_lpd_mar %>% dplyr::filter(year > start_comparison - 1)
temp_lpd_mar3 <- temp_lpd_mar2 %>% dplyr::filter(year < end_year + 1)

temp_lpd_mar3 <- temp_lpd_mar3 %>% mutate(timing = case_when(year < end_comparison + 1 ~ "Before monitoring",
                                                             year > start_year - 1 ~ "During monitoring"))

temp_lpd_mar3b <- na.omit(temp_lpd_mar3)

# Keep only time series with more than 4 years of data
temp_lpd_mar3b <- temp_lpd_mar3b %>%
  filter(duration > 4) %>% distinct()

temp_models_lpd_mar <- temp_lpd_mar3b %>% group_by(timeseries_id, timing) %>%
  do(broom::tidy(lm(mean_temp ~ year, data = .)))

temp_models_lpd_mar <- temp_models_lpd_mar %>%
  spread(term, estimate)

lpd_temp_models_spread_mar <- temp_models_lpd_mar %>% drop_na(year) %>%
  dplyr::select(timeseries_id, timing, year) %>%
  spread(timing, year) %>%
  drop_na(`Before monitoring`)

# Note that there are some time series for which climate data were not available
# (too far back in the past) thus why there were time series for which temperature
# change before the monitoring began could not be calculated

(lpd_climate2 <- ggplot(lpd_temp_models_spread_mar, aes(y = `Before monitoring`, x = `During monitoring`)) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon",
                    alpha = 0.8) +
    labs(x = "\nTemperature change\nduring monitoring (slope)",
         y = "Temperature change\nbefore monitoring (slope)\n", fill = "Density") +
    scale_fill_gradient(low = "#50b2b2", high = "#235051") +
    ggtitle("\ne     Climate warming over time (Living Planet marine)") +
    #coord_cartesian(y = c(-0.12, 0.12), xlim = c(-0.12, 0.12)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme_classic() +
    theme(legend.position = c(0.12, 0.8),
          legend.title = element_text(color = "black", size = 10),
          text = element_text(color = "#22211d"),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.line.x = element_line(size = 0.3),
          axis.line.y = element_line(size = 0.3),
          axis.ticks.length = unit( .1, "cm"),
          #plot.margin = unit(c(2, 2, 2, 2),"cm"),
          plot.title = element_text(size = 14, hjust = 0.05, 
                                    color = "grey20", face = "bold")))

lpd_temp_models_spread_mar$test <- lpd_temp_models_spread_mar$`During monitoring` >  lpd_temp_models_spread_mar$`Before monitoring`
summary(lpd_temp_models_spread_mar$test)
# Mode   FALSE    TRUE 
# logical     573     465 

465/(465+573)

# 45% of LPD marine time series experienced more warming during time series than before

# Marine BioTIME
colnames(sst_sites_long)[2] <- "year"
sst_sites_long$year <- parse_number(sst_sites_long$year)

temp_bt_mar <- left_join(sst_sites_long, bt_simpler, by = "rarefyID")

str(temp_bt_mar)

temp_bt_mar2 <- temp_bt_mar %>% dplyr::filter(year > start_comparison -1)
temp_bt_mar3 <- temp_bt_mar2 %>% dplyr::filter(year < end_year + 1)

temp_bt_mar3 <- temp_bt_mar3 %>% mutate(timing = case_when(year < end_comparison + 1 ~ "Before monitoring",
                                                           year > start_year - 1 ~ "During monitoring"))

temp_bt_mar3b <- na.omit(temp_bt_mar3)
temp_bt_mar3b <- distinct(temp_bt_mar3)
temp_bt_mar3b <- temp_bt_mar3b %>% drop_na(temperature)

temp_models_bt_mar <- temp_bt_mar3b %>% group_by(rarefyID, timing) %>%
  do(broom::tidy(lm(temperature ~ year, data = .)))

temp_models_bt_mar <- temp_models_bt_mar %>%
  spread(term, estimate)

bt_temp_models_spread_mar <- temp_models_bt_mar %>% drop_na(year) %>%
  dplyr::select(rarefyID, timing, year) %>%
  spread(timing, year)

(bt_climate2 <- ggplot(bt_temp_models_spread_mar, aes(y = `Before monitoring`, x = `During monitoring`)) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon",
                    alpha = 0.8) +
    labs(x = "\nTemperature change\nduring monitoring (slope)",
         y = "Temperature change\nbefore monitoring (slope)\n", fill = "Density") +
    scale_fill_gradient(low = "#50b2b2", high = "#235051") +
    ggtitle("\nb     Climate warming over time (Living Planet marine)") +
    #coord_cartesian(y = c(-0.12, 0.12), xlim = c(-0.12, 0.12)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme_classic() +
    theme(legend.position = c(0.12, 0.8),
          legend.title = element_text(color = "black", size = 10),
          text = element_text(color = "#22211d"),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.line.x = element_line(size = 0.3),
          axis.line.y = element_line(size = 0.3),
          axis.ticks.length = unit( .1, "cm"),
          #plot.margin = unit(c(2, 2, 2, 2),"cm"),
          plot.title = element_text(size = 14, hjust = 0.05, 
                                    color = "grey20", face = "bold")))

bt_temp_models_spread_mar$test <- bt_temp_models_spread_mar$`During monitoring` >  bt_temp_models_spread_mar$`Before monitoring`
summary(bt_temp_models_spread_mar$test)

# Mode   FALSE    TRUE 
# logical    3013    4388 

4388/(4388+3013)

# 0.5928929

# 59% of BioTIME mar time series experienced more warming during time series than before

warming_terr_mar_panel <- grid.arrange(lpd_climate1, bt_climate1, lpd_climate2, bt_climate2, ncol = 4)
warming_terr_mar_panel <- grid.arrange(lpd_climate1, bt_climate1, lpd_climate2, ncol = 3)

ggsave(warming_terr_mar_panel, filename = "figures/warming_panel3.png", height = 6, width = 19)
ggsave(warming_terr_mar_panel, filename = "figures/warming_panel3.pdf", height = 6, width = 19)

# Model testing climate warming before and during time series moniitoring

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set priors
prior <- c(set_prior(prior = 'normal(0,6)', class='Intercept', coef=''))	

# turn data into long form
bt_temp_models_long_mar <- bt_temp_models_spread_mar %>%
  gather(period, estimate, 2:3)

bt_mar_temp_m <- brm(bf(estimate ~ period), 
                     data = bt_temp_models_long_mar, 
                     prior = prior, iter = 6000,
                     warmup = 2000,
                     inits = '0',
                     control = list(adapt_delta = 0.80),
                     cores = 2, chains = 2)

# Check model and save output
summary(bt_mar_temp_m)
plot(bt_mar_temp_m)
save(bt_mar_temp_m, file = "data/output/bt_mar_temp_m.RData")

# turn data into long form
lpd_temp_models_long_mar <- lpd_temp_models_spread_mar %>%
  gather(period, estimate, 2:3)

lpd_mar_temp_m <- brm(bf(estimate ~ period), 
                      data = lpd_temp_models_long_mar, 
                      prior = prior, iter = 6000,
                      warmup = 2000,
                      inits = '0',
                      control = list(adapt_delta = 0.80),
                      cores = 2, chains = 2)

# Check model and save output
summary(lpd_mar_temp_m)
plot(lpd_mar_temp_m)
save(lpd_mar_temp_m, file = "data/output/lpd_mar_temp_m.RData")

# turn data into long form
bt_temp_models_long_terr <- temp_models_spread %>%
  gather(period, estimate, 2:3)

bt_terr_temp_m <- brm(bf(estimate ~ period), 
                      data = bt_temp_models_long_terr, 
                      prior = prior, iter = 6000,
                      warmup = 2000,
                      inits = '0',
                      control = list(adapt_delta = 0.80),
                      cores = 2, chains = 2)

# Check model and save output
summary(bt_terr_temp_m)
plot(bt_terr_temp_m)
save(bt_terr_temp_m, file = "data/output/bt_terr_temp_m.RData")

# turn data into long form
lpd_temp_models_long_terr <- lpd_temp_models_spread %>%
  gather(period, estimate, 2:3)

lpd_terr_temp_m <- brm(bf(estimate ~ period), 
                       data = lpd_temp_models_long_terr, 
                       prior = prior, iter = 6000,
                       warmup = 2000,
                       inits = '0',
                       control = list(adapt_delta = 0.80),
                       cores = 2, chains = 2)

# Check model and save output
summary(lpd_terr_temp_m)
plot(lpd_terr_temp_m)
save(lpd_terr_temp_m, file = "data/output/lpd_terr_temp_m.RData")
