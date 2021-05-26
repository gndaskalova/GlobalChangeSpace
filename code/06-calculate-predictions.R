# Gergana Daskalova
# May 2021
# gndaskalova@gmail.com

# The goal of this script is to calculate model predictions for how 
# global change intensity varies across sites represented by the Living Planet,
# BioTIME and PREDICTS databases and across the world

# Libraries
library(tidyverse)
library(rstan)
library(brms)
library(modelr)
library(bayesplot)
library(tidybayes)
library(sjstats)
library(parameters)
library(sjPlot)
library(ggeffects)

# Load models, output from from 05-run-models.R
load("data/output/terr_driver_m.RData")
load("data/output/mar_driver_m.RData")

# Calculate model predictions
predictions_terr <- ggpredict(terr_driver_m, terms = c("driver", "sampling2"))

# Note that a warning message appears about a deprecated function
# Warning message:
# posterior_linpred(transform = TRUE) is deprecated. Please use pp_expect() instead, without the 'transform' argument. 

# The code still works

predictions_random <- predictions_terr %>% filter(group == "Random sampling")
colnames(predictions_random) <- c("x", "random_predicted",
                                  "random_conf.low",
                                  "random_conf.high", "group")

predictions_random2 <- bind_rows(predictions_random,
                                 predictions_random,
                                 predictions_random)
predictions_random2 <- arrange(predictions_random2, x)
predictions_random2 <- predictions_random2 %>% dplyr::select(-group, -x)
predictions_terr <- predictions_terr %>% filter(group != "Random sampling")
predictions_terr <- arrange(predictions_terr, x)
predictions_terr <- bind_cols(predictions_terr, predictions_random2)
predictions_terr$deviation <- predictions_terr$predicted - predictions_terr$random_predicted
predictions_terr$deviation_lower <- predictions_terr$conf.low - predictions_terr$random_conf.low
predictions_terr$deviation_upper <- predictions_terr$conf.high - predictions_terr$random_conf.high

predictions_mar <- ggpredict(mar_driver_m, terms = c("driver", "sampling2"))
predictions_random <- predictions_mar %>% filter(group == "Random sampling")
colnames(predictions_random) <- c("x", "random_predicted",
                                  "random_conf.low",
                                  "random_conf.high", "group")

predictions_random2 <- bind_rows(predictions_random,
                                 predictions_random)
predictions_random2 <- arrange(predictions_random2, x)
predictions_random2 <- predictions_random2 %>% dplyr::select(-group, -x)
predictions_mar <- predictions_mar %>% filter(group != "Random sampling")
predictions_mar <- arrange(predictions_mar, x)
predictions_mar <- bind_cols(predictions_mar, predictions_random2)
predictions_mar$deviation <- predictions_mar$predicted - predictions_mar$random_predicted
predictions_mar$deviation_lower <- predictions_mar$conf.low - predictions_mar$random_conf.low
predictions_mar$deviation_upper <- predictions_mar$conf.high - predictions_mar$random_conf.high
predictions_mar$realm <- "Marine"
predictions_terr$realm <- "Terrestrial"

predictions_combo <- bind_rows(predictions_terr, predictions_mar)
predictions_combo$data_realm <- paste0(predictions_combo$group, predictions_combo$realm)
predictions_combo$data_realm <- factor(predictions_combo$data_realm,
                                       levels = c("Living PlanetTerrestrial",
                                                  "BioTIMETerrestrial",
                                                  "PREDICTSTerrestrial",
                                                  "Living PlanetMarine",
                                                  "BioTIMEMarine"),
                                       labels = c("Living PlanetTerrestrial",
                                                  "BioTIMETerrestrial",
                                                  "PREDICTSTerrestrial",
                                                  "Living PlanetMarine",
                                                  "BioTIMEMarine"))

save(predictions_combo, file = "data/output/predictions_combo.RData")