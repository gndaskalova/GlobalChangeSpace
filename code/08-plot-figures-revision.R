library(tidyverse)
library(data.table)

# Data viz theme
change_theme <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica"),
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(0.82, 0.27), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

# Load population, biodiversity and driver data
load("data/output/popbio2022.RData") # Living Planet and BioTIME databases
popbio <- popbio %>% filter(biome != "NULL") # removing blank duplicates

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

# Colours
# #56b3b1 BioTIME terr
# #4d5f72 BioTIME mar
# #b6af3c LPD terr
# #72886e LPD mar
# #eb5a65 PREDICTS

# Fig 1 Global change spatial representation ----
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

combined_pop <- rbind(random_drivers[, c(1:6, 9, 11)], lpd.coords[,c(11:16, 3, 19)])

combined_pop$sampling2 <- factor(combined_pop$sampling2, 
                                 levels = c("MarineRandom sampling", "TerrestrialRandom sampling",
                                            "MarineLiving Planet Database", "TerrestrialLiving Planet Database"),
                                 labels = c("Marine Random sampling", "Terrestrial Random sampling",
                                            "Marine Living Planet Database", "Terrestrial Living Planet Database"))

combined_pop <- na.omit(combined_pop)
combined_pop_marine <- combined_pop %>% filter(realm == "Marine")
combined_pop_terr <- combined_pop %>% filter(realm == "Terrestrial")

combined_pop_marine_simple <- combined_pop_marine %>% 
  dplyr::select(-cumulative, - realm) %>% distinct()

pop_mar_pca <- prcomp(combined_pop_marine_simple[, -6], scale = TRUE)

(combined_pop_pca_mar <- factoextra::fviz_pca_biplot(pop_mar_pca, label ="var", 
                                                     select.var= list(name = c("human_use", "climate_change",
                                                                               "human_population")),
                                                     habillage = combined_pop_marine_simple$sampling2, col.var = "grey20",
                                                     addEllipses = TRUE, ellipse.level = 0.95,
                                                     repel = TRUE, pointsize = 2, arrowsize = 2,
                                                     labelsize = 0) +
    scale_colour_manual(values = c("grey80", "#419191")) +
    scale_fill_manual(values = c("grey80", "#419191")) +
    guides(colour = F, shape = F, fill = F) +
    theme_void() +
    annotate(geom = "text", x = -2.5, y = -2.8, label = "Climate change", size = 5,
             color="grey20") +
    #  annotate(geom = "text", x = -4.5, y = 0.2, label = "Human density", size = 5,
    #          color="grey20") +
    labs(title = NULL) +
    coord_cartesian(xlim = c(5, -5),
                    ylim = c(3, -3)))

combined_pop_terr_simple <- combined_pop_terr %>% 
  dplyr::select(-cumulative, - realm) %>%
  distinct()

pop_terr_pca <- prcomp(combined_pop_terr_simple[, 1:5], scale = TRUE)

(combined_pop_pca_terr <- fviz_pca_biplot(pop_terr_pca, label ="var", 
                                          select.var= list(name = c("human_use", "climate_change",
                                                                    "human_population")),
                                          habillage = combined_pop_terr_simple$sampling2,
                                          col.var = "grey20",
                                          addEllipses = TRUE, ellipse.level = 0.95,
                                          repel = TRUE, pointsize = 2, arrowsize = 2, labelsize = 0) +
    scale_colour_manual(values = c("grey80", "#9a6b36")) +
    scale_fill_manual(values = c("grey80", "#9a6b36")) +
    guides(colour = F, shape = F, fill = F) +
    theme_void() +
    annotate(geom = "text", x = -2.9, y = -3.5, label = "Climate change", size = 5,
             color="grey20") +
    # annotate(geom = "text", x = -6.62, y = 0.5, label = "Human density", size = 5,
    #         color="grey20") +
    labs(title = NULL) +
    coord_cartesian(xlim = c(7, -7),
                    ylim = c(3.5, -3.5)))

pop_pca_panel <- grid.arrange(combined_pop_pca_terr, combined_pop_pca_mar, ncol = 2)

ggsave(pop_pca_panel, filename = "figures/pop_pca_2022.pdf",
       height = 5, width = 10)

ggsave(pop_pca_panel, filename = "figures/pop_pca_2022.png",
       height = 5, width = 10)

bt.coords$sampling <- "BioTIME"
bt.coords$sampling2 <- paste0(bt.coords$realm, bt.coords$sampling)

combined_bio <- rbind(random_drivers[, c(1:6, 9, 11)], bt.coords[,c(11:16, 3, 19)])

combined_bio$sampling2 <- factor(combined_bio$sampling2, 
                                 levels = c("MarineRandom sampling", "TerrestrialRandom sampling",
                                            "MarineBioTIME", "TerrestrialBioTIME"),
                                 labels = c("Marine Random sampling", "Terrestrial Random sampling",
                                            "Marine BioTIME", "Terrestrial BioTIME"))
combined_bio <- na.omit(combined_bio)
combined_bio_marine <- combined_bio %>% filter(realm == "Marine")
combined_bio_terr <- combined_bio %>% filter(realm == "Terrestrial")

combined_bio_marine_simple <- combined_bio_marine %>% 
  dplyr::select(-cumulative, - realm) %>% distinct()

bio_mar_pca <- prcomp(combined_bio_marine_simple[, -6], scale = TRUE)

(combined_bio_pca_mar <- fviz_pca_biplot(bio_mar_pca, label ="var", 
                                         select.var= list(name = c("human_use", "climate_change",
                                                                   "human_population")),
                                         habillage = combined_bio_marine_simple$sampling2, col.var = "grey20",
                                         addEllipses = TRUE, ellipse.level = 0.95,
                                         repel = TRUE, pointsize = 2, arrowsize = 2, labelsize = 0) +
    scale_colour_manual(values = c("grey80", "#3b738f")) +
    scale_fill_manual(values = c("grey80", "#3b738f")) +
    guides(colour = F, shape = F, fill = F) +
    theme_void() +
    annotate(geom = "text", x = 2.5, y = -2.8, label = "Climate change", size = 5,
             color="grey20") +
    # annotate(geom = "text", x = 4.5, y = 0.2, label = "Human density", size = 5,
    #         color="grey20") +
    labs(title = NULL) +
    coord_cartesian(xlim = c(-5, 5),
                    ylim = c(3, -3)))

combined_bio_terr_simple <- combined_bio_terr %>% 
  dplyr::select(-cumulative, - realm) %>% distinct()

bio_terr_pca <- prcomp(combined_bio_terr_simple[, 1:5], scale = TRUE)

(combined_bio_pca_terr <- fviz_pca_biplot(bio_terr_pca, label ="var", 
                                          select.var= list(name = c("human_use", "climate_change",
                                                                    "human_population")),
                                          habillage = combined_bio_terr_simple$sampling2,
                                          col.var = "grey20",
                                          addEllipses = TRUE, ellipse.level = 0.95,
                                          repel = TRUE, pointsize = 2, arrowsize = 2, labelsize = 0) +
    scale_colour_manual(values = c("grey80", "#91367e")) +
    scale_fill_manual(values = c("grey80", "#91367e")) +
    guides(colour = F, shape = F, fill = F) +
    theme_void() +
    annotate(geom = "text", x = -2.9, y = -3.5, label = "Climate change", size = 5,
             color="grey20") +
    #  annotate(geom = "text", x = -6.62, y = 0.5, label = "Human density", size = 5,
    #          color="grey20") +
    labs(title = NULL) +
    coord_cartesian(xlim = c(7, -7),
                    ylim = c(3.5, -3.5)))

bio_pca_panel <- grid.arrange(combined_bio_pca_terr, combined_bio_pca_mar, ncol = 2)
ggsave(bio_pca_panel, filename = "figures/bio_pca_2022.pdf",
       height = 5, width = 10)
ggsave(bio_pca_panel, filename = "figures/bio_pca_2022.png",
       height = 5, width = 10)

combined_predicts <- rbind(random_drivers[,c(1:6, 9, 11)], predicts_drivers[,c(1:6, 10, 9)])

combined_predicts$sampling2 <- factor(combined_predicts$sampling2, 
                                      levels = c("TerrestrialRandom sampling", "Terrestrial PREDICTS"),
                                      labels = c("Terrestrial Random sampling", "Terrestrial PREDICTS"))

combined_predicts <- na.omit(combined_predicts)
combined_predicts <- filter(combined_predicts, realm == "Terrestrial")

combined_predicts_terr_simple <- combined_predicts %>% dplyr::select(-cumulative, - realm) %>%
  distinct()

predicts_terr_pca <- prcomp(combined_predicts_terr_simple[, 1:5], scale = TRUE)

(combined_predicts_pca_terr <- fviz_pca_biplot(predicts_terr_pca, label ="var", 
                                               select.var= list(name = c("human_use", "climate_change",
                                                                         "human_population")),
                                               habillage = combined_predicts_terr_simple$sampling2,
                                               col.var = "grey20",
                                               addEllipses = TRUE, ellipse.level = 0.95,
                                               repel = TRUE, pointsize = 2, arrowsize = 2, labelsize = 0) +
    scale_colour_manual(values = c("grey80", "#7d8141")) +
    scale_fill_manual(values = c("grey80", "#7d8141")) +
    guides(colour = F, shape = F, fill = F) +
    theme_void() +
    annotate(geom = "text", x = -2.9, y = -3.5, label = "Climate change", size = 5,
             color="grey20") +
    # annotate(geom = "text", x = -6.62, y = 0.5, label = "Human density", size = 5,
    #         color="grey20") +
    labs(title = NULL) +
    coord_cartesian(xlim = c(7, -7),
                    ylim = c(3.5, -3.5)))

ggsave(combined_predicts_pca_terr, filename = "figures/predicts_gc_space2022.pdf",
       height = 5, width = 5)

ggsave(combined_predicts_pca_terr, filename = "figures/predicts_gc_space2022.png",
       height = 5, width = 5)

rds_combo <- list.files(path = "data/output/HPC_hypervolumes", pattern = "*.rds", full.names = TRUE ) %>%
  map_dfr(readRDS)

all_paths <-
  list.files(path = "data/output/HPC_hypervolumes",
             pattern = "*.rds",
             full.names = TRUE)

all_content <-
  all_paths %>%
  lapply(readRDS)

all_filenames <- all_paths %>%
  basename() %>%
  as.list()

all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)

all_results <- rbindlist(all_lists, fill = T)
# change column name
names(all_results)[6] <- "path"

all_results <- all_results %>%
  separate(path, c("a", "b", "c", "d", "e", "f", "g"), sep = "_")

all_results$d <- as.factor(all_results$d)
all_results$e <- as.factor(all_results$e)
all_results$g <- as.factor(all_results$g)

mean_overlaps <- all_results %>%
  dplyr::group_by(d, e, g) %>%
  dplyr::summarise(mean_sorensen = mean(sorensen))

mean_overlaps$data <- paste(mean_overlaps$d, mean_overlaps$e, sep = "_")
mean_overlaps$g <- as.numeric(as.character(mean_overlaps$g))

mean_overlaps$data <- factor(mean_overlaps$data,
                             levels = c("pop_terr", "bio_terr", 
                                        "predicts_terr", "pop_marine",
                                        "bio_marine"),
                             labels = c("Living Planet (terrestrial)",
                                        "BioTIME (terrestrial)",
                                        "PREDICTS (terrestrial)",
                                        "Living Planet (marine)",
                                        "BioTIME (marine)"))

(overlap_accumulation <- ggplot(mean_overlaps, aes(x = g, y = mean_sorensen, 
                          colour = data, group = data, fill = data)) +
  geom_point() +
  scale_colour_manual(values = c("#a26929", "#91367e", "#7a8235", "#009392", "#3b738f")) +
  scale_fill_manual(values = c("#a26929", "#91367e", "#7a8235", "#009392", "#3b738f")) +
  scale_x_continuous(breaks = c(0.01, 0.25, 0.50, 0.75, 1)) +
  geom_smooth() +
    guides(fill = F, colour = F) +
  change_theme() +
  labs(x = "\nProportion of data included",
       y = "Mean overlap with\nglobal change space\n"))

ggsave(overlap_accumulation, filename ="figures/overlap_accumulation.png",
       height = 6, width = 8)


full_gc_panel <- grid.arrange(combined_pop_pca_terr, combined_bio_pca_terr, combined_predicts_pca_terr,
                              combined_pop_pca_mar, combined_bio_pca_mar, overlap_accumulation,
                              ncol = 3)

ggsave(full_gc_panel, filename = "figures/full_gc_panel2022.pdf",
       height = 10, width = 16)

ggsave(full_gc_panel, filename = "figures/full_gc_panel2022.png",
       height = 10, width = 16)

# Fig 2 Distributions of driver intensities ----
# Standardise between 0 and 1
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
lpd.coords <- lpd.coords %>% gather(driver, intensity, 11:16)

bt.coords$database <- "BioTIME"
bt.coords <- bt.coords %>% gather(driver, intensity, 11:16)

colnames(predicts_drivers)
colnames(lpd.coords)
colnames(bt.coords)
colnames(bt.coords)[2] <- "study_id"
colnames(lpd.coords)[2] <- "study_id"

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

(driver_density_plot_mar <- ggplot(drivers_combined_r[drivers_combined_r$realm == "Marine",], 
                                   aes(x = intensity, y = fct_rev(driver), 
                                       colour = fct_rev(sampling2),
                                       fill = fct_rev(sampling2))) +
    geom_density_ridges(alpha = 0, scale = 0.95, size = 1) +
    #change_theme() +
    theme_void() +
    scale_fill_manual(values = c("grey80","#009392", "#3b738f")) +
    scale_colour_manual(values = c("grey80","#009392", "#3b738f")) +
    facet_grid(~fct_rev(database)) +
    labs(y = NULL, x = "\nDriver intensity") +
    scale_x_continuous(limits = c(0, 1),
                       breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    theme(strip.background = element_blank(),
          axis.text.y = element_text(size = 14, hjust = 0.95, vjust = 0.2),
          axis.text.x = element_text(size = 14, hjust = 0.95, vjust = 0.2),
          axis.title.x = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          plot.title = element_text(size = 14, hjust = 0.05, 
                                    color = "grey20", face = "bold")) +
    guides(fill = F, colour = F, linetype = F))

drivers_combined_r_pr_r <- drivers_combined_r %>% filter(sampling2 == "TerrestrialRandom sampling") 
drivers_combined_r_pr_r$database <- "PREDICTS"

drivers_combined_r <- bind_rows(drivers_combined_r, drivers_combined_r_pr_r)

# reorder levels
summary(as.factor(drivers_combined_r$database))
drivers_combined_r$database <- factor(drivers_combined_r$database,
                                      levels = c("LPD", "BioTIME", "PREDICTS"),
                                      labels = c("LPD", "BioTIME", "PREDICTS"))

(driver_density_plot_terr <- ggplot(drivers_combined_r[drivers_combined_r$realm == "Terrestrial",], 
                                    aes(x = intensity, y = fct_rev(driver), 
                                        colour = fct_rev(sampling2),
                                        fill = fct_rev(sampling2))) +
    geom_density_ridges(alpha = 0, scale = 0.95, size = 1) +
    theme_void() +
    scale_fill_manual(values = c("grey80","#a26929", "#91367e", "#7a8235")) +
    scale_colour_manual(values = c("grey80","#a26929", "#91367e", "#7a8235")) +
    facet_grid(~database) +
    labs(y = NULL, x = "\nDriver intensity") +
    scale_x_continuous(limits = c(0, 1),
                       breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    theme(strip.background = element_blank(),
          axis.text.y = element_text(size = 14, hjust = 0.95, vjust = 0.2),
          axis.text.x = element_text(size = 14, hjust = 0.95, vjust = 0.2),
          axis.title.x = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          plot.title = element_text(size = 14, hjust = 0.05, 
                                    color = "grey20", face = "bold")) +
    guides(fill = F, colour = F, linetype = F))

# Saving the graphs
ggsave(driver_density_plot_terr, filename = "figures/distributions2022_terr.pdf", height = 7, width = 8.3)
ggsave(driver_density_plot_terr, filename = "figures/distributions2022_terr.png", height = 7, width = 8.3)

ggsave(driver_density_plot_mar, filename = "figures/distributions2022_mar.pdf", height = 7, width = 6)
ggsave(driver_density_plot_mar, filename = "figures/distributions2022_mar.png", height = 7, width = 6)

# Load model predictions
# Models were run in 05-run-models
# Predictions were calculated in 06-calculate-predictions
load("data/output/predictions_combo2022.RData")

(predictions_deviation_plot <- ggplot(predictions_combo, aes(x = deviation, y = x)) +
    geom_bar(aes(fill = data_realm), stat = "identity", 
             position = position_dodge2(width = 0.9, preserve = "single"),
             alpha = 0.8, width = 7) +
    #geom_errorbar(aes(xmin = deviation_lower, xmax = deviation_upper),
    #             width=.2, position = position_dodge2(width = 1.2, preserve = "single")) +
    facet_grid(~fct_rev(realm), scales = "free", space='free') +
    labs(y = NULL, x = "Predicted deviation\nfrom global mean\n") +
    geom_vline (xintercept = 0, linetype = "dotted") +
    scale_fill_manual(values = c("#a26929", "#91367e", "#7a8235", "#009392", "#3b738f")) +
    scale_colour_manual(values = c("#a26929", "#91367e", "#7a8235", "#009392", "#3b738f")) +
    theme_classic() +
    coord_flip() +
    guides(fill = F, colour = F) +
    theme(legend.position = c(0.14, 0.12),
          legend.title = element_text(color = "black", size = 10),
          text = element_text(color = "#22211d"),
          axis.text.y = element_text(size = 14, hjust = 0.95, vjust = 0.2),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.line.x = element_line(size = 0.3),
          axis.line.y = element_line(size = 0.3),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(2, 2, 2, 2),"cm"),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 14, hjust = 0.05, 
                                    color = "grey20", face = "bold")))

ggsave(predictions_deviation_plot, filename = "figures/predictions_deviation_plot2022.pdf", height = 4, width = 14.3)
ggsave(predictions_deviation_plot, filename = "figures/predictions_deviation_plot2022.png", height = 4, width = 14.3)

# Fig 3 Global change temporal representation ----
load("data/output/CRU_BioTIME_mean2022.RData")
load("data/output/CRU_LPD_mean2022.RData")
load("data/output/NOAA_BioTIME_mean2022.RData")
load("data/output/NOAA_LPD_mean2022.RData")

# Filter out the temperature records for periods
# equal to the duration of each time-series (before they started)

# e.g. if a time-series is from 1980 to 1990, the matching before period
# is 1970 to 1980

bt <- popbio %>% filter(type == "Biodiversity")

bt$start_comparison <- bt$start_year - bt$duration
bt$end_comparison <- bt$start_year - 1

bt_simpler <- bt %>% dplyr::select(timeseries_id, start_year, end_year, start_comparison,
                                   end_comparison, duration)
colnames(bt_simpler)[1] <- "rarefyID"

temp_bt <- left_join(tmp_sites_long, bt_simpler, by = "rarefyID")

str(temp_bt)

temp_bt2 <- temp_bt %>% dplyr::filter(year > start_comparison -1)
temp_bt3 <- temp_bt2 %>% dplyr::filter(year < end_year + 1)

temp_bt3 <- temp_bt3 %>% mutate(timing = case_when(year < end_comparison + 1 ~ "Before monitoring",
                                                   year > start_year - 1 ~ "During monitoring"))

temp_bt3b <- na.omit(temp_bt3)

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
    ggtitle("\nd     Climate warming over time (BioTIME terrestrial)\n") +
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
# logical     659    1254       2 

1254/(1254+659)

# 65.5% of BioTIME terr time series experienced more warming during time series than before

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

# # Keep only time series with duration of more than 4 years
 temp_lpd3b <- temp_lpd3b %>%
   filter(duration > 3)

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
    ggtitle("\nc     Climate warming over time (Living Planet)\n") +
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
# logical    1390    2384 

2384/(2384+1390)

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
    ggtitle("\ne     Climate warming over time (Living Planet marine)\n") +
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
    scale_fill_gradient(low = "#6a97ba", high = "#293c4b") +
    ggtitle("\nb     Climate warming over time (BioTIME marine)\n") +
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
# logical    4252    3022 

3022/(3022+4152)

# 0.4212434

# 46.7% of BioTIME mar time series experienced more warming during time series than before

warming_terr_mar_panel <- grid.arrange(lpd_climate1, bt_climate1, lpd_climate2, bt_climate2, ncol = 2)

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

# Rate of forest cover change
load("data/input/luh_polys_longUpdated.RData")
forest_bt <- left_join(luh_polys_longUPdated, bt_simpler, by = "rarefyID")

str(forest_bt)

forest_bt2 <- forest_bt %>% dplyr::filter(year > start_comparison -1)
forest_bt3 <- forest_bt2 %>% dplyr::filter(year < end_year + 1)

forest_bt3 <- forest_bt3 %>% mutate(timing = case_when(year < end_comparison + 1 ~ "Before monitoring",
                                                       year > start_year - 1 ~ "During monitoring"))

forest_bt3_filt <- forest_bt3 %>% filter(primf > 0.5)

amount1 <- forest_bt3 %>% filter(year == start_comparison |
                                   year == end_comparison)

amount1 <- amount1 %>% group_by(rarefyID) %>% 
  mutate(forest_loss_before = lag(primf, default = first(primf)) - primf) %>%
  filter(row_number() == 2L)

amount2 <- forest_bt3 %>% filter(year == start_year |
                                   year == end_year)

amount2 <- amount2 %>% group_by(rarefyID) %>% 
  mutate(forest_loss_during = lag(primf, default = first(primf)) - primf) %>%
  filter(row_number() == 2L)

amount2_simple <- amount2 %>% dplyr::select(rarefyID, forest_loss_during)

amount1 <- left_join(amount1, amount2_simple, by = "rarefyID")

# Line plots

str(forest_bt)

forest_bt$year <- as.numeric(as.character(forest_bt$year))

bt_forest_filt <- forest_bt %>% filter(primf != 0)

bt_forest_filt1 <- bt_forest_filt %>% filter(year < start_year)
bt_forest_filt2 <- bt_forest_filt %>% filter(year > start_year - 1)

(bt_forest_time <- ggplot() +
    geom_line(data = bt_forest_filt1, aes(x = year, y = primf*100, group = rarefyID),
              alpha = 0.03, colour = "grey30") +
    geom_line(data = bt_forest_filt2, aes(x = year, y = primf*100, group = rarefyID),
              alpha = 0.8, colour = "#91367e") +
    labs(x = "\nYear\n", y = "Primary forest cover (%)") +
    ggtitle("\nb     Forest loss over time (BioTIME)") +
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

ggsave(bt_forest_time, filename = "figures/bt_forest_time_2021.pdf", height = 4, width = 10)
ggsave(bt_forest_time, filename = "figures/bt_forest_time_2021.png", height = 4, width = 10)

ggplot(forest_bt3_filt, aes(x = primf)) +
  geom_density(aes(fill = timing), alpha = 0.6)

load("data/input/LUH_polys_long_lpi.RData")
luh_polys_lpi_long <- ungroup(luh_polys_lpi_long)
names(lpd_simpler)[1] <- "id"
lpd_simpler$id <- as.character(lpd_simpler$id)
luh_polys_lpi_long$id <- as.character(luh_polys_lpi_long$id)

forest_lpd <- left_join(luh_polys_lpi_long, lpd_simpler, by = "id")

str(forest_lpd)

forest_lpd2 <- forest_lpd %>% dplyr::filter(year > start_comparison -1)
forest_lpd3 <- forest_lpd2 %>% dplyr::filter(year < end_year + 1)

forest_lpd3 <- forest_lpd3 %>% mutate(timing = case_when(year < end_comparison + 1 ~ "Before monitoring",
                                                         year > start_year - 1 ~ "During monitoring"))

forest_lpd3_filt <- forest_lpd3 %>% filter(primf > 0.5)

amount1 <- forest_lpd3 %>% filter(year == start_comparison |
                                    year == end_comparison)

amount1 <- amount1 %>% group_by(id) %>% 
  mutate(forest_loss_before = lag(primf, default = first(primf)) - primf) %>%
  filter(row_number() == 2L)

amount2 <- forest_lpd3 %>% filter(year == start_year |
                                    year == end_year)

amount2 <- amount2 %>% group_by(id) %>% 
  mutate(forest_loss_during = lag(primf, default = first(primf)) - primf) %>%
  filter(row_number() == 2L)

amount2_simple <- amount2 %>% dplyr::select(id, forest_loss_during)

amount1 <- left_join(amount1, amount2_simple, by = "id")

# Line plots

str(forest_lpd)

forest_lpd$year <- as.numeric(as.character(forest_lpd$year))

lpd_forest_filt <- forest_lpd %>% filter(primf != 0)

lpd_forest_filt1 <- lpd_forest_filt %>% filter(year < start_year)
lpd_forest_filt2 <- lpd_forest_filt %>% filter(year > start_year - 1)

(lpd_forest_time <- ggplot() +
    geom_line(data = lpd_forest_filt1, aes(x = year, y = primf*100, group = id),
              alpha = 0.03, colour = "grey30") +
    #geom_line(data = bt_forest_filt2, aes(x = year, y = primf*100, group = rarefyID),
    #          alpha = 0.8, colour = "#91367e") +
    geom_line(data = lpd_forest_filt2, aes(x = year, y = primf*100, group = id),
              alpha = 0.8, colour = "#a26929") +
    labs(x = "\nYear\n", y = "Primary forest cover (%)") +
    ggtitle("\na     Forest loss over time (Living Planet)") +
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

ggsave(lpd_forest_time, filename = "figures/lpd_forest_time_2021.pdf", height = 4, width = 10)
ggsave(lpd_forest_time, filename = "figures/lpd_forest_time_2021.png", height = 4, width = 10)

forest_panel <- grid.arrange(lpd_forest_time, bt_forest_time, ncol = 2)
ggsave(forest_panel, filename = "figures/lpd_bt_forest_panel_2021.pdf", height = 4, width = 12)

terr_panel <- grid.arrange(forest_panel, warming_terr_panel, nrow = 1, widths = c(0.65, 0.35))

ggsave(terr_panel, filename = "figures/lpd_bt_terr_forest_time_2021.pdf", height = 10, width = 16)
ggsave(terr_panel, filename = "figures/lpd_bt_terr_forest_time_2021.png", height = 10, width = 16)
