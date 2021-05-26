# Figures for pre-registration
# Gergana Daskalova
# Updated 1st March 2019

# devtools::install_github("wilkox/treemapify")
library(treemapify)
library(viridis)
library(ggalt)
library(ggthemes)
library(tidyverse)
library(stargazer)
library(gridExtra)
library(dggridR)
library(raster)
library(rgdal)
library(ggfortify)

# Load population, biodiversity and driver data
load("data/output/popbio.RData")
cumulative <- stack("data/input/Cumulative_4_GD.grd")


driver_theme <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 16), 
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
          legend.position = c(0.9, 0.9), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}


# Taxonomic representation ----
pop_sum <- popbio %>%
  filter(type == "Population") %>%
  filter(taxa %in% c("Aves", "Amphibia", "Mammalia", "Reptilia",
                      "Actinopterygii", "Elasmobranchii")) %>%
  group_by(taxa) %>% tally()

pop_sum$taxa <- factor(pop_sum$taxa,
                        levels = c("Elasmobranchii", "Actinopterygii",
                                   "Amphibia", "Aves", "Mammalia", "Reptilia"),
                        labels = c("Elasmobranchii", "Actinopterygii",
                                   "Amphibia", "Aves", "Mammalia", "Reptilia"))

(pop_area <- ggplot(pop_sum, aes(area = n, fill = taxa, label = n,
                    subgroup = taxa)) +
         geom_treemap() +
    geom_treemap_subgroup_border(colour = "white", size = 1) +
    geom_treemap_text(colour = "white", place = "center", reflow = T) +
    scale_fill_manual(values = c("#59BAC0", "#3f517a", "#25662c", "#d8b70a", "#972d15", "#dc863b")) +
    guides(fill = FALSE))

bio_sum <- popbio %>%
  filter(type == "Biodiversity") %>%
  distinct() %>%
  group_by(taxa) %>% tally()

bio_sum$taxa <- factor(bio_sum$taxa,
                        levels = c("Fish",
                                   "Amphibians",
                                   "Birds",
                                   "Mammals",
                                   "Terrestrial invertebrates",
                                   "Terrestrial plants",
                                   "Marine invertebrates",      
                                   "Marine invertebrates/plants",
                                   "Benthos",
                                   "All"),
                        labels = c("Fish",
                                   "Amphibians",
                                   "Birds",
                                   "Mammals",
                                   "Terrestrial invertebrates",
                                   "Terrestrial plants",
                                   "Marine invertebrates",      
                                   "Marine invertebrates/plants",
                                   "Benthos",
                                   "Multiple taxa"))

(bio_area <- ggplot(bio_sum, aes(area = n, fill = taxa, label = n,
                    subgroup = taxa)) +
    geom_treemap() +
    geom_treemap_subgroup_border(colour = "white", size = 1) +
    geom_treemap_text(colour = "white", place = "center", reflow = T) + 
    scale_fill_manual(values = c("#3f517a", "#25662c", "#d8b70a", "#972d15",
                               "darkcyan", "#919c4c", "deepskyblue4", "thistle3",
                               "violetred4", "#828585")) +
    guides(fill = FALSE))

samples_panel <-  grid.arrange(pop_area, bio_area, ncol = 2)
ggsave(samples_panel, filename = "figures/samples.pdf", dpi = 300, device = "pdf",
       height = 3, width = 10)

# Geographic representation ----
world <- map_data("world")

mus.coords <- popbio %>%
  filter(type == "Population")

mus.coords$taxa <- factor(mus.coords$taxa,
                        levels = c("Elasmobranchii", "Actinopterygii",
                                   "Amphibia", "Aves", "Mammalia", "Reptilia"),
                        labels = c("Elasmobranchii", "Actinopterygii",
                                   "Amphibia", "Aves", "Mammalia", "Reptilia"))

(pop.map <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    coord_proj("+proj=eck4") +
    theme_map() +
    geom_point(data = mus.coords[mus.coords$taxa == "Aves",], 
               aes(x = long, y = lat),
               alpha = 0.6, colour = "#25662c", size = 0.5) +
    geom_point(data = mus.coords[mus.coords$taxa != "Aves",], 
               aes(x = long, y = lat, colour = taxa),
               alpha = 0.6, size = 0.5) +
    scale_colour_manual(values = c("#59BAC0", "#3f517a", "#d8b70a", "#972d15", "#dc863b")) +
    scale_y_continuous(limits = c(-80, 80)) +
    guides(colour = FALSE) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 10),
          legend.justification = "top"))

bt.coords <- popbio %>%
  filter(type == "Biodiversity")
  
bt.coords$taxa <- factor(bt.coords$taxa,
                         levels = c("Fish",
                                    "Amphibians",
                                    "Birds",
                                    "Mammals",
                                    "Terrestrial invertebrates",
                                    "Terrestrial plants",
                                    "Marine invertebrates",      
                                    "Marine invertebrates/plants",
                                    "Benthos",
                                    "All"),
                         labels = c("Fish",
                                    "Amphibians",
                                    "Birds",
                                    "Mammals",
                                    "Terrestrial invertebrates",
                                    "Terrestrial plants",
                                    "Marine invertebrates",      
                                    "Marine invertebrates/plants",
                                    "Benthos",
                                    "Multiple taxa"))

(bio.map <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    coord_proj("+proj=eck4") +
    theme_map() +
    geom_point(data = bt.coords[bt.coords$taxa == "Birds",],
               aes(x = long, y = lat), colour = "#d8b70a",
               alpha = 0.6, size = 0.5) +
    geom_point(data = bt.coords[bt.coords$taxa != "Birds",],
               aes(x = long, y = lat, colour = taxa),
               alpha = 0.6, size = 0.5) +
    scale_colour_manual(values = c("#3f517a", "#25662c", "#972d15",
                                   "darkcyan", "#919c4c", "deepskyblue4", "thistle3",
                                   "violetred4", "#828585")) +
    scale_y_continuous(limits = c(-80, 80)) +
    guides(colour = FALSE) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 10),
          legend.justification = "top"))

map.panel <- grid.arrange(pop.map, bio.map, ncol = 2)
ggsave(map.panel, filename = "figures/map_panel.pdf", dpi = 300, device = "pdf",
       height = 5, width = 10)

# Temporal representation ----
# Create a sorting variable
mus.coords$sort <- mus.coords$taxa
mus.coords$sort <- factor(mus.coords$sort, levels = c("Elasmobranchii",
                                                      "Actinopterygii",
                                                      "Amphibia",
                                                      "Aves",
                                                      "Mammalia",
                                                      "Reptilia"),
                          labels = c(1, 2, 3, 4, 5, 6))

mus.coords$sort <- paste0(mus.coords$sort, mus.coords$start_year)
mus.coords$sort <- as.numeric(as.character(mus.coords$sort))

(timeline_lpd <- ggplot() +
    geom_linerange(data = mus.coords, aes(ymin = start_year, ymax = end_year, 
                                         colour = taxa,
                                         x = fct_reorder(timeseries_id, desc(sort))),
                  size = 0.3) +
    scale_colour_manual(values = c("#59BAC0", "#3f517a", "#25662c", "#d8b70a", "#972d15", "#dc863b")) +
    labs(x = "\nYear", y = NULL) +
    theme_bw() +
    coord_flip() +
    guides(colour = F) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(),
          axis.ticks = element_blank(),
          legend.position = "bottom", 
          panel.border = element_blank(),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0),
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 20)))

bt.coords$sort <- bt.coords$taxa
bt.coords$sort <- factor(bt.coords$sort, levels = c("Fish",
                                                      "Amphibians",
                                                      "Birds",
                                                      "Mammals",
                                                      "Terrestrial invertebrates",
                                                      "Terrestrial plants",
                                                      "Marine invertebrates",
                                                      "Marine invertebrates/plants",
                                                      "Benthos",
                                                      "Multiple taxa"),
                          labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

bt.coords$sort <- paste0(bt.coords$sort, bt.coords$start_year)
bt.coords$sort <- as.numeric(as.character(bt.coords$sort))

(timeline_bt <- ggplot() +
    geom_linerange(data = bt.coords, aes(ymin = start_year, ymax = end_year, 
                                          colour = taxa,
                                          x = fct_reorder(timeseries_id, desc(sort))),
                   size = 0.3) +
    scale_colour_manual(values = c("#3f517a", "#25662c", "#d8b70a", "#972d15",
                                   "darkcyan", "#919c4c", "deepskyblue4", "thistle3",
                                   "violetred4", "#828585")) +
    scale_y_continuous(breaks = seq(1850, 2020, 30)) +
    labs(x = "\nYear", y = NULL) +
    theme_bw() +
    coord_flip() +
    guides(colour = F) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(),
          axis.ticks = element_blank(),
          legend.position = "bottom", 
          panel.border = element_blank(),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0),
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 20)))

timeline_panel <-  grid.arrange(timeline_lpd, timeline_bt, ncol = 2)
ggsave(timeline_panel, filename = "figures/timelines.pdf", dpi = 300, device = "pdf",
       height = 3, width = 10)

# Global change representation ----

popbio_long <- popbio %>% gather(driver, intensity, select = 15:20)

popbio_long$driver <- factor(popbio_long$driver,
                             levels = c("climate_change", "human_use", "human_population",
                                        "pollution", "invasions", "cumulative"),
                             labels =c("climate_change", "human_use", "human_population",
                                       "pollution", "invasions", "cumulative"))

popbio_long$realm <- factor(popbio_long$realm,
                            levels = c("Terrestrial", "Marine"),
                            labels = c("Terrestrial", "Marine"))

popbio_long$type <- factor(popbio_long$type,
                            levels = c("Population", "Biodiversity"),
                            labels = c("Population", "Biodiversity"))

(hist1 <- ggplot(popbio_long, aes(x = intensity, colour = type, fill = type)) +
    geom_histogram(alpha = 0.6) +
    facet_grid(type + realm ~ driver, scales = "free") +
    driver_theme() +
    scale_colour_manual(values = c("#CC9145", "#6eafae")) +
    scale_fill_manual(values = c("#CC9145", "#6eafae")) +
    theme(strip.background = element_blank(), 
          strip.text = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          panel.spacing = unit(1.5, "lines")) +
    guides(colour = F, fill = F) +
    #scale_x_continuous(limits = c(0, 13),
    #                   breaks = c(1, 4, 7, 10, 13)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "\nIntensity", y = "Number of time-series\n"))

ggsave(hist1, filename = "figures/histograms.pdf", device = "pdf",
       dpi = 300, height = 10, width = 15)

# Driver combination representation ----
names(cumulative)
cc <- subset(cumulative, 1)
hu <- subset(cumulative, 2)
hp <- subset(cumulative, 3)
po <- subset(cumulative, 4)
inv <- subset(cumulative, 5)
all <- subset(cumulative, 6)

# PCAs
N <- 7000  # How many cells to sample

dggs <- dgconstruct(res = 12, metric = FALSE, resround = 'nearest')

maxcell <- dgmaxcell(dggs)                     #Get maximum cell id
cells <- sample(1:maxcell, N, replace = FALSE) #Choose random cells
grid <- dgcellstogrid(dggs, cells, frame = TRUE, wrapcells = TRUE) #Get grid

#Get polygons for each country of the world
countries <- map_data("world")

#Plot everything on a flat map
(p <- ggplot() + 
    geom_polygon(data = countries, aes(x = long, y = lat, group = group), 
                 fill = "grey60", color = "grey60") +
    coord_proj("+proj=eck4") +
    theme_map() +
    geom_point(data = grid, aes(x = long, y = lat, group = group), 
               colour = "deepskyblue4", alpha = 0.05))

ggsave(p, filename = "figures/map_ideal_sampling.pdf", dpi = 300, device = "pdf",
       height = 5, width = 5)

# The world is covered!

world_coords <- grid %>% dplyr::select(cell, long, lat) %>% distinct() %>%
  group_by(cell) %>% slice(1)

colnames(world_coords) <- c("cell", "long", "lat")

coords_again <- SpatialPoints(world_coords[, c("long", "lat")])
world_spatial_df <- SpatialPointsDataFrame(coords_again, world_coords)

random_coords <- SpatialPoints(cbind(world_spatial_df$long, world_spatial_df$lat),
                          proj4string = CRS("+proj=longlat"))

random_coords <- spTransform(random_coords, CRS("+proj=eck4 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
random <- as.data.frame(random_coords@coords)
colnames(random) <- c("long", "lat")

df7 <- as.data.frame(raster::extract(cc, random))
colnames(df7) <- "climate_change"
df8 <- as.data.frame(raster::extract(hu, random))
colnames(df8) <- "human_use"
df9 <- as.data.frame(raster::extract(hp, random))
colnames(df9) <- "human_population"
df10 <- as.data.frame(raster::extract(po, random))
colnames(df10) <- "pollution"
df11 <- as.data.frame(raster::extract(inv, random))
colnames(df11) <- "invasions"
df12 <- as.data.frame(raster::extract(all, random))
colnames(df12) <- "cumulative"

random_drivers <- cbind(df7, df8, df9, df10, df11, df12)
random_drivers <- na.omit(random_drivers)

(random_pca <- autoplot(prcomp(random_drivers), loadings = TRUE,
                        loadings.label = TRUE, loadings.label.size = 3,
                        colour = "grey80") +
    driver_theme())

mus.coords$sampling <- "Population time-series"
random_drivers$sampling <- "Random sampling"

combined_pop <- rbind(random_drivers, mus.coords[,c(15:20, 23)])

combined_pop$sampling <- factor(combined_pop$sampling, 
                                levels = c("Random sampling", "Population time-series"),
                                labels = c("Random sampling", "Population time-series"))

combined_pop <- na.omit(combined_pop)

(combined_pop_pca <- autoplot(prcomp(combined_pop[1:6]), data = combined_pop, colour = 'sampling',
                              loadings = TRUE, loadings.colour = 'black', scale = 1,
                              loadings.label.colour = 'black', alpha = 0.3,
                              #          ) +
                              #,
                              loadings.label = TRUE, loadings.label.size = 5,
                              loadings.label.repel = TRUE, loadings.label.vjust = 0.2) +
    scale_colour_manual(values = c("grey80", "#CC9145")) +
    scale_x_continuous(breaks = c(-0.01, 0, 0.01),
                       labels = c("-0.01", "0", "0.01")) +
    scale_y_continuous(breaks = c(-0.01, 0, 0.01),
                       labels = c("-0.01", "0", "0.01")) +
    guides(colour = FALSE) +
    theme_LPI())

bt.coords$sampling <- "Biodiversity time-series"

combined_bio <- rbind(random_drivers, bt.coords[,c(15:20, 23)])

combined_bio$sampling <- factor(combined_bio$sampling, 
                            levels = c("Random sampling", "Biodiversity time-series"),
                            labels = c("Random sampling", "Biodiversity time-series"))

combined_bio <- na.omit(combined_bio)

(combined_bio_pca <- autoplot(prcomp(combined_bio[1:6]), data = combined_bio, colour = 'sampling',
                          loadings = TRUE, loadings.colour = 'black', scale = 1,
                          loadings.label.colour = 'black', alpha = 0.3,
                          loadings.label = TRUE, loadings.label.size = 5,
                          loadings.label.repel = TRUE, loadings.label.vjust = 0.2) +
    scale_colour_manual(values = c("grey80", "#6eafae")) +
    scale_x_continuous(breaks = c(-0.01, 0, 0.01),
                       labels = c("-0.01", "0", "0.01")) +
    scale_y_continuous(breaks = c(-0.01, 0, 0.01),
                       labels = c("-0.01", "0", "0.01")) +
    guides(colour = FALSE) +
    driver_theme())

pca_panel <- grid.arrange(combined_pop_pca, combined_bio_pca, ncol = 2)
ggsave(pca_panel, filename = "figures/lpd_pca_no_labels.pdf", height = 7, width = 14,
       dpi = 300, device = "pdf")

(combined_bio_pca_five <- autoplot(prcomp(combined_bio[1:5]), data = combined_bio, colour = 'sampling',
                              loadings = TRUE, loadings.colour = 'black', scale = 1,
                              loadings.label.colour = 'black', alpha = 0.3) +
                              #,
                              #loadings.label = TRUE, loadings.label.size = 5,
                             # loadings.label.repel = TRUE, loadings.label.vjust = 0.2) +
    scale_colour_manual(values = c("grey80", "#6eafae")) +
    scale_x_continuous(breaks = c(-0.01, 0, 0.01),
                       labels = c("-0.01", "0", "0.01")) +
    scale_y_continuous(breaks = c(-0.01, 0, 0.01),
                       labels = c("-0.01", "0", "0.01")) +
    guides(colour = FALSE) +
    driver_theme())

ggsave(combined_bio_pca_five, filename = "figures/bio_pca_drivers_no_labels.pdf", height = 7, width = 7,
       dpi = 300, device = "pdf")

ggsave(combined_bio_pca, filename = "figures/bio_pca_drivers_all.pdf", height = 7, width = 7,
       dpi = 300, device = "pdf")

ggsave(combined_bio_pca, filename = "figures/bio_pca_drivers_all_no_labels.pdf", height = 7, width = 7,
       dpi = 300, device = "pdf")


# Sample size table
popbio$taxa <- as.factor(popbio$taxa)
sample_size <- popbio %>% group_by(type, realm, taxa) %>% tally()
write.csv(sample_size, file = "data/output/sample_size.csv")

data_samples <- rbind(head(popbio), popbio[5000:5006,])
write.csv(data_samples, file = "data/output/data_samples.csv")

# Maps ----

library(rasterVis)

levelplot(cumulative)
levelplot(cumulative, layers = c(3,5))
bwplot(cumulative)

library(sdmvspecies)

the_five <- subset(cumulative, 1:5)
test <- sdmvspecies::rescale(the_five)
levelplot(test)

test2 <- subset(cumulative, 6)
values(test2) <- scales::rescale(values(test2), to = c(0, 4))

new_cumulative <- stack(test, test2)

pdf("name.pdf")
print(levelplot(new_cumulative))
dev.off()

maps <- levelplot(new_cumulative)
maps

popbio <- arrange(popbio, type)

(map <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    coord_proj("+proj=eck4") +
    theme_map() +
    #geom_point(data = popbio, 
               aes(x = long, y = lat, colour = type),
               alpha = 0.6, size = 0.5) +
    scale_colour_manual(values = c("#6eafae", "#CC9145")) +
    scale_y_continuous(limits = c(-80, 80)) +
    guides(colour = FALSE) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 10),
          legend.justification = "top"))

ggsave(map, filename = "figures/map.pdf", dpi = 300, device = "pdf",
       height = 5, width = 5)
