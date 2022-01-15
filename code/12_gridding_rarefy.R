################################################################
### BodySize
### October 2021
### 
###  Full dataset: Gridding & Rerefy
###
### Ines S. Martins (iism1@st-andrews.ac.uk) 
################################################################

### BASED ON:
# https://zenodo.org/record/3367444: sChange-workshop-BioGeo-BioDiv-Change-51fa6f7
# from Shane Blowes and Sarah Supp - 23 Feb 2017
# 01_Study_to_Grid
# 02_rarefy_griddedData_clusterVersion.R


# clear workspace
rm(list=ls())

# set working directory
setwd("/iism1/Dropbox (iDiv)/BodySize/0_data/")

# Load necessary packages
library(readr)
library("tidyr")
library("tidyverse")
library("rgeos")
library(dplyr)
library(purrr)
library(dggridR)
library(geosphere)
library(iNEXT)
library(lazyeval)
library("vegan")
library("scales")
library("rnaturalearth")
library("rnaturalearthdata")


#### 0 - Import data =========================================================
bt_meta <-bioTIMEMetadataJune2021 <- read_csv("bioTIMEMetadataJune2021.csv")
bt_abu_bio<-readRDS("BioTIMEQ2021.rds")

### for testing
bt_meta<-subset(bt_meta, TAXA=="Reptiles")

## merge information
bt<-inner_join(bt_abu_bio,bt_meta,by="STUDY_ID")


#### 1. Gridding and Grouping observations =========================================================
##### 1.0. Distinguishing SL and ML studies =========================================================

# look at frequency distribution of SL studies, and check for any extreme outliers
bt_meta <- bt_meta %>%
  mutate(StudyMethod = ifelse(NUMBER_LAT_LONG == 1, "SL", NA))

# Calculate the extent and mean for SL studies, without the outlier[s]
SL_extent<- bt %>% 
  filter(NUMBER_LAT_LONG == 1,
         AREA_SQ_KM<=500) %>% 
  summarise(extent_mean = mean(AREA_SQ_KM),    
            extent_sd = sd(AREA_SQ_KM))

SL_extent_sd   <- SL_extent$extent_sd
SL_extent_mean <- SL_extent$extent_mean

# Change the MLs into SLs that satisfy the criterion (< mean + sd)
bt <- bt %>%
  mutate(StudyMethod = ifelse(NUMBER_LAT_LONG == 1, "SL", "ML"),
         StudyMethod = ifelse(AREA_SQ_KM < (SL_extent_mean+SL_extent_sd), "SL", StudyMethod))

# We'll want to 'grid' SL and ML studies differently, add new coords to dataframe
# If StudyMethod=='SL', want to use the central lat and long, else
# ML uses observed coords for each observation
bt<- bt %>%
  mutate(lon_to_grid = if_else(StudyMethod=="SL", CENT_LONG, LONGITUDE), #if TRUE=CENT_LONG, if FALSE=LONGITUDE
         lat_to_grid = if_else(StudyMethod=="SL", CENT_LAT, LATITUDE))

bt %>% group_by()%>% 
  summarise(n_distinct(STUDY_ID),n_distinct(SPECIES),n_distinct(YEAR),min(YEAR),max(YEAR),n())-> bt_stats

##### 1.1. Gridding Observations =========================================================
## Create Grid
# package "dggridR" is not available in R 4.0, use 
#install_github('r-barnes/dggridR', vignette=TRUE) #instead
#	create a global grid with cells approximately equal to extent +/- sd of the 'true' SL studies?

dgg <- dgconstruct(res=12)  #resolution '12'= average cell area (96 km$^2$)

# #	determine the resolution closest to our cutoff point for SL vs ML studies
# res <- dg_closest_res_to_area(dgg, SL_extent_mean+SL_extent_sd)
# #	set the resolution
# dgg <- dgsetres(dgg, res)
#dgg <- dgsetres(dgg, res = 12)
dggetres(dgg) %>% filter( res == 12)

# Grid Multi-location studies - get the corresponding grid cells for all observations
#The 'dgtransform(dggs,lat,lon)' function has been deprecated. Please use 'dgGEO_to_SEQNUM(dggs,lon,lat)' instead! NOTE THE ARGUMENT ORDER!
#bt <- bt %>% mutate(cell = dgtransform(dgg, lat=lat_to_grid, lon=lon_to_grid))
bt <- bt %>% mutate(cell = dgGEO_to_SEQNUM(dgg,lon_to_grid,lat_to_grid)$seqnum)

##	Checks - what just happened?
check <- bt %>% group_by(StudyMethod, STUDY_ID) %>% summarise(n_cell = n_distinct(cell),
                                                              n_years= n_distinct(YEAR),max_YEAR= max(YEAR), min_YEAR= min(YEAR))
#	do all SL studies have one grid cell?
if (sum(dplyr::filter(check, StudyMethod=='SL') %>% .$n_cell != 1)==0) { 
  print("all SL studies have 1 grid cell") 
} else { print("ERROR: some SL studies have > 1 grid cell") }

# ok, how many cells/year/study are there? (e.g. how spread out were samples in a given study and year?)
check2 <- bt %>% 
  group_by(StudyMethod, STUDY_ID, YEAR) %>% 
  summarise(n_cell = n_distinct(cell))

range(check2$n_cell)


##### 1.2. Calculate new centers  (only ML) =========================================================
#	in order to calculate new centres we need this new bt object
#	add rarefyID:
bt <- bt %>% unite(col=rarefyID, STUDY_ID, cell, sep="_", remove=FALSE)
#save(bt, dgg, file='biotime_cells_68s.Rdata')
#load('~/Desktop/current/BioTime/data/biotime_cells_res11.Rdata')

bt %>% group_by()%>% 
  summarise(n_distinct(STUDY_ID),n_distinct(GENUS_SPECIES),n_distinct(rarefyID),n_distinct(paste(rarefyID,GENUS_SPECIES)),n_distinct(YEAR),min(YEAR),max(YEAR),n())-> bt_stats

##	want to calculate a single coordinate (i.e., the centre) for each rarefyID to place rarefyID's
##	within ecoregions....
##	also calculate cell_extent = area in convex hull (polygon) of points within a rarefyID. We can use this
##	to calculate a new 'extent' for each study AreaSum = sum(cell_extent) Issue #7 github
rarefyID_coords_nest <- ungroup(bt) %>%  ##	we don't need to do anything with the SL studies
  filter(StudyMethod!='SL') %>%  ##	select columns
  select(STUDY_ID, rarefyID, LONGITUDE, LATITUDE) %>%  ##	retain only uniqe locations within rarefyIDs (may want to change this if we want to weight the calculation of the centre)
  distinct(rarefyID, LONGITUDE, LATITUDE, .keep_all=TRUE) %>% ##	there are also some rarefyID's with only one unique geographic
  group_by(rarefyID) %>%
  mutate(n_locations = n_distinct(LONGITUDE,LATITUDE)) %>%
  ungroup() %>%  ##	drop rarefyIDs with only one location
  filter(n_locations > 1) %>% ##	drop our location counter
  select(-n_locations) %>% ##	group & nest 
  group_by(STUDY_ID, rarefyID) %>%
  nest()

##	i can't get purrr::map and chull to play nice so a loop it is! This takes awhile....
cell_extent <- numeric()
centre_rarefyID_x <- numeric()
centre_rarefyID_y <- numeric()
vertices_check <- data.frame()

for(i in 1:nrow(rarefyID_coords_nest)){ ##	sanity check 1.5h
  print(paste('rarefyID', i, 'out of', length(unique(rarefyID_coords_nest$rarefyID)))) ##	put a convex hull around the coords  
  hull = chull(x=unlist(rarefyID_coords_nest$data[[i]][,'LONGITUDE']), y=unlist(rarefyID_coords_nest $data[[i]][,'LATITUDE']))	 ##	get the vertices of the convex hull
  vertices = rarefyID_coords_nest$data[[i]][hull,c('LONGITUDE', 'LATITUDE')] ##	put some metadata together for checking later
  info = cbind.data.frame(Realm=rep(rarefyID_coords_nest$STUDY_ID[i], times=nrow(vertices)), rarefyID=rep(rarefyID_coords_nest$rarefyID[i], times=nrow(vertices)), vertices)
  vertices_check = rbind.data.frame(vertices_check, info)	# this could be used to check all these convex hulls
  ##	calculate the extent and centres (NB cell_extent==0 if there are only two points)
  cell_extent[i] = geosphere::areaPolygon(data.frame(x=vertices$LONGITUDE, y=vertices$LATITUDE))	# km2
  centre_rarefyID_x[i] = geosphere::geomean(cbind(x=vertices$LONGITUDE, y=vertices$LATITUDE))[1]
  centre_rarefyID_y[i] = geosphere::geomean(cbind(x=vertices$LONGITUDE, y=vertices$LATITUDE))[2]
}

##	combine STUDY_ID, rarefyID and the new cell_extent, and geographic centres
rarefyID_cell_centre <- cbind.data.frame(rarefyID_coords_nest[,1:2], cell_extent, rarefyID_x=centre_rarefyID_x, rarefyID_y=centre_rarefyID_y)		
rarefyID_cell_centre <- as_tibble(rarefyID_cell_centre)				  	


##### 1.3. Combine the SL + ML studies with only 1 location per rarefyID =========================================================

SL_coords <- ungroup(bt) %>% ##	get the SL studies
  filter(StudyMethod=='SL') %>% ##	get the SL studies
  select(STUDY_ID, rarefyID, CENT_LONG, CENT_LAT) %>% ##	select columns (use the central coords)
  mutate(cell_extent = 0, ##	create cell_extent column, and rename coords for joining with rarefyID_cell_centre
         rarefyID_x = CENT_LONG,
         rarefyID_y = CENT_LAT) %>%
  select(-CENT_LONG, -CENT_LAT)

ML_coords <- ungroup(bt) %>%
  filter(StudyMethod!='SL') %>%	 
  select(STUDY_ID, rarefyID, LONGITUDE, LATITUDE) %>% ##	select columns
  distinct(rarefyID, LONGITUDE, LATITUDE, .keep_all=TRUE) %>%  ##	retain only uniqe locations within rarefyIDs (may want to change this if we want to weight the calculation of the centre)
  group_by(rarefyID) %>% ##	there are also some rarefyID's with only one unique geographic
  mutate(n_locations = n_distinct(LONGITUDE,LATITUDE)) %>%
  ungroup() %>%  
  filter(n_locations == 1) %>% ##	retain rarefyIDs with only one location
  mutate(cell_extent = 0,##	create cell_extent column, and rename coords for joining with rarefyID_cell_centre
         rarefyID_x = LONGITUDE,
         rarefyID_y = LATITUDE) %>%
  select(-LONGITUDE, -LATITUDE, -n_locations)

#	put them together
rarefyID_cell_centre <- bind_rows(rarefyID_cell_centre, SL_coords, ML_coords)
#	not sure why but I have multiple entries for each rarefyID!
rarefyID_cell_centre <- rarefyID_cell_centre %>% distinct(STUDY_ID, rarefyID, cell_extent, rarefyID_x, rarefyID_y)
#save(rarefyID_cell_centre, file='rarefyID_cell_centre_181021.Rdata')

#### 2.  Prep data for rarefying  =========================================================
#	reduce to data required for rarefying, and rename a couple of columns
names(bt)[15]<-"Species"
bt$sum.allrawdata.BIOMASS<-as.numeric(bt$sum.allrawdata.BIOMASS)

#	NB: we will rarefy to the smallest number of ObsEventID's within studies
bt_grid <- bt %>% 
  dplyr::select(CLIMATE, REALM, TAXA, StudyMethod, ABUNDANCE_TYPE, BIOMASS_TYPE, STUDY_ID, YEAR, PLOT,
                cell, Species, sum.allrawdata.ABUNDANCE, sum.allrawdata.BIOMASS)
names(bt_grid)[12:13] <- c('Abundance', 'Biomass')

# add column for rarefyID (STUDY_ID + cell = samples within a study that are spatially grouped)
bt_grid <- bt_grid %>% unite(col=rarefyID, STUDY_ID, cell, sep="_", remove=FALSE)

# add column for ObsEventID (STUDY_ID + PLOT + YEAR = discrete sampling events within a year)
bt_grid <- bt_grid %>% unite(col=ObsEventID, rarefyID, PLOT, YEAR, sep="_", remove=FALSE)

#	combine with bt_grid (and we are ready to rarefy)
# bt_grid2 <- inner_join(ungroup(bt_grid), select(SampEvent_per_cell, STUDY_ID, YEAR, cell, n_samps))

#save(dgg, bt_grid, file='BioTIME_grid_181021.Rdata')


#### 3. Filter data on Coverage ==========================================
#	first collate taxa within cells, to calculate coverage for count data (see below for incidence) 
##### 3.1. *COUNT data coverage data coverage  ========================================

bt_grid_collate <- bt_grid %>%
  group_by(CLIMATE, REALM, TAXA, StudyMethod, ABUNDANCE_TYPE, rarefyID, STUDY_ID, YEAR, cell, Species) %>%
  summarise(
    Abundance = sum(as.numeric(Abundance), na.rm=TRUE),
    Biomass = sum(as.numeric(Biomass), na.rm=TRUE)) %>% 
  ungroup()

abund_coverage <- bt_grid_collate %>%  # get only rows representing count data (abundance)
  filter(ABUNDANCE_TYPE!='Presence/Absence' & ABUNDANCE_TYPE!='<NA>') %>%  # remove zeroes and NAs(some studies are designated count, but record Biomass estimates)
  filter(Abundance > 0 & !is.na(Abundance)) %>%
  group_by(CLIMATE, REALM, TAXA, StudyMethod, ABUNDANCE_TYPE, rarefyID, STUDY_ID, YEAR, cell) %>%
  summarise(
    singletons = sum(Abundance==1),# how many singletons
    doubletons = sum(Abundance==2),# how many doubletons
    N = sum(Abundance),# how many individuals in total sample
    Chat = 1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2*doubletons)),    # eqn 4a in Chao & Jost 2012 Ecology (==eqn 12 in Chao et al 2014 Ecol Monogr)
    # with fix from Chao et al 2014 Subroutine for iNEXT code (appendix)
    # correction for communities with no doubletons (this prevents NaN results for either singletons==0 or doubletons==0)
    Chat_corrected = ifelse(doubletons>0, 
                            1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2*doubletons)), 
                            1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2))),
    # the iNEXT coverage calculation has some extras in it that I don't understand
    # e.g., it corrects using f0.hat, the number of unobserved species (Chao1 estimator)? Why? Exactly how? Ref?
    f1_iNEXT = DataInfo(Abundance)$f1,
    f2_iNEXT = DataInfo(Abundance)$f2,
    #n_iNEXT = DataInfo(Abundance)$n,
    coverage = DataInfo(Abundance)$SC) %>%
  ungroup()

# are my estimates of singletons and doubletons equal to the iNEXT estimates?
sum(abund_coverage$singletons!=abund_coverage$f1_iNEXT)		# yes, good
sum(abund_coverage$doubletons!=abund_coverage$f2_iNEXT)		# yes, good

mn=mean(abund_coverage$Chat_corrected, na.rm=TRUE)
sd=sd(abund_coverage$Chat_corrected, na.rm=TRUE)
#abundance-based coverage of each(annual) sample (meaan=0.9795134, sd=0.0.05725734)

##### 3.2. *PRESENCE/ABSENCE data coverage (for INCIDENCE DATA) ========================================
#PRESENCE/ABSENCE coverage (for INCIDENCE DATA)... by rarefyID (year as samples)
bt_grid_collate_incidence <- bt_grid %>%
  # incidence data and remove NAs
  filter(ABUNDANCE_TYPE=='Presence/Absence' & ABUNDANCE_TYPE!='<NA>') %>%
  filter(!is.na(Abundance)) %>%
  #took out ObsEventID and PLOT as groups b/c no distinct plots for P/A data (Each year has only 1 sample)
  group_by(CLIMATE, REALM, TAXA, StudyMethod, ABUNDANCE_TYPE, rarefyID, STUDY_ID, YEAR, cell, Species) %>%
  summarise(
    Abundance = sum(Abundance, na.rm=TRUE)) %>%
  ungroup() %>%
  dplyr::select(-CLIMATE, -REALM, -TAXA, -StudyMethod, -ABUNDANCE_TYPE) %>%
  # create incidence column
  mutate(incidence = ifelse(Abundance==0, 0, 1))

pa_coverage <- data.frame()
for (r in unique(bt_grid_collate_incidence$rarefyID)){
  dat <- bt_grid_collate_incidence[bt_grid_collate_incidence$rarefyID == r,]
  if(length(unique(dat$YEAR))<2){ next }
  dat_pa <- dat %>%
    dplyr::select(YEAR, Species, incidence) %>%
    complete(YEAR, Species, fill = list(incidence = 0)) %>%
    #each row is a unique cell and Species, each column is a year
    spread(YEAR, incidence)
  # how many singletons (present in only 1 year)    #FIXME: This is a *really* ugly way to count
  singletons <- length(which(rowSums(dat_pa[,2:ncol(dat_pa)])==1))
  # how many doubletons (present in only 2 years)
  doubletons <- length(which(rowSums(dat_pa[,2:ncol(dat_pa)])==2))
  # how many incidences in the whole matrix
  N <- length(which(dat_pa[,2:ncol(dat_pa)] == 1))
  # how many samples (years) in the matrix
  T <- ncol(dat_pa) - 1
  # Chat for incidence: eqn Chat_sample(T) from Table 2, Chao et al. 2014, similar to eqn 4a in Chao & Jost Ecology
  Chat_i <- 1 - (singletons/N) * (((T-1)*singletons)/((T-1)*singletons + 2*doubletons))
  # with fix from Chao et al 2014 Subroutine for iNEXT code (appendix)
  # correction for communities with no doubletons (this prevents NaN results for
  # either singletons==0 or doubletons==0)
  if (doubletons > 0) {
    Chat_corrected_i <- Chat_i  }
  else {
    Chat_corrected_i <- 1 - (singletons/N) * (((T-1)*singletons)/((T-1)*singletons + 2))  }
  calcs <- data.frame(rarefyID=r, singletons, doubletons, N, Chat_i, Chat_corrected_i)
  pa_coverage <- rbind(pa_coverage, calcs)
  rm(dat_pa)
}

#what is the mean and sd of the presence data?
mnp=mean(pa_coverage$Chat_corrected_i)
sdp=sd(pa_coverage$Chat_corrected_i)

##### 3.3. *BIOMASS data coverage data coverage (for INCIDENCE DATA) ========================================
# BIOMASS DATA coverage (as INCIDENCE)... by rarefyID (year as samples)
bt_grid_collate_biomass <- bt_grid %>%
  # incidence data and remove NAs
  filter(is.na(ABUNDANCE_TYPE)) %>%
  filter(!is.na(Biomass)) %>%
  #took out ObsEventID and PLOT as groups b/c no distinct plots for biomass data (Few years have > 1 sample)
  group_by(CLIMATE, REALM, TAXA, StudyMethod, BIOMASS_TYPE, rarefyID, STUDY_ID, YEAR, cell, Species) %>%
  summarise(
    Biomass = sum(Biomass, na.rm=TRUE)) %>%
  ungroup() %>%
  dplyr::select(-CLIMATE, -REALM, -TAXA, -StudyMethod, -BIOMASS_TYPE) %>%
  # create incidence column
  mutate(incidence = ifelse(Biomass==0, 0, 1))

bm_coverage <- data.frame()
for (r in unique(bt_grid_collate_biomass$rarefyID)){
  dat <- bt_grid_collate_biomass[bt_grid_collate_biomass$rarefyID == r,]
  if(length(unique(dat$YEAR))<2){ next }
  dat_pa <- dat %>%
    dplyr::select(YEAR, Species, incidence) %>%
    complete(YEAR, Species, fill = list(incidence = 0)) %>%
    #each row is a unique cell and Species, each column is a year
    spread(YEAR, incidence)
  # how many singletons (present in only 1 year)    #FIXME: This is a *really* ugly way to count
  singletons <- length(which(rowSums(dat_pa[,2:ncol(dat_pa)])==1))
  # how many doubletons (present in only 2 years)
  doubletons <- length(which(rowSums(dat_pa[,2:ncol(dat_pa)])==2))
  # how many incidences in the whole matrix
  N <- length(which(dat_pa[,2:ncol(dat_pa)] == 1))
  # how many samples (years) in the matrix
  T <- ncol(dat_pa) - 1
  # Chat for incidence: eqn Chat_sample(T) from Table 2, Chao et al. 2014, similar to eqn 4a in Chao & Jost Ecology
  Chat_i <- 1 - (singletons/N) * (((T-1)*singletons)/((T-1)*singletons + 2*doubletons))
  # with fix from Chao et al 2014 Subroutine for iNEXT code (appendix)
  # correction for communities with no doubletons (this prevents NaN results for
  # either singletons==0 or doubletons==0)
  if (doubletons > 0) {
    Chat_corrected_i <- Chat_i  }
  else {
    Chat_corrected_i <- 1 - (singletons/N) * (((T-1)*singletons)/((T-1)*singletons + 2))  }
  calcs <- data.frame(rarefyID=r, singletons, doubletons, N, Chat_i, Chat_corrected_i)
  bm_coverage <- rbind(bm_coverage, calcs)
  rm(dat_pa)
}

mnb=mean(bm_coverage$Chat_corrected_i)
sdb=sd(bm_coverage$Chat_corrected_i)
# #plot the presence coverage data using the abundance-based mean and standard deviation as a guide
# ggplot(bm_coverage, aes(Chat_corrected_i)) + geom_histogram(aes(fill=Chat_corrected_i>=mn-sd), binwidth=0.05)

##### 3.4. Summary section ========================================
## Make a list of the rarefyIDs to keep for analysis based on coverage threshold (>= mn-sd of count data)
# NOTE: Count data will drop individual years that don't meet the criteria, and Presence and Biomass data will drop entire rarefyIDs 
countkeep <- unique(abund_coverage[abund_coverage$Chat_corrected>=mn-sd,c('rarefyID', 'YEAR')])
countkeep <- unite(countkeep, col=keep, rarefyID, YEAR, sep="_")
countkeep <- as.vector(countkeep$keep)

pakeep <- unique(pa_coverage[pa_coverage$Chat_corrected_i>=mnp-sdp,'rarefyID']) # presence
bmkeep <- unique(bm_coverage[bm_coverage$Chat_corrected_i>=mnb-sdb,'rarefyID']) # biomass

# Filter gridded studies for the coverage cutoff, prior to rarefaction
bt_count_filtered <- bt_grid %>%
  unite(col=keep, rarefyID, YEAR, sep="_", remove=FALSE) %>%
  filter(keep %in% countkeep) %>%
  dplyr::select(-keep)

bt_pabm_filtered <- bt_grid %>%
  filter(rarefyID %in% pakeep | rarefyID %in% bmkeep)

# DATASET filtered for coverage (>= mean-sd of count data) to be used for rarefaction
#   Note: count data was filtered at rarefyID + YEAR, and presence and biomass data was filtered at rarefyID
bt_grid_filtered <- rbind(bt_count_filtered, bt_pabm_filtered)

# #NOTE: bt_pabm_filtered does not exist because pa_coverage and bm_coverage have ZERO records
# bt_grid_filtered <- bt_count_filtered

#save(dgg, bt_grid_filtered2, file='BioTIME_grid_filtered2_181021.Rdata')

bt_grid_filtered %>% group_by()%>% 
  summarise(n_distinct(STUDY_ID),n_distinct(Species),n_distinct(rarefyID),n_distinct(paste(rarefyID,Species)),n_distinct(YEAR),min(YEAR),max(YEAR),n())-> bt_stats


#### 4.  Rarefy data ========================================

load('./BioTIME_grid_filtered.Rdata')
load('./BioTIME_grid.Rdata')

print(Sys.time())

##	use the job and task id as a counter for resampling and to set the seed
##	of the random number generator
uniq_id <- Sys.getenv('UNIQ_ID')
seed <- uniq_id
#seed <- substr(seed, 1, 4)
set.seed(seed)

# grid=bt_grid_abund
# # grid=bt_grid_abund[]
# type="count"
# resamples=1
# trimsamples=FALSE
# i=1
# j=1

#FUNCTION TO RAREFY DATA
rarefy_diversity <- function(grid, type=c("count", "presence", "biomass"), resamples=100, trimsamples=FALSE){

#	restrict calculations to where there is abundance>0 AND
#	following removal of NAs and 0's there are still more than 2 years

# grid: input dataset
# type: count, presence, or biomass
# resamples: the number of bootstrap resampling events desired (default is 100)
# trimsamples: TRUE means that years with < 1/2 the average number of samples should be removed to avoid excessive information loss. Default is FALSE.
# calculate the number of sampling events per year, and find the minimum
# resample the data to rarefy the diversity metrics
# output a new dataframe

if(type == "count" | type == "presence") { field = "Abundance" 
} else { field = "Biomass" }

# Get the sample-size to rarefy to. How many sampling events per cell per year?
# This is handled differently depending on the data type

# define a filter for 'field' to be >0 and not NA 
zero_NA_filter <- interp(~y > x & !is.na(y), 
        .values = list(y = as.name(field), x = 0))

#	define a function to calculate the sum(field) for use on the rarefied sample
sum_field <- interp(~sum(as.numeric(var), na.rm=T),
                    var= as.name(field))

nsamples <- ungroup(grid) %>%
  group_by(rarefyID, YEAR) %>%
  # remove 0's or NA's (this is to catch any places where abundance wasn't actually recorded 
  # or the species is indicated as absent)
  filter_(.dots=zero_NA_filter) %>%
  # calculate how many observations per year per study
  dplyr::summarise(nsamples = n_distinct(ObsEventID)) 

# Check if you wanted to remove years with especially low samples (< 1/2 the average number of samples)
if(trimsamples) {
  # Calculate the mean number of samples per cell
  mean_samp <- ungroup(nsamples) %>% group_by(rarefyID) %>%
    mutate(mean_samp = mean(nsamples), 
           lower_bound = mean(nsamples)/2) %>%
    filter(nsamples >= lower_bound)
  
  min_samp <- ungroup(mean_samp) %>% group_by(rarefyID) %>%
    mutate(min_samp = min(nsamples)) %>%
    # retain only the rows with the minimum sample size for a given cell
    filter(nsamples==min_samp) %>%
    distinct(rarefyID, min_samp, .keep_all=TRUE)
  
  # join the data and filter out years with  nsamples less than 1/2 the mean number of samples
  grid <- inner_join(grid, dplyr::select(min_samp, - YEAR, -nsamples)) %>%
    filter(n_samps >= lower_bound)
  rm(mean_samp,min_samp)
  print("trimsamples==TRUE: Removing years with < 1/2 the average number of samples for a given ID")
  
} else {
  # Calculate the minimum number of samples per cell
  min_samp <- ungroup(nsamples) %>% group_by(rarefyID) %>%
    mutate(min_samp = min(nsamples)) %>%
    # retain only the rows with the minimum sample size for a given cell
    filter(nsamples==min_samp) %>%
    distinct(rarefyID, min_samp, .keep_all=TRUE)
  
  #	Add the min_samp to the data and tidy a little
  grid <- inner_join(grid, dplyr::select(min_samp, - YEAR, -nsamples))
  rm(min_samp)
}

# Re-calculate metadata
new_meta2 <- ungroup(grid) %>%
  # remove 0's or NA's (this is to catch any places where abundance wasn't actually recorded 
  # or the species is indicated as absent)
  filter_(.dots=zero_NA_filter) %>%
  group_by(rarefyID,CLIMATE,REALM,TAXA,StudyMethod,ABUNDANCE_TYPE,BIOMASS_TYPE,STUDY_ID) %>%
  summarise(
    # total number of species in rarefyID time-series 
    SamplePool = n_distinct(Species),
    # total number of individuals
    SampleN = ifelse(type=='count', sum(as.numeric(Abundance)),
                     NA),
    # number of years sampled
    num_years = n_distinct(YEAR),
    # duration of time series, start and end points
    duration = max(YEAR) - min(YEAR) + 1,
    startYear = min(YEAR),
    endYear = max(YEAR)) 

#	Create dataframe where unique observations (i.e., the data of an
#	ObsEventID's [individual species abundances])
#	are nested within cells within years within studies 
bt_grid_nest <- ungroup(grid) %>%
  group_by(rarefyID, ObsEventID, cell, YEAR, min_samp) %>% #+biomass? - sp with same biomass get grouped
  # remove 0's or NA's (to catch any places where abundance wasn't actually recorded or the species is indicated as absent)
  filter_(.dots=zero_NA_filter) %>%
  # depending on type: nest(Species, Abundance) OR nest(Species, Biomass)
  nest_legacy(nest_cols=c("Species", "Abundance","Biomass")) %>%
  # reduce to studies that have more than two time points for a given cell
  group_by(rarefyID) %>%
  #keeps all study_cells with 2 or more years of data
  filter(n_distinct(YEAR)>=2) %>%  
  ungroup()

##	initialise df to store all data 
rare_samp_all_db <- data.frame()
#rare_samp_mean_db <-data.frame()

new<-Sys.time()
# ##	rarefy rarefy_resamps times
for(i in 1:resamples){
  #   
  uniq_id=i
  ## loop to do rarefaction for each study
  for(j in 1:length(unique(bt_grid_nest$rarefyID))){
  #for(j in 1:2){
    print(paste('rarefaction', i, 'out of', resamples, 'for study_cell', j, '(', unique(bt_grid_nest$rarefyID)[j], ')',  'in', length(unique(bt_grid_nest$rarefyID))))
    
    ##	get the jth study_cell
    study <- bt_grid_nest %>%
      filter(rarefyID==unique(bt_grid_nest$rarefyID)[j]) 
    # get minimum sample size for rarefaction	
    min_samp <- study %>% distinct(min_samp) %>% .$min_samp
    
    # check that there is only one cell represented (This shouldn't be a problem)
    if(length(unique(study$cell))>1) { 
      stop(paste0("ERROR: ", unique(study$rarefyID), " contains more than one grid cell")) }
    
    # check there there is more than one year in the cell
    if(length(unique(study$YEAR))<2) {
      print(paste0("ERROR: ", unique(study$rarefyID), " does not have more than one year"))
      next }
       
    rare_samp_all <- study %>%
      # rarefy to min_samp 
      group_by(rarefyID, YEAR) %>%
      sample_n(size=min_samp) %>%
      # unpack and collate taxa from rarefied sample	
      unnest(data) %>%
      # add unique counter for a resampling event
      mutate(rarefy_resamp = uniq_id) %>%  #
      ungroup()	
    
    rare_samp_all$resample<-uniq_id
    
      # add to dataframe for all studies
    rare_samp_all_db <- bind_rows(rare_samp_all_db, rare_samp_all)
    #rare_samp_mean_db <-bind_rows(rare_samp_mean_db, rare_samp_mean)
    
    #counter2 <- 1
    
    
  }	# rarefyID loop (STUDY_CELL ID)
  
}

print("rarefy done")

rare_samp_all_db_meta<- inner_join(new_meta2[,1:8], rare_samp_all_db, by='rarefyID')

old<-Sys.time()
print(old-new)

#write.csv(rare_samp_mean_db,file="./rare_samp_mean_68s3.csv")

  # rarefaction loop
  return(rare_samp_all_db_meta)
} # END function

## Calculate mean rarefied diversity for each data type=======================================
## Separate true abundance (count) data from presence and biomass/cover data
bt_grid_abund <- bt_grid %>%
  filter(ABUNDANCE_TYPE %in% c("Count", "Density", "MeanCount"))

bt_grid_pres <- bt_grid %>%
  filter(ABUNDANCE_TYPE == "Presence/Absence")

bt_grid_bmass <- bt_grid %>%
  filter(is.na(ABUNDANCE_TYPE)) #only want to calculate on biomass data when abundance data is also not available


## Get rarefied resample for each type of measurement (all data)
#bt_grid gets separated (for optimization)
bt_grid_abund_set1<-bt_grid_abund[bt_grid_abund$STUDY_ID<=120,]
bt_grid_abund_set2<-bt_grid_abund[bt_grid_abund$STUDY_ID>=121 & bt_grid_abund$STUDY_ID<=215,]
bt_grid_abund_set3<-bt_grid_abund[bt_grid_abund$STUDY_ID>=216,]

rarefy_abund1 <- rarefy_diversity(grid=bt_grid_abund_set1, type="count", resamples=1) #18204 - 18.58642 mins
rarefy_abund2 <- rarefy_diversity(grid=bt_grid_abund_set2, type="count", resamples=1) #31217 - 1.083978 hours
rarefy_abund3 <- rarefy_diversity(grid=bt_grid_abund_set3, type="count", resamples=1) #14200 - 28.80629 mins
rarefy_pres <- rarefy_diversity(grid=bt_grid_pres, type="presence", resamples=1) #4844 - 4.055046 mins
rarefy_biomass <- rarefy_diversity(grid=bt_grid_bmass, type="biomass", resamples=1) #4842 - 4.689915 mins

rare_samp_all<-rbind(rarefy_abund1,rarefy_abund2,rarefy_abund3,rarefy_pres,rarefy_biomass)

#export data
write.csv(rare_samp_all,file="./rare_samp_all.csv")

rare_samp_all %>% group_by(Species)%>% 
  summarise(n_distinct(TAXA),n())-> bt_stats

write.csv(bt_stats,file="./rare_samp_all_sp_list.csv")


