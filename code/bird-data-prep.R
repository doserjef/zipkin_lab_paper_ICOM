# bird-data-prep.R: this script prepares the eBird and BBS data for analysis
#                   in an ICOM framework. 
rm(list = ls())
library(tidyverse)
library(lubridate)
library(sf)
# For reading in eBird data
library(auk)
# For US states data as sf objects
library(spData)
# For linking bird codes with species info. 
library(wildlifeR)
# For some exploratory plots
library(ggnewscale)
# For grabbing NED and NLCD data
library(FedData)
# Need the old guys for NED and NLCD data
library(sp)
library(rgdal)
library(raster)

# Get BBS Data ------------------------------------------------------------
# Read in BBS Data from NE states
# bbs.full.dat <- list.files(path = "data/BBS/ne-states", full.names = TRUE) %>%
#   lapply(read.csv) %>%
#   bind_rows()
# Only use PA for testing purposes. 
bbs.full.dat <- read.csv("data/BBS/ne-states/Pennsyl.csv")
# Get associated route data
route.dat <- read.csv("data/BBS/routes.csv")
# Get associated weather data
weather.dat <- read.csv("data/BBS/weather.csv")
# Only grab weather data from NE states
weather.dat <- weather.dat %>%
  filter(StateNum %in% unique(bbs.full.dat$StateNum))
# Join BBS data with route data
bbs.dat <- left_join(bbs.full.dat, route.dat, by = c('Route', 'CountryNum', 'StateNum'))
# Only grab data from one year (using 2017 for now, it's totally arbitrary)
bbs.2017 <- bbs.dat %>%
  filter(Year == 2017)
# Select columns of interest
bbs.dat <- bbs.dat %>%
  dplyr::select(Route, Year, AOU, starts_with('Stop'), Latitude, Longitude)
# Select columns of interest
bbs.2017 <- bbs.2017 %>%
  dplyr::select(RouteDataID, Latitude, Longitude, AOU, starts_with("Count")) %>%
  dplyr::select(-CountryNum)
# Fill in implicit zeros
bbs.2017 <- bbs.2017 %>%
  complete(AOU, nesting(RouteDataID, Latitude, Longitude))
# Replace NAs with 0s for all columns at once.
bbs.2017 <- bbs.2017 %>%
  replace(is.na(.), 0)
# Link AOU codes with species info. The AOU_species_codes object is 
# from wildlifeR. 
aou.info <- AOU_species_codes %>%
  mutate(AOU = spp.num)
bbs.2017 <- left_join(bbs.2017, aou.info, by = c("AOU"))

# Read in guild information -----------------------------------------------
guild.info <- read.table("data/guild-info.txt")
interior.for <- guild.info %>%
  dplyr::filter(Response_Guild == 'InteriorForestObligate') %>%
  pull(AOU_Code)

# Get detection covariates in proper format -------------------------------
weather.dat <- weather.dat %>% 
  unite('date', sep = '-', Year, Month, Day, remove = FALSE)
weather.dat$date <- as.Date(weather.dat$date, tz = "America/New_York")
# Get julian date of each survey
weather.dat$julian <- as.numeric(format(weather.dat$date, '%j'))
# Join detection covariates with observations
bbs.2017 <- left_join(bbs.2017, weather.dat, by = c('RouteDataID'))

# Filter out for only interior forest obligates ---------------------------
bbs.2017 <- bbs.2017 %>%
  separate(spp, c('orderName', 'sppName'), sep = ' ', remove = FALSE) %>%
  mutate(orderName = tolower(orderName))
bbs.2017 <- bbs.2017 %>%
  filter(alpha.code %in% interior.for)
# Get observation data in long format
y.bbs.long <- bbs.2017 %>%
  dplyr::select(alpha.code, starts_with('Count'), Route, Longitude, Latitude, 
		ObsN, julian) %>%
  mutate(sp = as.character(alpha.code)) %>%
  dplyr::select(-alpha.code, -CountryNum)
n.route <- n_distinct(y.bbs.long$Route)
n.sp <- n_distinct(y.bbs.long$sp)
n.visit <- 5
y.bbs.long <- y.bbs.long %>%
  mutate(across(starts_with('Count'), function(a) {ifelse(a > 0, 1, 0)}))
# Add summary binomial count. 
y.bbs.long$binom <- apply(y.bbs.long[, c('Count10', 'Count20', 'Count30', 
					 'Count40', 'Count50')], 1, sum)
# Raw occurrence probabilities at the route level. 
occ.raw.bbs <- y.bbs.long %>%
  mutate(prop = binom / n.visit) %>%
  group_by(sp) %>%
  summarize(avg.prop = mean(prop)) %>%
  arrange(desc(avg.prop)) %>%
  print(n = nrow(.))

# Get eBird data ----------------------------------------------------------
# This is for testing right now. Just PA for now. 
ebird.ne <- read_ebd('data/eBird-test/ebd_US-PA_201705_201707_relSep-2021.txt')
# Filter for records in May and June
ebird.ne <- ebird.ne %>%
  mutate(month = month(observation_date)) %>%
  filter(month %in% c(5, 6))

# Filter eBird data set to only those that have auxiliary information
# and not extreme information
ebird.dat <- ebird.ne %>%
  filter(month %in% c(5, 6), # May and June closure is assumed
	 protocol_type %in% c('Stationary', 'Traveling'), # stationary or traveling protocols
	 all_species_reported == TRUE, # complete checklists
	 number_observers < 10, # restrict number of observers
	 duration_minutes < 300) %>% # remove super long surveys
  mutate(effort_distance_km = ifelse(protocol_type == 'Stationary', 0,
				     effort_distance_km)) %>%
  filter(effort_distance_km < 5) # restrict spatial scale. This should inform the size of grid cell.

# Select columns of interest, create new columns of interest
ebird.dat <- ebird.dat %>%
  dplyr::select(checklist_id, common_name, observation_count, latitude, longitude,
	 observation_date, time_observations_started, observer_id,
	 duration_minutes, effort_distance_km, month,
	 number_observers, all_species_reported) %>%
  mutate(year = year(observation_date),
	 week = week(observation_date))

# Grid up study area ------------------------------------------------------
data(us_states)
# Full data
# ne.states <- us_states %>% 
#   filter(NAME %in% c('Connecticut', 'Delaware', 'Maine', 'Maryland', 
# 		     'Massachusetts', 'New Hampshire', 'New Jersey', 
# 		     'New York', 'Pennsylvania', 'Rhode Island', 
# 		     'Vermont'))
pa.states <- us_states %>% 
  filter(NAME %in% c('Pennsylvania'))
pa.states <- pa.states %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Can change the grid cell size as needed. 
grid <- pa.states %>%
  st_make_grid(cellsize = c(5000, 5000)) # in meters

# Connnect eBird to the gridded area --------------------------------------
ebird.sf <- st_as_sf(ebird.dat, coords = c('longitude', 'latitude'), 
		     crs = st_crs(us_states))
# Convert ebird locations to albers equal area
ebird.sf <- ebird.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Convert grid object to sf data frame
grid.sf <- grid %>% st_as_sf()
# Get cell each eBird observation falls in
cell.eBird <- st_contains(grid.sf, ebird.sf)
# Get number of eBird observations in each cell
cell.counts <- sapply(cell.eBird, length)
# Number of cells
n.cells <- length(cell.counts)
# Cell corresponding to each value in cells.eBird.vec
cell.by.val <- unlist(sapply(1:n.cells, function(a) {rep(a, cell.counts[a])}))
# Rows in eBird that are in the grid. 
cells.eBird.vec <- unlist(cell.eBird)
grid.sf$counts <- cell.counts
grid.pa <- grid.sf %>%
  st_intersection(pa.states)
# Assign each eBird observation to a cell
ebird.dat$cell <- NA
ebird.dat$cell[cells.eBird.vec] <- cell.by.val
# Remove records not in the study area
ebird.dat <- ebird.dat %>%
  filter(!is.na(cell))

# Connect BBS data to the gridded area ------------------------------------
# Convert BBS data to sf object.
bbs.sf <- st_as_sf(y.bbs.long, coords = c('Longitude', 'Latitude'),
		     crs = st_crs(us_states))
# Convert to albers equal area
bbs.sf <- bbs.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# Get cell each BBS data point falls in
cell.bbs <- st_contains(grid.sf, bbs.sf)
bbs.sf$cell <- NA
for (i in 1:n.cells) {
  if (length(cell.bbs[[i]]) > 0) {
    bbs.sf$cell[cell.bbs[[i]]] <- i    
  }
}
# Total possible counts (denominator of binomial)
bbs.sf$denom <- 5

# BBS cells regardless of species
bbs.sf.all <- bbs.sf %>%
  group_by(Route) %>% 
  summarize(all.count = sum(binom))

# Some exploratory plots --------------------------------------------------
# eBird sampling bias. 
# This is obviously related to urbanization. Should include urbanization as a 
# covariate on the preferential sampling model and the detection model. 
# This is currently for all species, as you need to retain all species 
# to create the preferential sampling part of the model. 
# Doesn't really show anything. 
ggplot() + 
  geom_sf(data = grid.pa, aes(fill = counts)) + 
  geom_sf(data = pa.states, alpha = 0, col = 'black') + 
  scale_fill_viridis_c() + 
  labs(fill = "eBird Counts") +
  new_scale('fill') +
  geom_sf(data = bbs.sf.all, pch = 21, aes(fill = all.count), size = 4) + 
  scale_fill_viridis_c() + 
  labs(fill = "BBS Counts") + 
  theme_bw()

# Get covariates within each grid cell ------------------------------------
# Split up grabbing the covariates into multiple data sources. 
coords.sp <- grid.sf %>% 
  st_coordinates() %>%
  as.data.frame()
J <- n_distinct(coords.sp$L2)
n.coords <- nrow(coords.sp)
coordinates(coords.sp) <- ~X + Y
proj4string(coords.sp) <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'
# Loop through all the sites
vals <- split(1:J, ceiling(seq_along(1:J)/6))
# Average elevation within a cell
elev.btbw <- rep(0, J)
# Proportion of forest cover within the cell
for.btbw <- rep(0, J)
# Proportion of developed land in the cell
devel.btbw <- rep(0, J)
# Get proportion of forest
props <- function(a, na.rm = TRUE) {
  my.sum <- sum(!is.na(a))	
  prop.for <- sum(a %in% c(41, 42, 43), na.rm = na.rm) / my.sum
  return(prop.for)
}
props.developed <- function(a, na.rm = TRUE) {
  my.sum <- sum(!is.na(a))	
  prop.for <- sum(a %in% c(21, 22, 23, 24), na.rm = na.rm) / my.sum
  return(prop.for)
}
ned.dat <- list()
nlcd.dat <- list()
# Takes a day or so to run this. 
# for (i in 1:length(vals)) {
#   print(paste("Currently on iteration ", i, " out of ", length(vals), sep = ''))
#   indx <- which(coords.sp$L2 %in% vals[[i]])
#   ned.dat[[i]] <- get_ned(template = coords.sp[indx, ], label = paste('birds', i))
#   nlcd.dat[[i]] <- get_nlcd(template = coords.sp[indx, ], label = paste('birds', i), year = 2016)
#   elev.btbw[vals[[i]]] <- raster::extract(ned.dat[[i]], grid.sf[vals[[i]], ], fun = mean)
#   for.btbw[vals[[i]]] <- raster::extract(nlcd.dat[[i]], grid.sf[vals[[i]], ], fun = props)
#   devel.btbw[vals[[i]]] <- raster::extract(nlcd.dat[[i]], grid.sf[vals[[i]], ], fun = props.developed)
# }
# 
# # Save covariates so you don't have to rerun all the time.
# #save(elev.btbw, for.btbw, devel.btbw, file = "data/spatial-covariates.rda")

# Load in the spatial covariates
load("data/spatial-covariates.rda")
# This brings in elev.btbw, for.btbw, and devel.btbw
# 
# # Get BBS data in analysis format -----------------------------------------
y.bbs <- bbs.sf %>% 
  as.data.frame() %>%
  dplyr::select(binom, sp, cell, obs = ObsN, julian) %>%
  arrange(sp)
# Keep data in long format. Will be faster for analysis. 
# 
# Get eBird Data in analysis format ---------------------------------------
# Number of species
N <- n_distinct(y.bbs$sp)
# Julian Date. 
ebird.dat$julian <- as.numeric(format(ebird.dat$observation_date, '%j'))
# Time survey started in minutes since midnight.  
ebird.dat$time.started <- ebird.dat$time_observations_started %>%
  hms() %>%
  period_to_seconds() / 60

# Spatio-temporal subset the eBird data -----------------------------------
# Select one checklist per week/cell combination to reduce the effects
# of preferential sampling on the resulting estimates
ebird.dat <- ebird.dat %>%
  unite(col = "weekCell", c(week, cell), remove = FALSE)

n.combos <- n_distinct(ebird.dat$weekCell)
week.cell.combos <- unique(ebird.dat$weekCell)
curr.indx <- vector(mode = 'list', length = n.combos)
# Takes a minute or so
for (i in 1:n.combos) {
  tmp <- ebird.dat %>%
    filter(weekCell == week.cell.combos[i])
  check.tmp <- unique(tmp$checklist_id)
  curr.val <- sample(1:length(check.tmp), 1)
  curr.indx[[i]] <- which(ebird.dat$checklist_id == check.tmp[curr.val])
} # i

# Filter out eBird data
ebird.filt.dat <- ebird.dat[unlist(curr.indx), ]

checklist.data <- ebird.filt.dat %>%
  group_by(cell) %>%
  summarize(n.lists = n_distinct(checklist_id)) %>%
  ungroup()

# Get unique checklist id in each cell. 
checklist.by.cell <- ebird.filt.dat %>%
  group_by(cell) %>%
  summarize(listID = unique(checklist_id)) %>%
  ungroup()

checklist.by.cell <- checklist.by.cell %>%
  slice(rep(row_number(), N))
# Total number of data set rows (species x rep/site combos)
n.ebird <- nrow(checklist.by.cell)

# Species in alphabetical order of 4-letter code
sp <- unique(y.bbs$sp)
sp.names <- aou.info %>%
  filter(alpha.code %in% sp) %>%
  arrange(alpha.code) %>%
  pull(name) %>%
  as.character()

# Create data frame to hold data and auxiliary information
ebird.df <- checklist.by.cell
ebird.df$sp <- rep(sp, times = sum(checklist.data$n.lists))
ebird.df$day <- NA
ebird.df$time <- NA
ebird.df$length <- NA
ebird.df$dist <- NA
ebird.df$obsv <- NA
ebird.df$week <- NA
ebird.df$y <- NA

# Get all eBird data for analysis. This takes a couple of minutes to run. 
for (i in 1:(n.ebird/N)) {
  if(i %% 1000 == 0) print(i)
  tmp <- ebird.filt.dat %>%
    filter(cell == checklist.by.cell$cell[i],
	   checklist_id == checklist.by.cell$listID[i])
  low <- (i - 1) * N + 1
  high <- i * N
  ebird.df$day[low:high] <- max(tmp$julian)
  ebird.df$time[low:high] <- max(tmp$time.started)
  ebird.df$week[low:high] <- max(tmp$week)
  ebird.df$length[low:high] <- max(tmp$duration_minutes)
  ebird.df$dist[low:high] <- max(tmp$effort_distance_km)
  ebird.df$obsv[low:high] <- max(tmp$number_observers)
  tmp.2 <- tmp %>% filter(str_to_upper(common_name) %in% str_to_upper(sp.names))
  curr.indices <- which(str_to_upper(sp.names) %in% str_to_upper(tmp.2$common_name))
  curr.vals <- ifelse(1:N %in% curr.indices, 1, 0)
  ebird.df$y[low:high] <- curr.vals
}


# Save results ------------------------------------------------------------
occ.covs <- data.frame(elev = elev.btbw, 
		       pf = for.btbw, 
		       devel = devel.btbw)
save(y.bbs, ebird.df, occ.covs, grid.sf, file = "data/pa-data-bundle.rda")

