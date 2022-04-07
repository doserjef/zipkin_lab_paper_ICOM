# summary.R: this file summarizes results from the ICOM fit to a community
#            of interior forest obligate birds using eBird and BBS data. 
rm(list = ls())
library(coda)
library(nimble)
library(tidyverse)
# For plotting
library(sf)
library(stars)
library(viridis)
library(RColorBrewer)
library(ggpubr)

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}


# Load data and results ---------------------------------------------------
load("data/ne-data-bundle.rda")
load("results/ne-icom-1-chain-2022-04-06.R")
samples.1 <- samples
load("results/ne-icom-2-chain-2022-04-06.R")
samples.2 <- samples
load("results/ne-icom-3-chain-2022-04-06.R")
samples.3 <- samples

samples.list <- mcmc.list(samples.1, samples.2, samples.3)

# Assess convergence ------------------------------------------------------
# Number of species
N <- n_distinct(y.bbs$sp)
param.names <- unlist(dimnames(samples)[2])
gelman.diag(samples.list)

# Estimate species richness post-hoc --------------------------------------
# Covariate data used to fit the model
unique.cells <- sort(unique(c(unique(y.bbs$cell), unique(ebird.df$cell))))
J <- length(unique.cells)
occ.covs.small <- occ.covs[unique.cells, ]
# Extract species-specific occurrence parameters
curr.indx <- which(substr(param.names, 1, 5) == 'beta[')
beta.samples <- samples.list[, curr.indx]
# Intercept
beta.1.samples <- do.call('rbind', beta.samples[, 1:N])
# Linear elevation
beta.2.samples <- do.call('rbind', beta.samples[, (N + 1):(N * 2)])
# Quadratic elevation
beta.3.samples <- do.call('rbind', beta.samples[, (N * 2 + 1):(N * 3)])
# Forest cover
beta.4.samples <- do.call('rbind', beta.samples[, (N * 3 + 1):(N * 4)])
# Predict across all locations (including fitted)
occ.covs.pred <- occ.covs
# Get mean and sd for covariates used to fit the model
mean.elev.fit <- mean(occ.covs.small$elev)
sd.elev.fit <- sd(occ.covs.small$elev)
mean.pf.fit <- mean(occ.covs.small$pf)
sd.pf.fit <- sd(occ.covs.small$pf)
# Design matrix for model predictions
X.0 <- cbind(1, (occ.covs.pred$elev - mean.elev.fit) / sd.elev.fit, 
	     ((occ.covs.pred$elev - mean.elev.fit) / sd.elev.fit)^2, 
	     (occ.covs.pred$pf - mean.pf.fit) / sd.pf.fit)
# Set missing values to 0
X.0[which(is.na(X.0))] <- 0
colnames(X.0) <- c('intercept', 'elev', 'elev.2', 'forest')
# Number of prediction locations
J.0 <- nrow(X.0)
# Number of posterior samples
n.post <- nrow(beta.1.samples)
# Temporary occurrence probability for each species (not saved)
z.tmp <- matrix(NA, N, J.0)
# Species richness samples
rich.samples <- matrix(NA, n.post, J.0)
# Predict Richness --------------------
for (a in 1:n.post) {
  print(paste("Currently on iteration ", a, " out of ", n.post, sep = ''))
  for (i in 1:N) {
    z.tmp[i,] <- rbinom(J.0, 1, logit.inv(beta.1.samples[a, i] + 
			                  beta.2.samples[a, i] * X.0[, 2] + 
			                  beta.3.samples[a, i] * X.0[, 3] +
			                  beta.4.samples[a, i] * X.0[, 4]))
  } # i (species)
  rich.samples[a, ] <- apply(z.tmp, 2, sum)
} # a (iteration)

# Extract species richness estimates --------------------------------------
rich.mean <- apply(rich.samples, 2, mean, na.rm = TRUE)
rich.sd <- apply(rich.samples, 2, sd, na.rm = TRUE)
rich.low <- apply(rich.samples, 2, quantile, 0.025, na.rm = TRUE)
rich.high <- apply(rich.samples, 2, quantile, 0.975, na.rm = TRUE)

# Put in richness values directly into grid
grid.ne$rich.mean <- rich.mean
grid.ne$rich.sd <- rich.sd
grid.ne$ci.width <- rich.high - rich.low
grid.ne$elev <- occ.covs$elev
grid.ne$forest <- occ.covs$pf

# Number of eBird checklists in each cell
grid.ne$eBirdCheck <- 0
# Maximum number of checklists in a cell would be 9.
tmp.eb <- ebird.df %>%
  group_by(cell) %>%
  summarize(n.check = n_distinct(listID))
for (i in 1:nrow(grid.ne)) {
  if (i %in% tmp.eb$cell) {
    grid.ne$eBirdCheck[i] <- tmp.eb$n.check[which(tmp.eb$cell == i)]
  }
}
grid.ne$eBirdCheck <- as.integer(grid.ne$eBirdCheck)

# Get map of US -----------------------------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Full data
ne.states <- usa %>% 
  filter(ID %in% c('connecticut', 'delaware', 'maine', 'maryland', 
		     'massachusetts', 'new hampshire', 'new jersey', 
		     'new york', 'pennsylvania', 'rhode island', 
		     'vermont'))

# Map of the data used to fit the model -----------------------------------
# Load BBS locations
bbs.locs <- bbs.sf %>%
  filter(sp == 'HAWO')
grid.plot <- st_intersection(grid.ne, st_make_valid(ne.states))
counts.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = eBirdCheck, col = eBirdCheck)) + 
  geom_sf(data = ne.states, col = 'gray', fill = NA) + 
  scale_fill_gradientn(colors = viridis(10), na.value = NA) +
  scale_color_gradientn(colors = viridis(10), na.value = NA) +
  guides(col = "none") + 
  geom_sf(data = bbs.locs, pch = 21, fill = 'white', size = 2.5, col = 'black') +
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "# eBird Lists") +
  theme(legend.position = c(0.83, 0.19), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))
counts.plot
ggsave(device = 'png', filename = 'figures/ne-data.png', height = 6, width = 6, 
       units = 'in')
ggsave(device = 'pdf', filename = 'figures/ne-data.pdf', height = 6, width = 6, 
       units = 'in')

# Plot covariates ---------------------------------------------------------
# Elevation
elev.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = elev, col = elev)) + 
  geom_sf(data = ne.states, col = 'black', fill = NA) + 
  scale_fill_gradientn(colors = viridis(10), na.value = NA) +
  scale_color_gradientn(colors = viridis(10), na.value = NA) +
  theme_bw(base_size = 23) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") + 
  guides(col = 'none') +
  theme(legend.position = c(0.84, 0.17), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1))
elev.plot

forest.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = forest, col = forest)) + 
  geom_sf(data = ne.states, col = 'black', fill = NA) + 
  scale_fill_gradientn(colors = viridis(10), na.value = NA) +
  scale_color_gradientn(colors = viridis(10), na.value = NA) +
  theme_bw(base_size = 23) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") + 
  guides(col = 'none') +
  theme(legend.position = c(0.84, 0.17), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1))
forest.plot

# Species Richness Map ----------------------------------------------------
rich.mean.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = rich.mean, col = rich.mean)) + 
  geom_sf(data = ne.states, col = 'black', fill = NA) + 
  # scale_fill_gradientn(colors = viridis(10), na.value = NA) +
  # scale_color_gradientn(colors = viridis(10), na.value = NA) +
  scale_fill_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  scale_color_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  theme_bw(base_size = 23) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") + 
  guides(col = 'none') +
  theme(legend.position = c(0.84, 0.17), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1))
rich.mean.plot
ggsave(device = 'png', filename = 'figures/ne-richness.png', height = 7, width = 10, 
       units = 'in')
ggsave(device = 'pdf', filename = 'figures/ne-richness.pdf', height = 7, width = 10, 
       units = 'in')
# Map of CI width
rich.ci.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = ci.width, col = ci.width)) + 
  geom_sf(data = ne.states, col = 'black', fill = NA) + 
  # scale_fill_gradientn(colors = viridis(10), na.value = NA) +
  # scale_color_gradientn(colors = viridis(10), na.value = NA) +
  scale_fill_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  scale_color_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  theme_bw(base_size = 25) +
  guides(col = 'none') + 
  labs(x = "Longitude", y = "Latitude", fill = "CI Width") +
  theme(legend.position = c(0.84, 0.17), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1))
rich.ci.plot
ggsave(device = 'png', filename = 'figures/ne-richness-ci.png', height = 7, width = 10, 
       units = 'in')
ggsave(device = 'pdf', filename = 'figures/ne-richness-ci.pdf', height = 7, width = 10, 
       units = 'in')


# Plot of covariate effects -----------------------------------------------
param.names <- attr(samples, 'dimnames')[[2]]
curr.indx <- which(substr(param.names, 1, 4) == 'beta')
# Samples of species and community-level regression coefficients
beta.samples <- do.call("rbind", samples.list[, curr.indx])
# Order: elevation, elevation^2, forest cover. 
# Hardcoding because I'm lazy
# Includes species-specific effects of forest cover as well as community-level. 
# beta.for.samples <- beta.samples[, c(55:81, ncol(beta.samples))]
beta.for.samples <- beta.samples[, (N * 3 + 1):(N * 4)]
beta.comm.for.samples <- beta.samples[, ncol(beta.samples)]
# Species codes, plus the community code (COMM)
sp.code.factor <- factor(unique(y.bbs$sp), 
		         levels = unique(y.bbs$sp))
sp.codes <- as.numeric(sp.code.factor)
N <- length(sp.codes)
cov.plot.df <- data.frame(for.mean = apply(beta.for.samples, 2, mean), 
			  for.low = apply(beta.for.samples, 2, quantile, 0.25), 
			  for.lowest = apply(beta.for.samples, 2, quantile, 0.025), 
			  for.high = apply(beta.for.samples, 2, quantile, 0.75), 
			  for.highest = apply(beta.for.samples, 2, quantile, 0.975), 
			  sp = sp.codes)
# Rearrange and add things to get the plot to display effects in increasing
# order. 
cov.plot.df <- cov.plot.df %>%
  arrange(for.mean)
cov.plot.df$sp.factor <- as.character(sp.code.factor[cov.plot.df$sp])
cov.plot.df$sort.sp <- 1:N

# Add in the community level covariate
comm.plot.df <- data.frame(for.mean = mean(beta.comm.for.samples),
			   for.low = quantile(beta.comm.for.samples, 0.25), 
			   for.lowest = quantile(beta.comm.for.samples, 0.025), 
			   for.high = quantile(beta.comm.for.samples, 0.75), 
			   for.highest = quantile(beta.comm.for.samples, 0.975), 
			   sp = 28, 
			   sp.factor = 'COMM',
                           sort.sp = N + 1)
cov.plot.df <- rbind(cov.plot.df, comm.plot.df)
cov.plot.df$sp.factor <- factor(cov.plot.df$sp.factor, levels = c(unique(cov.plot.df$sp.factor)))


for.cov.plot <- ggplot(data = cov.plot.df, aes(x = sort.sp, fill = for.mean, group = sp.factor)) + 
  geom_hline(yintercept = 0, col = 'black', size = 0.75, lty = 2) + 
  geom_vline(xintercept = 27.5, col = 'black', size = 0.5, lty = 1) + 
  geom_boxplot(aes(ymin = for.lowest, lower = for.low, middle = for.mean, 
		   upper = for.high, ymax = for.highest), stat = 'identity', col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_classic(base_size = 18) + 
  guides(fill = "none") + 
  labs(x = "Species", y = "Effect of Forest Cover") + 
  scale_x_continuous(breaks = 1:(N+1), labels = cov.plot.df$sp.factor) +  
  coord_flip()
for.cov.plot
ggsave(device = 'png', filename = 'figures/ne-forest-cover.png', height = 6, width = 6, 
       units = 'in')
ggsave(device = 'pdf', filename = 'figures/ne-forest-cover.pdf', height = 6, width = 6, 
       units = 'in')

ggarrange(rich.mean.plot, for.cov.plot, nrow = 1, ncol = 2, 
          labels = c('(A)', '(B)'), 
          font.label = list(size = 18))
ggsave(device = 'pdf', filename = 'figures/icom-main-figure.pdf', height = 7, width = 12, 
       units = 'in')
ggsave(device = 'png', filename = 'figures/icom-main-figure.png', height = 7, width = 12, 
       units = 'in')
