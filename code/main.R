# main.R: main script for calling NIMBLE to run the ICOM for a community of 
#         interior forest obligates using eBird and BBS data. 
# Author: Jeffrey W. Doser

rm(list = ls())
library(tidyverse)
library(coda)
library(sf)
library(nimble)

load("data/pa-data-bundle.rda")
source("code/icom-ne-birds-nimble.R")

# Data prep ---------------------------------------------------------------
# Total number of cells in the study area. 
J <- nrow(grid.sf)
# Number of cells with eBird data
J.ebird <- n_distinct(ebird.df$cell)
# Number of cells with BBS data
J.bbs <- n_distinct(y.bbs$cell)
# Number of species
N <- n_distinct(ebird.df$sp)
# Cells between the eBird and BBS data don't currently match up. 

# Get chain number from command line run ----------------------------------
# This is used to save the chain number in the resulting file name.
# chain <- as.numeric(commandArgs(trailingOnly = TRUE))
# For testing
chain <- 1
if(length(chain) == 0) base::stop('Need to tell NIMBLE the chain number')

# Some EDA ----------------------------------------------------------------
grid.sf$elevation <- occ.covs$elev
grid.sf$forest <- occ.covs$pf
grid.sf$devel <- occ.covs$devel

# ebird.df %>%
#   group_by(sp) %>%
#   summarize(n.obs = sum(y)) %>%
#   print(n = nrow(.))
# 
# y.bbs %>%
#   group_by(sp) %>%
#   summarize(n.obs = sum(ifelse(binom > 0, 1, 0))) %>%
#   print(n = nrow(.))

# There is way more eBird data than BBS data. 
# Get Data Prepped for NIMBLE ---------------------------------------------
# Occurrence design matrix
X <- cbind(c(scale(occ.covs$elev)), 
	   c(scale(occ.covs$elev)^2), 
	   c(scale(occ.covs$pf)), 
           c(scale(occ.covs$devel)))
# Number of occurrence parameters
p.occ <- ncol(X)
# Total number of cells
J <- nrow(grid.sf)
# Prep BBS Data------------------------
bbs.df <- y.bbs
# The actual data
y.bbs <- bbs.df$binom
# Unique species names
sp.names <- unique(bbs.df$sp)
# Total number of species
N <- length(sp.names)
# Species index
sp.indx.bbs <- as.numeric(factor(bbs.df$sp))
# Cell index
cell.bbs <- bbs.df$cell
# Observer index
obs.indx <- as.numeric(factor(bbs.df$obs))
n.obs.bbs <- n_distinct(obs.indx)
# Fixed effects detection design matrix
X.bbs <- cbind(c(scale(bbs.df$julian)), c(scale(bbs.df$julian)^2))
# Number of BBS detection fixed effects
p.det.bbs <- ncol(X.bbs)
# Total number of BBS data points
n.vals.bbs <- nrow(X.bbs)
# Prep eBird Data ---------------------
# The actual data
y.eb <- ebird.df$y
# Species index
sp.indx.eb <- as.numeric(factor(ebird.df$sp))
# Cell index
cell.eb <- ebird.df$cell
# Fixed effects detection design matrix
X.eb <- cbind(c(scale(ebird.df$day)), 
	      c(scale(ebird.df$day)^2), 
	      c(scale(ebird.df$time)), 
	      c(scale(ebird.df$length)), 
	      c(scale(ebird.df$dist)), 
	      c(scale(ebird.df$obsv))) 
# Number of eBird detection fixed effects
p.det.eb <- ncol(X.eb)
# Total number of eBird data points
n.vals.eb <- nrow(X.eb)

# Constants ---------------------------------------------------------------
icom.consts <- list(p.occ = p.occ, p.det.bbs = p.det.bbs, p.det.eb = p.det.eb, 
		    N = N, n.obs.bbs = n.obs.bbs, J = J, X = X, X.bbs = X.bbs, 
		    sp.indx.bbs = sp.indx.bbs, obs.indx = obs.indx, 
		    cell.bbs = cell.bbs, n.vals.bbs = n.vals.bbs, X.eb = X.eb,
		    sp.indx.eb = sp.indx.eb, cell.eb = cell.eb, n.vals.eb = n.vals.eb)
# Data --------------------------------------------------------------------
icom.data <- list(y.bbs = y.bbs, y.eb = y.eb)
# Initial values ----------------------------------------------------------
z.init <- matrix(1, N, J)
icom.inits <- list(z = z.init, beta.comm = rnorm(p.occ), 
		   sigma.sq.beta = runif(p.occ, 0.5, 3),
		   alpha.comm.bbs = rnorm(p.det.bbs), 
		   sigma.sq.bbs = runif(p.det.bbs + 1, 0.5, 3),
		   alpha.comm.eb = rnorm(p.det.eb), 
		   sigma.sq.eb = runif(p.det.eb, 0.5, 3))
# Create the model --------------------------------------------------------
icom.model <- nimbleModel(code = icom.code, name = 'icom', constants = icom.consts,
		          data = icom.data, inits = icom.inits)
# Configure MCMC ----------------------------------------------------------
icom.conf <- configureMCMC(icom.model, monitors = c('beta.comm', 'sigma.sq.beta', 
						    'alpha.comm.bbs', 'sigma.sq.bbs', 
						    'alpha.comm.eb', 'sigma.sq.eb', 
						    'beta', 'alpha.bbs', 'alpha.eb'))
# Create an MCMC function -------------------------------------------------
icom.mcmc <- buildMCMC(icom.conf)
# Compile model
icom.c.model <- compileNimble(icom.model)
icom.c.mcmc <- compileNimble(icom.mcmc, project = icom.model)
# Number of iterations --------------------------------------------------
n.iter <- 3000
n.burn <- 1000
n.thin <- 2
n.chain <- 1
samples <- runMCMC(icom.c.mcmc, niter = n.iter, nburnin = n.burn,
	           thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE)

# Save results ------------------------------------------------------------
save(samples, file = paste("results/pa-icom-", chain, "-chain-", 
		           Sys.Date(), ".R", sep = ''))

