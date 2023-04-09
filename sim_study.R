#!/usr/bin/env Rscript

# 1. Initialization ----
time.start <- Sys.time()
set.seed(1)

## i) Load libraries ----
library(INLA)
library(raster)
library(rgdal)
library(TMB)
library(sf)
library(geoR)
library(haven)
library(geosphere)
library(scoringRules)
library(locfit)

## ii) Get utilities ----
setwd("~/Github/donut/") # Your working directory
source("functions.R")
source("makeIntegrationPoints.R")
source("donut_integration.R")

# 2. Setup ----
## i) Sub-integration points ----
nSubRPerPoint <- 10
nSubAPerPoint <- 10

## ii) Load transformed covariate rasters ----
pop <- raster("cov_rasters/pop_min_max.tif")
topo <- raster("cov_rasters/topo_min_max.tif")
water <- raster("cov_rasters/water_min_max.tif")

cov.raster.list <- list(pop, topo, water) # list of covariates
cov.raster.list <- NULL # Set to NULL if no covaraites

## iii) Specify jittering schemes as a dataframe ----

# Only donut jittering with inner radius 1 km and outer radius 3 km.
jittering.schemes <- data.frame( type = c("d"),
                                 R1 = c(1),
                                 R2 = c(3)
)


# 3. Simulate data ----
prep.start <- Sys.time()

## i) Get  Geographic data ----

# Masterframe of enumeration areas (EAs) which is used to sample clusters
master.frame <- readRDS("EAs_pop4.Rda")

# Administrative boudnaries
nigeria.adm2 <- readOGR("gadm41_NGA_shp/gadm41_NGA_2.shp")


## ii) Create meshes ----
sample <- readRDS("EAs_pop4_second.Rda") # Sample locations used to construct mesh
mesh <- inla.mesh.2d(loc.domain = cbind(sample$east, sample$north),
                                n = 1000, 
                                max.n = 10000,
                                offset = -.08,
                                cutoff = 40, 
                                max.edge = c(60, 120))

## iii) Simulate response and create inputs ----

# Set parameters
range <- 180 # km
sample.size <- 1400 # Number of observations
beta.scale <- 4 # The magnitude of the covariates, which all have same magnitude
beta.vec <- beta.scale * rep(1, length(cov.raster.list) + 1)
sigma2.S = 1 # Spatial variance
alpha.sp = 2
sigma2.N <- 0.1 # Nugget variance
num.datasets <- 1 # Number of datasets for each jittering scheme
true.params <- c(range, sigma2.S, sigma2.N, beta.scale)

# Import prediction grid (Constant across scenarios)
pred.locs <- readRDS("pred_grid_km_v2.Rda")
n.preds <- nrow(pred.locs)


# Simulate observation locs. and responses.
obs.data.list <- list()
obs.locs.list <- list()
pred.data.list <- list()
true.locs.list <- list()

# Make list of inputs also 
inputs <- list()

# Set priors
matern.pri <- c(range, 0.50, 1., 0.50) # c(a, b, c, d) P(range < a) = b; P(sigma > c) = d
beta.pri <- c(0, 10) # N(mean, sigma)
nugget.pri <- c(0.968, 0.01) # c(U, a), P(sigma_N > U) = a
prior.list <- list(
  matern.pri = matern.pri,
  beta.pri = beta.pri,
  nugget.pri = nugget.pri
)
model.params <- list(alpha = alpha.sp)

# Simulate datasets
for(i in 1:num.datasets){
  obs.locs.temp <- list()
  inputs.temp <- list()
  
  # Sample EAs uniformly
  EA.idxs <- sample.int(nrow(master.frame), sample.size, replace = FALSE)
  true.locs.lat.lon <- master.frame[EA.idxs, c("x", "y")]
  true.locs.km <- convert_deg_to_km(true.locs.lat.lon)
  colnames(true.locs.km)[1:2] <- c("east", "north") 
  true.locs.list[[i]] <- true.locs.km
  
  # Simulate responses
  sim.locs <-rbind(true.locs.list[[i]], pred.locs)
  pred.idxs <- (nrow(true.locs.list[[i]]) + 1):(nrow(true.locs.list[[i]]) + nrow(pred.locs))
  sim.data <- spde.sim(sim.locs, kappa = sqrt(8)/range, variance = sigma2.S, alpha = alpha.sp, mesh = mesh, 
                       beta.vec = beta.vec, cov.raster.list, nug.std = sqrt(sigma2.N), test.idxs = pred.idxs)
  obs.data <- sim.data[[1]]
  pred.data <- sim.data[[2]]
  obs.data.list[[i]] = obs.data
  pred.data.list[[i]] = pred.data
  
  for(j in 1:nrow(jittering.schemes)){
    scheme <- jittering.schemes[j, ]
    true.lat.lon <- SpatialPoints(true.locs.lat.lon, proj4string = nigeria.adm2@proj4string, bbox = NULL)
    check1 <- over(true.lat.lon, nigeria.adm2, returnList = FALSE)
    if(scheme[1] == "rp"){
      loc.obs <- jitter.locs(scale = c(as.numeric(scheme[3])/2),
                             locKM = true.locs.km,
                             urbanRural = rep("U", nrow(true.locs.km)),
                             shapeFile = nigeria.adm2,
                             check1 = check1,
                             boundary = TRUE)
      loc.obs <- loc.obs[[1]] # Fix formatting.
    } else if( scheme[1] == "d"){
      loc.obs <- jitter.donut(true.locs.km, 
                              R1 = as.numeric(scheme[2]), 
                              R2 = as.numeric(scheme[3]), 
                              shapeFile = nigeria.adm2, 
                              check1 = check1, 
                              boundary = TRUE)
      
    } else{
      stop("Unknown jittering scheme.")
    }
    # Make TMB inputs
    donut.info <- NULL
    if(scheme[1] == "d"){donut.info <- list(num.pts.per.ring = 5,
                                            num.rings = 2,
                                            R1 = as.numeric(scheme[2]),
                                            R2 = as.numeric(scheme[3]))}
    tmb.input <- prepare.input(sim.data = obs.data,
                               loc.obs = as.matrix(loc.obs),
                               model.params = model.params, # Need this argument?
                               priors = prior.list,
                               jittering.scale = as.numeric(scheme[3])/2,
                               urban = rep("U", nrow(true.locs.km)),
                               admin.map = nigeria.adm2,
                               mesh.s = mesh,
                               cov.raster.list,
                               nSubRPerPoint = nSubRPerPoint, 
                               nSubAPerPoint = nSubAPerPoint,
                               beta.random = TRUE,
                               jittering.scheme = scheme[1],
                               donut.info = donut.info,
                               testMode = FALSE)
    inputs.temp[[j]] <- tmb.input
    obs.locs.temp[[j]] <- loc.obs
  }
  inputs[[i]] <- inputs.temp
  obs.locs.list[[i]] <- obs.locs.temp
}

## Save inputs ----
all.inputs <- list(inputs = inputs,
                   obs.data = obs.data.list,
                   obs.locs = obs.locs.list,
                   pred.data = pred.data.list,
                   true.locs = true.locs.list,
                   jittering.schemes = jittering.schemes,
                   true.params = true.params)

saveRDS(all.inputs, file = "all_inputs.RData")

prep.end <- Sys.time()

# 4. Compile C++ files ----

## i. Model which accounts for random perturbation (disk) jittering ----
compile("jitter.cpp", "-O0 -g")
dyn.load(dynlib("jitter"))

## ii. Model which accounts for random perturbation jittering where only urban clusters are present ----
compile("no_rural.cpp", "-O0 -g")
dyn.load(dynlib("no_rural"))

## iii. Model which accounts for donut jittering ----
compile("jitter_donut.cpp", "-O0 -g")
dyn.load(dynlib("jitter_donut"))

## iv. Does not account for jittering ----
compile("no_jitter.cpp", "-O0 -g")
dyn.load(dynlib("no_jitter"))

# 5. Fit models ----

# Starting values in TMB
tmb_params <- list(beta = rep(0, length(cov.raster.list) + 1),
                   log_tau = 0,
                   log_kappa = 0,
                   log_nug_std = 0,
                   Epsilon_s = rep(0, mesh[['n']]))

# Random effects
rand_effs <-  c('Epsilon_s', 'beta')

# Fit models loop
result.list <- list()
time.fit.start <- Sys.time()
for(i in 1:num.datasets){
  temp.result <- list()
  
  print(paste("======= index i: ", i))
  for(j in 1:nrow(jittering.schemes)){
    print(paste("======= index j: ", j))
    
    scheme <- jittering.schemes[j, ]
    time.start.fit.curr <- Sys.time()
    result <- fit.sample.predict(nLoc = sample.size,
                                     data1 = inputs[[i]][[j]]$data_standard,
                                     data2 = inputs[[i]][[j]]$data_jittAccounted,
                                     parameters = tmb_params,
                                     random = rand_effs,
                                     flag2 = 0,
                                     pred.locs = pred.locs,
                                     data.test = pred.data.list[[i]],
                                     mesh.s = mesh,
                                     raster.list = cov.raster.list,
                                     true.params = true.params,
                                     jittering.scheme = scheme[1],
                                     return.Q = FALSE)
    time.end.fit.curr <- Sys.time()
    elaps.fit.curr <- time.end.fit.curr - time.start.fit.curr
    result$times <- list(start = time.start.fit.curr, end = time.end.fit.curr, elapsed = elaps.fit.curr)
    save(result, i, j, file = "temp_results.RData")
    temp.result[[j]] <- result
  }
  result.list[[i]] <- temp.result
}
time.end <- Sys.time()
elaps.fit.end <- time.end - time.fit.start
elaps.total <- time.end - time.start

result.list$times <- list(start = time.start, end = time.end, 
                          prep.time = prep.end - prep.start, fit.time = elaps.fit.end, total.time = elaps.total)

# Save results
saveRDS(result.list, file = "result.RData")












