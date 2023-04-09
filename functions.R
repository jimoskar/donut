# Functions used in simulation


# Convert from lat,lon to UTM37
convert.deg.to.km = function(loc, country = "nigeria"){
  if(country == "kenya"){
    locLatLon = SpatialPoints(loc, proj4string = CRS("+proj=longlat +datum=WGS84"))
    locKM = spTransform(locLatLon,
                        CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"))
  } else{
    locLatLon = SpatialPoints(loc, proj4string = CRS("+proj=longlat +datum=WGS84"))
    locKM = spTransform(locLatLon, CRS("+proj=tmerc +lat_0=4 +lon_0=12.5 +k=0.99975 +x_0=1110369.7 +y_0=0 +a=6378249.145 +rf=293.465 +towgs84=-92,-93,122,0,0,0,0 +units=km +no_defs +type=crs +datum=WGS84"))
  
  }                      
  return(locKM@coords[,c(1,2)])
}


convert.km.to.deg = function(loc, country = "nigeria") {
  if(country == "kenya"){
    locSP = SpatialPoints(loc, proj4string=CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"))
    lonLatCoords = spTransform(locSP, CRS("+proj=longlat +datum=WGS84"))
  } else{
    locSP = SpatialPoints(loc, proj4string=CRS("+proj=tmerc +lat_0=4 +lon_0=12.5 +k=0.99975 +x_0=1110369.7 +y_0=0 +a=6378249.145 +rf=293.465 +towgs84=-92,-93,122,0,0,0,0 +units=km +no_defs +type=crs"))
    lonLatCoords = spTransform(locSP, CRS("+proj=longlat +datum=WGS84"))
  }
  return(attr(lonLatCoords, "coords"))
}

# Takes in coordinats in km and converts to lat-lon
convert_km_to_deg <- function(locs){
  sp.pts.m <- SpatialPoints(locs*1000, CRS("EPSG:26392"))
  sp.pts.ll <- spTransform(sp.pts.m, CRS("EPSG:4326"))
  return(sp.pts.ll@coords)
}

convert_deg_to_km <- function(locs){
  sp.pts.ll <- SpatialPoints(locs, CRS("EPSG:4326"))
  sp.pts.m <- spTransform(sp.pts.ll, CRS("EPSG:26392"))
  return(sp.pts.m@coords / 1000)
}

################################################################################

# Read contraception
get.contraception = function(dat){
  # Convert the answers given to dat$v302a from yes/no into 1s and 0s :
  yesNo = rep(NA, length(dat[,1]))
  for (i in 1:nrow(dat)){
    if (dat$v302a[i] == 0){
      yesNo[i] = 0
    } else {
      yesNo[i] = 1
    }
  }
  
  # Associate with clusters
  contraceptUse = data.frame(yesNo = yesNo,
                             clustID = dat$v001)
  
  #sum 1s with respect to cluster IDs
  answers_x = aggregate(contraceptUse$yesNo,
                        by = list(clusterID = contraceptUse$clustID),
                        FUN = sum)
  answers_n= aggregate(contraceptUse$yesNo,
                       by = list(clusterID = contraceptUse$clustID),
                       FUN = length)
  answers_joint = merge(answers_x, answers_n, by="clusterID")
  colnames(answers_joint) = c("clusterID", "ys", "ns")
  
  return(answers_joint)
}


# Simulation  from SPDE----
# Function for simulating from SPDE representation
rspde<- function (coords, kappa, variance = 1, alpha = 2, n = 1, mesh,
                  verbose = FALSE, seed, return.attributes = FALSE, matern.priors = NULL)
{
  
  if(is.null(seed)) seed = 0 ## If ‘seed=0L’ then GMRFLib will set the seed
  ## intelligently/at 'random'
  
  t0 <- Sys.time()
  theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
  if (verbose)
    cat("theta =", theta, "\n")
  if (missing(mesh)) {
    mesh.pars <- c(0.5, 1, 0.1, 0.5, 1) * sqrt(alpha - ncol(coords)/2)/kappa
    if (verbose)
      cat("mesh.pars =", mesh.pars, "\n")
    attributes <- list(mesh = inla.mesh.2d(, coords[chull(coords),
    ],
    max.edge = mesh.pars[1:2],
    cutoff = mesh.pars[3],
    offset = mesh.pars[4:5]))
    if (verbose)
      cat("n.mesh =", attributes$mesh$n, "\n")
  }
  else attributes <- list(mesh = mesh)
  if(is.null(matern.priors)){
    attributes$spde <- inla.spde2.matern(attributes$mesh, alpha = alpha)
  } else {
    attributes$spde <- inla.spde2.pcmatern(mesh = mesh, alpha = alpha,
                                           prior.range = matern.priors[1:2],
                                           prior.sigma = matern.priors[3:4])
  }
  attributes$Q <- inla.spde2.precision(attributes$spde, theta = theta)
  attributes$A <- inla.mesh.project(mesh = attributes$mesh,
                                    loc = coords)$A
  if (n == 1)
    result <- drop(attributes$A %*% inla.qsample(Q = attributes$Q,
                                                 constr = attributes$spde$f$extraconstr))
  t1 <- Sys.time()
  result <- inla.qsample(n, attributes$Q, seed = ifelse(missing(seed),
                                                        0, seed),
                         constr = attributes$spde$f$extraconstr)
  if (nrow(result) < nrow(attributes$A)) {
    result <- rbind(result, matrix(NA, nrow(attributes$A) -
                                     nrow(result), ncol(result)))
    dimnames(result)[[1]] <- paste("x", 1:nrow(result), sep = "")
    for (j in 1:ncol(result)) result[, j] <- drop(attributes$A %*%
                                                    result[1:ncol(attributes$A),
                                                           j])
  }
  else {
    for (j in 1:ncol(result)) result[1:nrow(attributes$A),
                                     j] <- drop(attributes$A %*% result[, j])
    result <- result[1:nrow(attributes$A), ]
  }
  t2 <- Sys.time()
  attributes$cpu <- c(prep = t1 - t0, sample = t2 - t1, total = t2 -
                        t0)
  if (return.attributes)
    attributes(result) <- c(attributes(result), attributes)
  return(drop(result))
}

# TODO: Make this function return A and reuse them.
spde.sim <- function(loc, kappa, variance, alpha, mesh, 
                     beta.vec, cov.raster.list, nug.std, ns = 1, n = 1, seed = 0, test.idxs = NULL)
{
  set.seed(seed)
  # Simulate from SPDE
  theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
  spde <- inla.spde2.matern(mesh, alpha = alpha)
  Q <- inla.spde2.precision(spde, theta = theta)
  A <- inla.mesh.project(mesh, loc = as.matrix(loc))$A
  u.sim <- A %*% inla.qsample(Q = Q, constr = spde$f$extraconstr, seed = seed)
  u.sim <- as.vector(u.sim)
  
  # make design matrix
  X = make.design.mat(cov.raster.list, loc[,1], loc[,2])

  lin.pred = X %*% beta.vec + u.sim
  
  # Simulate Gaussian responses
  noise = rnorm(n = length(lin.pred), mean = 0, sd = nug.std)
  y.gauss = lin.pred + noise 
  
  # Simulate Binomial responses
  ps = expit(lin.pred)
  y.binom = rbinom(n = length(lin.pred),
                   size = ns,
                   prob = ps)
  
  u.sim.test <- lin.pred.test <- y.gauss.test <- y.binom.test <- NULL
  if(!is.null(test.idxs)){
    y.gauss.test <- y.gauss[test.idxs]
    y.binom.test <- y.binom[test.idxs]
    lin.pred.test <- lin.pred[test.idxs]
    u.sim.test <- u.sim[test.idxs]
    noise.test <- noise[test.idxs]
    
    y.gauss <- y.gauss[-test.idxs]
    y.binom <- y.binom[-test.idxs]
    u.sim <- u.sim[-test.idxs]
    lin.pred <- lin.pred[-test.idxs]
    noise <- noise[-test.idxs]
  }
  
  return(list(data.obs = data.frame(
                u.sim = u.sim,
                lin.pred = lin.pred,
                y.gauss = y.gauss,
                y.binom = y.binom,
                noise = noise),
              data.test = data.frame(
                u.sim.test = u.sim.test,
                lin.pred.test = lin.pred.test,
                y.gauss.test = y.gauss.test,
                y.binom.test = y.binom.test),
              Q = Q)
  )
}

sample.locations <- function(Nloc, Npred = 0, sample.df, pred.weight = FALSE){
  result <- list()
  if(Npred & pred.weight){
    N <- Nloc + Npred
    sampled.locs <- sample.df[sample.int(nrow(sample.df), N, replace = FALSE, prob = sample.df$nga_ppp_2018), c(2,3)]
    result$pred.locs <- tail(sampled.locs, Npred) 
    result$pred.idxs <- as.numeric(rownames(result$pred.locs))
    result$obs.locs <- head(sampled.locs, Nloc) 
    result$obs.idxs <- as.numeric(rownames(result$obs.locs))
  } else if(Npred & !(pred.weight)) {
    result$obs.locs <- sample.df[sample.int(nrow(sample.df), Nloc, replace = FALSE, prob = sample.df$nga_ppp_2018), c(2,3)]
    result$obs.idxs <- as.numeric(rownames(result$obs.locs))
    sample.df.pred <- sample.df[-result$obs.idxs]
    result$pred.locs <- sample.df.pred[sample.int(nrow(sample.df), Npred, replace = FALSE), c(2,3)]
    result$pred.idxs <- as.numeric(rownames(result$pred.locs))
  } else {
    N <- Nloc
    sampled.locs <- sample.df[sample.int(nrow(sample.df), N, replace = FALSE, prob = sample.df$nga_ppp_2018), c(2,3)]
    result$obs.locs <- sampled.locs
    result$obs.idxs <- as.numeric(rownames(sampled.locs))
  }
  return(result)
}

simulate.locations <- function(N, p.urban = 0.4, n.preds = 0, uniform = FALSE){
  pc.df <- readRDS("cov_data/pop4_grid_cleaned.Rda")
  if(uniform){
    locs.idxs <- sample.int(nrow(pc.df), N, replace = FALSE)
  } else{
    locs.idxs <- sample.int(nrow(pc.df), N, replace = FALSE, prob = pc.df$nga_ppp_2018)
  }
  locs <- pc.df[locs.idxs, c(2,3)]
  locs[,c("x", "y")] <- convert.deg.to.km(locs[,c("x", "y")])
  # 1 is for urban, 0 is for rural if numbers are used.
  urban <- sample(c("U", "R"), nrow(locs), replace = TRUE, prob = c(p.urban, 1 - p.urban))
  locs$urban <- urban
  locs <- list(obs.locs = locs)
  
  if(n.preds){
    pred.set <- 1:nrow(pc.df)
    pred.set <- pred.set[! pred.set %in% locs.idxs]
    pred.locs <- pc.df[sample(pred.set, n.preds, replace = FALSE), c(2,3)]
    pred.locs <- convert.deg.to.km(pred.locs)
    locs[["pred.locs"]] <- pred.locs
  }
  
  return(locs)
}

# Simualte from Matern ----

# Matérn covariance function
covMatern = function(dMat, range, stdDev){
  Sig = inla.matern.cov(nu = 1,
                        kappa = sqrt(8*1)/range,
                        x = dMat,
                        corr = TRUE)
  Sig = stdDev^2*Sig
}

# Simulate from model
simulate.responses = function(loc, 
                             ns, 
                             beta.vec,
                             cov.raster.list, 
                             range.sp, 
                             sigma.sp, 
                             nug.std){
  
  # Spatial covariance matrix
  covMat = covMatern(dMat = as.matrix(dist(loc)),
                     range = range.sp,
                     stdDev = sigma.sp)
  
  # Add small noise
  covMat = covMat + sigma.sp^2*diag(x = 1e-6^2, 
                                       nrow = dim(covMat)[1], 
                                       ncol = dim(covMat)[2])
  
  # Simulate spatial effect
  L = t(chol(covMat))
  u.sim = L%*%rnorm(dim(L)[1])
  
  # make design matrix
  X = as.matrix(rep(1, rep(nrow(loc))))
  if(!is.null(cov.raster.list)){
    X = make.design.mat(cov.raster.list, loc[,1], loc[,2])
  }
  
  lin.pred = X %*% beta.vec + u.sim
  
  # Simulate Gaussian responses
  y.gauss = rnorm(n = length(lin.pred),
                  mean = lin.pred,
                  sd = nug.std)
  
  # Simulate Binomial responses
  ps = expit(lin.pred)
  y.binom = rbinom(n = length(lin.pred),
                   size = ns,
                   prob = ps)
  
  return(data.frame(u.sim = u.sim,
                    lin.pred = lin.pred,
                    y.gauss = y.gauss,
                    y.binom = y.binom))
}


################################################################################
# Generate random distances with respect to the selected jittering scale
# type : a vector of location types : U for urban, R for rural
# s = scaling factor (1, 3 ,5, 10)
random.distance<- function(type, s){      
  distance<- rep(0, length(type))
  for (i in 1:length(type)){
    if (type[[i]]=="U"){
      distance[[i]]=runif(1, min = 0, max = 2*s)}
    else {
      if (runif(1) < 0.01){
        distance[[i]]=runif(1, min = 0, max = 10*s)
      } else{
        distance[[i]]=runif(1, min = 0, max = 5*s)
      }
      list(rand.dist=distance)
    }}
  return(distance)
}

################################################################################
#Generate random displacement angles (between 0 - 2pi) 

#Argument: length: Number of observations.
random.angle<- function(length){
  angle<- (runif(length, min = 0, max = 2*pi))
  return(angle)
}

################################################################################
# Displace locations w.r.t distance and angle generated by two functions above
displace <- function(east, north, angle, distance){
  locx=rep(0, length(north))
  locy=rep(0, length(north))	
  for (i in 1:length(north)){
    (locy[i]=north[i]+((distance[i])*sin(angle[i])))&
      (locx[i]=east[i]+((distance[i])*cos(angle[i])))}
  results=data.frame(locx=locx, locy=locy)
  return(results)
}

################################################################################
# Bring 3 functions above together. Jitter coordinate sets by respecting/not respecting admin1 areas.
# scale-->jittering scale,  locKM--> true locations as east/north (km)
# urbanRural--> should be a vector of U/R, boundary--> TRUE if we respect admin1 boundaries

jitter.locs = function(scale, locKM, urbanRural, shapeFile, check1, boundary, adm.level = 2, min.disp = 0){
  # Jitter each true coordinate one by one and then check the administrative areas they are in now
  #If they landed into a different area then their previous one, jitter that location again
  #Stop when it is jittered and stayed in the same area as befor. Continue jittering with the next location.
  eastOriginal = locKM[, "east"]
  northOriginal = locKM[,"north"]
  nLoc = length(eastOriginal)
  jitteredCoords = list()
  for (j in 1:length(scale)){
    newLocationSet=data.frame(east = rep(NA, nLoc), north = rep(NA, nLoc))
    for (i in 1:nLoc){
      repeat{
        #arguments to be used in jittering:
        east = eastOriginal[i]; north = northOriginal[i]; angle = random.angle(1)
        distance = random.distance(type = urbanRural[i], s = scale[j])
        if(min.disp > 0){
          repeat{
            distance = random.distance(type = urbanRural[i], s = scale[j])
            if(distance < min.disp){next} else{break}
          }
        }
        #jitter location i with respect to scale j  (in east/north format)
        newPoint_eastNorth = displace(east = east, north = north, angle = angle, distance = distance)
        # If we respect admin1 boundaries, check the new point against initial point in polygon table (check1)
        if (boundary == "TRUE"){
          #convert jittered location i to a spatialPoints object
          newPoint_spatialPointsObject = SpatialPoints(cbind(newPoint_eastNorth[,1], newPoint_eastNorth[,2])*1000, proj4string = CRS("EPSG:26392"), bbox = NULL)
          # Transform it also into longitude/latitude format
          newPoint_longLat <- sp::spTransform(newPoint_spatialPointsObject, shapeFile@proj4string)
          #see which admin1 area it landed in:
          check2 <- over(newPoint_longLat, shapeFile, returnList = FALSE)
          #compare it with its previous admin1 area, keep jittering the same point until two areas match
          if(adm.level == 1){
            if((is.na(check2[,"NAME_1"][[1]]) == FALSE) & (check2[,"NAME_1"][[1]] == check1[,"NAME_1"][[i]])){
              break
            } else{next}
          } else if(adm.level == 2) {
            if((is.na(check2[,"NAME_2"][[1]]) == FALSE) & (check2[,"NAME_2"][[1]] == check1[,"NAME_2"][[i]])){
              break
            } else{next}
          } else{stop("Unknown admin level")}
        }else{break}
      } #fill in the jittered location i
      newLocationSet[[i,1]] = newPoint_eastNorth[[1,1]]
      newLocationSet[[i,2]] = newPoint_eastNorth[[1,2]]
    }
    jitteredCoords[[j]] = newLocationSet
  }
  return(jitteredCoords)
}

# Jitter locations according to the donut jittering scheme 
# R1 is the inner radius while R2 us the outer radius
jitter.donut = function(locs, R1, R2, shapeFile, check1, boundary, adm.level = 2){
  # Jitter each true coordinate one by one and then check the administrative areas they are in now
  #If they landed into a different area then their previous one, jitter that location again
  #Stop when it is jittered and stayed in the same area as befor. Continue jittering with the next location.
  eastOriginal = locs[, "east"]
  northOriginal = locs[,"north"]
  nLoc = length(eastOriginal)
  jitteredCoords = list()
  newLocationSet = data.frame(east = rep(NA, nLoc), north = rep(NA, nLoc))
  for (i in 1:nLoc){
    repeat{
      #arguments to be used in jittering:
      east = eastOriginal[i]; north = northOriginal[i]; 
      angle = runif(1, 0, 2*pi)
      radius = runif(1, R1, R2)
      
      #jitter location i with respect to scale j  (in east/north format)
      newPoint_eastNorth = displace(east = east, north = north, angle = angle, distance = radius)
      
      # If we respect admin2 boundaries, check the new point against initial point in polygon table (check1)
      if (boundary == "TRUE"){
        #convert jittered location i to a spatialPoints object
        newPoint_spatialPointsObject = SpatialPoints(cbind(newPoint_eastNorth[,1], newPoint_eastNorth[,2])*1000, proj4string = CRS("EPSG:26392"), bbox = NULL)
        # Transform it also into longitude/latitude format
        newPoint_longLat <- sp::spTransform(newPoint_spatialPointsObject, shapeFile@proj4string)
        #see which admin1 area it landed in:
        check2 <- over(newPoint_longLat, shapeFile, returnList = FALSE)
        #compare it with its previous admin1 area, keep jittering the same point until two areas match
        if(adm.level == 1){
          if((is.na(check2[,"NAME_1"][[1]]) == FALSE) & (check2[,"NAME_1"][[1]] == check1[,"NAME_1"][[i]])){
            break
          } else{next}
        } else if(adm.level == 2) {
          if((is.na(check2[,"GID_2"][[1]]) == FALSE) & (check2[,"GID_2"][[1]] == check1[,"GID_2"][[i]])){
            break
          } else{next}
        } else{stop("Unknown admin level")}
      } else{break}
    } 
    
    #fill in the jittered location i
    newLocationSet[i, ] = newPoint_eastNorth
  }
  return(newLocationSet)
}

# Prepare inputs ----

# prepare TMB input w.r.t Gaussian/Binomial response
# flag2-->0! for Gaussian/Binomial response   
# sim.data--> contains spatial field, simulated Gaussian and Binomial responses

prepare.input = function(sim.data, loc.obs, model.params, priors, jittering.scale, urban, mesh.s, cov.raster.list = NULL, flag2 = 0, 
                         admin.map=NULL, nSubAPerPoint=10, nSubRPerPoint=10, beta.random = TRUE, jittering.scheme = "rp", donut.info = NULL, testMode=FALSE, return.grid = FALSE){
  
  if(!testMode){
    # extract arguments
    matern.pri <- priors$matern.pri
    beta.pri <- priors$beta.pri
    nugget.pri <- priors$nugget.pri
    alpha <- model.params$alpha
    
    # number of observed locations
    nLoc = length(loc.obs[,1])
    
    #response variable Gaussian/Binomial
    if (flag2 == 0){
      ys = sim.data$y.gauss[1:nLoc]
      ns = rep(1, nLoc)
    }else if(flag2 == 2){
      ys = sim.data$y.gauss[1:nLoc]
      ns = sim.data$n
    } else {
      ys = sim.data$y.binom[1:nLoc]
      ns = rep(100, nLoc)
    }
    
    
    # SPDE components
    
    #jittering the points a bit just to make a mesh
    #spde = getSPDEPrior(mesh.s, U=USpatial, alpha=alphaSpatial)
    spde <- inla.spde2.matern(mesh.s, alpha = alpha)
    A.proj = inla.spde.make.A(mesh = mesh.s, loc = cbind(loc.obs[,1], loc.obs[,2]))
    
    # Make design matrix
    X.standard <- as.matrix(rep(1, nLoc))
    if(!is.null(cov.raster.list)){
      X.standard <- make.design.mat(cov.raster.list, loc.obs[,1], loc.obs[,2])
    }
    num.fe <- ncol(X.standard)
    
    # Set options
    options <-  c(1, ## if 1, use normalization trick
                  1, ## if 1, run adreport
                  0) ## if beta is random
    if(beta.random){
      options[3] = 1
    }
  
    # TMB input for standard (jittering is not accounted for) model
    data_standard = list(num_i = nrow(loc.obs),  # Total number of observations
                         num_s = mesh.s$n,
                         num_fe = num.fe, # Number of fixed effects
                         y_i   =  ys, # num. of pos. obs in the cluster
                         n_i   = ns,  # num. of exposures in the cluster
                         X_standard  = X.standard, #des.mat
                         M0    = spde[['param.inla']][['M0']], # SPDE sparse matrix
                         M1    = spde[['param.inla']][['M1']], # SPDE sparse matrix
                         M2    = spde[['param.inla']][['M2']], # SPDE sparse matrix
                         Aproj = A.proj,             # Projection matrix
                         options = options,
                         flag1 = 1,# normalization flag.
                         flag2 = flag2, #(0/1 for Gaussian/Binomial)
                         beta_pri = beta.pri, ## normal
                         matern_pri = matern.pri,
                         nugget_pri = nugget.pri)
  }
  # TMB input for the model that accounts for jittering (Model-J)
  
  # convert urban U/R into TRUE/FALSE
  for (i in 1:length(urban)){
    if (urban[i]=='U'){
      urban[i]='TRUE'
    }else{
      urban[i]='FALSE'
    }
  }
  urbanVals=as.logical(urban)
  any.rural <- TRUE
  if(sum(urbanVals) == length(urban)){any.rural = FALSE}
  
  if(jittering.scheme == "rp"){
    intPointInfo = makeAllIntegrationPoints(coords = cbind(loc.obs[,1], loc.obs[,2]), urbanVals, 
                                            numPointsUrban=1+15*4, numPointsRural=1+15*9, 
                                            scalingFactor = jittering.scale, 
                                            JInnerUrban=5, JOuterUrban=0, 
                                            JInnerRural=5, JOuterRural=5, 
                                            adminMap=admin.map, 
                                            nSubAPerPoint=nSubAPerPoint, 
                                            nSubRPerPoint=nSubRPerPoint, 
                                            testMode=testMode,
                                            adm.level = 2)
  } else if( jittering.scheme == "d"){
    intPointInfo = make.integration.donut(cbind(loc.obs[,1], loc.obs[,2]),
                                          donut.info$R1, donut.info$R2, num.rings = donut.info$num.rings, 
                                          num.pts.per.ring = donut.info$num.pts.per.ring,
                                          admin.map = admin.map)
    
    x = intPointInfo$x
    y = intPointInfo$y
    w = intPointInfo$w
    
    X.mat <- make.design.mat(cov.raster.list, x = x, y = y)
    A.proj.donut <- inla.spde.make.A(mesh = mesh.s, loc = cbind(c(x), c(y)))
    
    data_jittAccounted <- list(num_i = nrow(loc.obs), # Total number of rural observations
                               num_s = mesh.s[['n']], # num. of vertices in SPDE mesh
                               num_fe = num.fe,
                               y_i   =  ys, # num. of pos. urban obs in the cluster
                               n_i = ns,
                               n_integrationPoints = ncol(w), 
                               w = w, 
                               X_mat = X.mat,
                               M0    = spde[['param.inla']][['M0']], # SPDE sparse matrix
                               M1    = spde[['param.inla']][['M1']], # SPDE sparse matrix
                               M2    = spde[['param.inla']][['M2']], # SPDE sparse matrix
                               Aproj = A.proj.donut,           # Projection matrix (rural)
                               options = options,
                               # normalization flag.
                               flag1 = 1,
                               flag2 = flag2, #(0/1 for Gaussian/Binomial)
                               beta_pri = beta.pri, ## normal
                               matern_pri = matern.pri,
                               nugget_pri = nugget.pri)
    if(return.grid){
      return(list(data_standard = data_standard, data_jittAccounted = data_jittAccounted, x = x, y = y))
    } else{
    return(list(data_standard = data_standard, data_jittAccounted = data_jittAccounted))
    }
  } else{
    stop("Unknown jittering scheme.")
  }
  if(testMode) {
    return(intPointInfo)
  }
  
  wUrban = intPointInfo$wUrban
  wRural = intPointInfo$wRural
  n_integrationPointsUrban = ncol(wUrban)
  n_integrationPointsRural = ncol(wRural)
  
  # Make design matrices
  X.urban <- make.design.mat(cov.raster.list, x = intPointInfo$xUrban, y = intPointInfo$yUrban)
  if(any.rural){
    X.rural <- make.design.mat(cov.raster.list, x = intPointInfo$xRural, y = intPointInfo$yRural)
  } else{
    X.rural <- NULL
  }
  
  
  # Construct projection matrices, and get other relevant info for TMB
  out = makeJitterDataForTMB(intPointInfo, ys , urbanVals, ns, spdeMesh=mesh.s)
  ysUrban = out$ysUrban
  ysRural = out$ysRural
  nsUrban = out$nsUrban
  nsRural = out$nsRural
  AUrban = out$AUrban
  ARural = out$ARural
  
  # Compile inputs for TMB
  data_jittAccounted <- list(num_iUrban = length(ysUrban),  # Total number of urban observations
                             num_iRural = length(ysRural),  # Total number of rural observations
                             num_s = mesh.s[['n']], # num. of vertices in SPDE mesh
                             num_fe = num.fe,
                             y_iUrban   = ysUrban, # num. of pos. urban obs in the cluster
                             y_iRural   = ysRural, # num. of pos. rural obs in the cluster
                             n_iUrban   = nsUrban,  # num. of urban exposures in the cluster
                             n_iRural   = nsRural,  # num. of rural exposures in the cluster
                             n_integrationPointsUrban = n_integrationPointsUrban, 
                             n_integrationPointsRural = n_integrationPointsRural, 
                             wUrban = wUrban, 
                             wRural = wRural, 
                             X_urban = X.urban,
                             X_rural = X.rural,
                             M0    = spde[['param.inla']][['M0']], # SPDE sparse matrix
                             M1    = spde[['param.inla']][['M1']], # SPDE sparse matrix
                             M2    = spde[['param.inla']][['M2']], # SPDE sparse matrix
                             AprojUrban = AUrban,             # Projection matrix (urban)
                             AprojRural = ARural,             # Projection matrix (rural)
                             options = options,
                             # normalization flag.
                             flag1 = 1,
                             flag2 = flag2, #(0/1 for Gaussian/Binomial)
                             beta_pri = beta.pri, ## normal
                             matern_pri = matern.pri,
                             nugget_pri = nugget.pri)
  
  return(list(data_standard = data_standard, data_jittAccounted = data_jittAccounted,
              x = intPointInfo$xUrban, y = intPointInfo$yUrban))
}

# Model fitting and prediction ----

# Simulate from multivariate normal
rmvnorm.prec <- function(mu, chol.prec, n.sims){
  z <- matrix(rnorm(length(mu) * n.sims), ncol = n.sims)
  L <- chol.prec #Cholesky(prec, super=TRUE)
  z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  mu + z
}

# Function for evaluating the performance of Model-S and Model-J
eval.model <- function(obj, pred.locs, A.pred, raster.list, 
                       num.betas, data.test, flag2, true.params, return.Q = FALSE, priv = FALSE){
  par <- obj$env$last.par
  mu <- par[names(par) %in% c("beta", "Epsilon_s")] # pick out RE
  prec <- obj$env$spHess(par, random = TRUE)

  L <- Cholesky(prec, super = T)
  
  t.draws <- rmvnorm.prec(mu = mu , chol.prec = L, n.sims = 10000)
  
  # Summarize the draws
  parnames <- names(mu)
  w.draws  <- t.draws[parnames == 'Epsilon_s',]
  beta.draws  <- matrix(t.draws[parnames == 'beta',], nrow = num.betas)
  
  # Make predictions
  X.pred <- make.design.mat(raster.list, pred.locs[,1], pred.locs[,2])
  lin.pred <- X.pred %*% beta.draws
  preds <- as.matrix(A.pred %*% w.draws)
  preds <- preds + lin.pred
  
  # Calculate prediction scores
  truth <-  data.test$y.gauss.test
  
  # Coverage
  n.preds = length(preds[,1])
  coverage = 0
  for (m in 1:n.preds){
    qnt = quantile(preds[m,], c(0.025, 0.975))
    if ((qnt[[1]] < truth[[m]])&(qnt[[2]] > truth[[m]])){
      coverage = coverage + 1
    }
  }
  coverage = coverage/n.preds
  
  #scores when jittering is not accounted for
  log.score = logs_sample(y = truth, dat = preds)
  log.score = mean(log.score)
  
  CRPS = crps_sample(y = truth, dat = preds, method = "edf")
  CRPS = mean(CRPS)
  
  if (flag2 != 0){
    # Convert to probability scale
    preds = expit(pred_simple)
  }
  
  # Find the median and sd across draws, as well as 95% intervals
  summ.preds <- cbind(mean = (apply(preds, 1, mean)),
                        median = (apply(preds, 1, median)),
                        sd     = (apply(preds, 1, sd)),
                        lower = (apply(preds, 1, quantile, .025)),
                        upper = (apply(preds, 1, quantile, .975)))
  
  # Use mean as prediction when calc. bias and RMSE
  preds.bias = mean(summ.preds[, 1] - truth)
  preds.RMSE = sqrt(mean((summ.preds[, 1] - truth)^2)) 
  
  # Summarize beta draws
  beta.median <- apply(beta.draws, 1, median)
  beta.mean <- apply(beta.draws, 1, mean)
  beta.bias <- beta.median - true.params[4]
  beta.lower <- apply(beta.draws, 1, quantile, 0.025)
  beta.upper <- apply(beta.draws, 1, quantile, 0.975)
  summ.beta <- cbind(median = beta.median, 
                     mean = beta.mean, 
                     lower = beta.lower, 
                     upper = beta.upper, 
                     length = beta.upper - beta.lower, 
                     bias = beta.bias) 
  
  
  # Summarize parameter estimates
  par.fixed <- c(range = sqrt(8)/exp(par[c("log_kappa")]),
                          sigma2.s = 1/(4.0 * 3.14159265359 * exp(2.0*par[c("log_tau")]) * exp(2.0 * par[c("log_kappa")])),
                          nug.var = exp(par[c("log_nug_std")])^2)
  summ.fixed <- cbind(par.fixed, bias = c(par.fixed[1] - true.params[1], 
                                          par.fixed[2] - true.params[2], 
                                          par.fixed[3] - true.params[3])) 
  rownames(summ.fixed) <- c("range", "sigma2.sp", "nug.var")
    
  
  if(return.Q){
    return(list(coverage = coverage,
                log.score = log.score,
                CRPS = CRPS,
                summ.preds = summ.preds,
                preds.bias = preds.bias,
                preds.RMSE = preds.RMSE,
                summ.fixed = summ.fixed,
                summ.beta = summ.beta,
                Q = prec,
                mu = mu))
  } else if(priv){ # Return things for calculation of privacy protection
    return(list(coverage = coverage,
                log.score = log.score,
                CRPS = CRPS,
                summ.preds = summ.preds,
                preds.bias = preds.bias,
                preds.RMSE = preds.RMSE,
                summ.fixed = summ.fixed,
                summ.beta = summ.beta))
  } else{
    return(list(coverage = coverage,
              log.score = log.score,
              CRPS = CRPS,
              summ.preds = summ.preds,
              preds.bias = preds.bias,
              preds.RMSE = preds.RMSE,
              summ.fixed = summ.fixed,
              summ.beta = summ.beta))
  }
}


# Fit Model-S and Model-J with TMB based on input data. Returns metrics and summaries of parameter estimates.
fit.sample.predict <- function(nLoc, data1, data2, parameters, random, flag2, 
                               pred.locs, data.test, raster.list, mesh.s, true.params, jittering.scheme = "rp", return.Q = FALSE, priv = FALSE, obs.locs = NULL){
  
  overall.sim.start = Sys.time()
  
  map=list(log_nug_std= as.factor(NA))
  
  fitStandard_start = Sys.time()
  # 1.MODEL-S
  # make the autodiff generated likelihood func & gradient
  if (flag2 != 0){
    obj.standard <- MakeADFun(data = data1,
                           parameters = parameters,
                           map = map,
                           random = random,
                           hessian = TRUE,
                           DLL = 'no_jitter')
  } else {
    obj.standard <- MakeADFun(data = data1,
                           parameters = parameters,
                           random = random,
                           hessian = TRUE,
                           DLL = 'no_jitter')
  }
  
  # We can normalize the GMRF outside of the nested optimization,
  # avoiding unnecessary and expensive Cholesky operations.
  obj.standard <- normalize(obj.standard, flag = "flag1", value = 0)
  
  #newtonOption(obj.standard, smartsearch=TRUE)
  
  print("====================== Fitting standard! =======================")
  opt0 <- nlminb(start  = obj.standard[['par']],
                 objective = obj.standard[['fn']],
                 gradient = obj.standard[['gr']],
                 lower = rep(-10, length(obj.standard[['par']])),
                 upper = rep( 10, length(obj.standard[['par']])),
                 control = list(trace = 0))   
 
  A.pred = inla.spde.make.A(mesh = mesh.s, loc = as.matrix(pred.locs)) # TODO: could take this as input from earlier.
  
  eval.standard <- eval.model(obj.standard, pred.locs, A.pred, raster.list, data1$num_fe, data.test, flag2, 
                              true.params, return.Q)
  
  # 2.MODEL-J ----
  fitJittAcc_start = Sys.time()
  
  if (flag2 != 0){
    paramsNew <- list(beta = unlist(opt0[["par"]])[names(unlist(opt0[["par"]])) == "beta"], # FE parameters
                      log_tau = opt0[["par"]][["log_tau"]], # Log tau (i.e. log spatial precision, Epsilon)
                      log_kappa = opt0[["par"]][["log_kappa"]], # SPDE parameter related to the range
                      Epsilon_s = rep(0, mesh.s[['n']]), # RE on mesh vertices
                      log_nug_std = log(sqrt(0.1))
    )
    
  } else {
    # What is this good for?
    paramsNew <- list(beta = obj.standard$env$last.par[names(obj.standard$env$last.par) == "beta"], #rep(0, data2$num_fe), # betas
                      log_tau = opt0[["par"]][["log_tau"]], # Log tau (i.e. log spatial precision, Epsilon)
                      log_kappa = opt0[["par"]][["log_kappa"]], # SPDE parameter related to the range
                      Epsilon_s = obj.standard$env$last.par[names(obj.standard$env$last.par) == "Epsilon_s"], #rep(0, mesh.s[['n']]), # RE on mesh vertices
                      log_nug_std = opt0[["par"]][["log_nug_std"]])
    
  }
  
  ## Make the autodiff generated likelihood func & gradient
  DLL = "jitter"
  if(jittering.scheme == "d"){
    DLL = "jitter_donut"
  } else if(data2$num_iRural == 0){
    DLL = "no_rural"
    data2$X_rural <-  NULL; data2$y_iRural <- NULL; data2$n_iRural <-  NULL; data2$n_integrationPointsRural <- NULL; 
    data2$n_iRural <- NULL; data2$wRural <- NULL; data2$AprojRural <- NULL; data2$num_iRural <- NULL;
  }

  if (flag2 != 0){
    obj.jitter <- MakeADFun(data = data2,
                     parameters = paramsNew,
                     map = map,
                     random = random,
                     hessian = TRUE,
                     DLL = DLL)
  } else {
    obj.jitter <- MakeADFun(data = data2,
                     parameters = paramsNew,
                     random = random,
                     hessian = TRUE,
                     DLL = DLL)
  }
  
  ## We can normalize the GMRF outside of the nested optimization,
  ## avoiding unnecessary and expensive cholesky operations.
  obj.jitter <- normalize(obj.jitter, flag="flag1", value = 0)
  
  # * Run TMB ----
  print("====================== Fitting jitter! =======================")
  opt0 <- nlminb(start = obj.jitter[['par']],
                 objective = obj.jitter[['fn']],
                 gradient = obj.jitter[['gr']],
                 lower = rep(-10, length(obj.jitter[['par']])),
                 upper = rep( 10, length(obj.jitter[['par']])),
                 control = list(trace=1))
  
  
  # Evalualte model
  eval.jitter <- eval.model(obj.jitter, pred.locs, A.pred, raster.list, data2$num_fe, data.test, flag2, 
                            true.params, return.Q, priv = TRUE)
  
  pp.result <-  NULL
  if(priv){
    pp.result <- calculate.privacy.protection(list(eval.jitter = eval.jitter), obs.locs, R2)
  }
  
  overall.sim.end <- Sys.time()
  
  return(list(overall.sim.start = overall.sim.start,
              overall.sim.end = overall.sim.end,
              eval.standard = eval.standard,
              eval.jitter = eval.jitter,
              pp.result = pp.result))
              #obj.s = obj.standard,
              #obj.j = obj.jitter))
}

fit.model <- function(n.loc, data1, data2, parameters, random, flag2, mesh.s, start.values = TRUE){
  map=list(log_nug_std= as.factor(NA))
  ## 1.MODEL-S ----
  # make the autodiff generated likelihood func & gradient
  if (flag2 != 0 && flag2 != 2){
    obj.standard <- MakeADFun(data = data1,
                              parameters = parameters,
                              map = map,
                              random = random,
                              hessian = TRUE,
                              DLL = 'no_jitter')
  } else {
    obj.standard <- MakeADFun(data = data1,
                              parameters = parameters,
                              random = random,
                              hessian = TRUE,
                              DLL = 'no_jitter')
  }
  
  # We can normalize the GMRF outside of the nested optimization,
  # avoiding unnecessary and expensive Cholesky operations.
  obj.standard <- normalize(obj.standard, flag = "flag1", value = 0)
  
  #newtonOption(obj.standard, smartsearch=TRUE)
  
  print("====================== Fitting standard! =======================")
  opt0 <- nlminb(start  = obj.standard[['par']],
                 objective = obj.standard[['fn']],
                 gradient = obj.standard[['gr']],
                 lower = rep(-10, length(obj.standard[['par']])),
                 upper = rep( 10, length(obj.standard[['par']])),
                 control = list(trace = 0))   
  
  ## 2.MODEL-J ----
  if(start.values){
    paramsNew <- list(beta = obj.standard$env$last.par[names(obj.standard$env$last.par) == "beta"], #rep(0, data2$num_fe), # betas
                      log_tau = opt0[["par"]][["log_tau"]], # Log tau (i.e. log spatial precision, Epsilon)
                      log_kappa = opt0[["par"]][["log_kappa"]], # SPDE parameter related to the range
                      Epsilon_s = obj.standard$env$last.par[names(obj.standard$env$last.par) == "Epsilon_s"], #rep(0, mesh.s[['n']]), # RE on mesh vertices
                      log_nug_std = opt0[["par"]][["log_nug_std"]])
  } else {
    paramsNew <- list(beta = rep(0, data2$num_fe), # betas
                      log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                      log_kappa = 0, # SPDE parameter related to the range
                      Epsilon_s = rep(0, mesh.s[['n']]), # RE on mesh vertices
                      log_nug_std = 0)
  }
  
  ## Make the autodiff generated likelihood func & gradient
  if (flag2 != 0 && flag2 != 2){
    obj.jitter <- MakeADFun(data = data2,
                            parameters = paramsNew,
                            map = map,
                            random = random,
                            hessian = TRUE,
                            DLL = 'jitter')
  } else {
    obj.jitter <- MakeADFun(data = data2,
                            parameters = paramsNew,
                            random = random,
                            hessian = TRUE,
                            DLL = 'jitter')
  }
  
  ## We can normalize the GMRF outside of the nested optimization,
  ## avoiding unnecessary and expensive cholesky operations.
  obj.jitter <- normalize(obj.jitter, flag="flag1", value = 0)
  
  # Run TMB 
  print("====================== Fitting jitter! =======================")
  opt0 <- nlminb(start = obj.jitter[['par']],
                 objective = obj.jitter[['fn']],
                 gradient = obj.jitter[['gr']],
                 lower = rep(-10, length(obj.jitter[['par']])),
                 upper = rep( 10, length(obj.jitter[['par']])),
                 control = list(trace=1))
  
  return(list(obj.j = obj.jitter, obj.s = obj.standard))
}

## Covariates ----
# make.cov.array <- function(raster.list, x.mat, y.mat){
#   l <- length(raster.list); n <- nrow(x.mat); m <- ncol(x.mat);
#   x <- as.vector(x.mat)
#   y <- as.vector(y.mat)
#   arr <- array(dim = c(l, n, m))
#   for(i in 1:l){
#     coords <- convert.km.to.deg(cbind(x, y))
#     vals <- extract(raster.list[[i]], coords)
#     val.mat <- matrix(vals, nrow = n, ncol = m)
#     arr[i, , ] = val.mat
#   }
#   return(arr)
# }

# Min-max scale a raster
min.max <- function(r){
  r.min <- cellStats(r, "min")
  r.max <- cellStats(r, "max")
  r.scale <- (r - r.min)/(r.max - r.min)
  return(r.scale)
}

make.design.mat <- function(raster.list, x.mat, y.mat){
  l <- length(raster.list); n <- nrow(x.mat); m <- ncol(x.mat);
  # if a vector is given instead of a matrix
  if(is.null(m)){
    m <- 1; n <- length(x.mat);
  }
  # if only intercept is used
  if(is.null(raster.list)){
    l <- 0
  }
  x <- as.vector(x.mat)
  y <- as.vector(y.mat)
  design.mat <- matrix(nrow = n*m, ncol = l + 1)
  design.mat[, 1] <- rep(1, n*m)
  if(!is.null(raster.list)){
      coords <- convert_km_to_deg(cbind(x, y))
    for(i in 1:l){
      vals <- terra::extract(raster.list[[i]], coords)
      vals <- array(unlist(vals))
      if(sum(is.na(vals))){
        vals[which(is.na(vals))] = 0
      }
      design.mat[, i + 1] = vals
    }
  }
  return(design.mat)
}

# This takes a dataframe with both x and y, not matrices x.mat and y.mat.
make.design.mat2 <- function(raster.list, locs){
  l <- length(raster.list); n <- nrow(locs);
  # if only intercept is used
  if(is.null(raster.list)){
    l <- 0
  }
  design.mat <- matrix(nrow = n, ncol = l + 1)
  design.mat[, 1] <- rep(1, n)
  if(!is.null(raster.list)){
    locs.lon.lat <- convert_km_to_deg(locs)
    for(i in 1:l){
      vals <- terra::extract(raster.list[[i]], locs.lon.lat)
      vals <- array(unlist(vals))
      if(sum(is.na(vals))){
        #browser()
        vals[which(is.na(vals))] = 0
      }
      design.mat[, i + 1] = vals
    }
  }
  return(design.mat)
}

make.design.mat3 <- function(spatrast, locs){
  l <- 0; n <- nrow(locs);
  # if only intercept is used
  if(!is.null(spatrast)){
    l <- spatrast@ptr$nlyr()
  }
  design.mat <- matrix(nrow = n, ncol = l + 1)
  design.mat[, 1] <- rep(1, n)
  if(!is.null(spatrast)){
    locs.lon.lat <- convert_km_to_deg(locs)
    vals <- as.matrix(extract(spatrast, locs.lon.lat))
    vals[is.nan(vals)] = 0
    design.mat[, 2:(l + 1)] = vals
  }
  return(design.mat)
}

get.adm.names <- function(coords, admin.map, admin.level = 2, non.spatial = FALSE){
  if(non.spatial){
    coords = SpatialPoints(coords, proj4string=admin.map@proj4string, bbox = NULL)
  }
  out = over(coords, admin.map, returnList = FALSE)
  if(admin.level == 1){
    return(out$NAME_1)
  }else{
    return(out$NAME_2)
  }
}

# Convert raster to dataframe 
raster.to.df <- function(raster){
  raster.pts <- rasterToPoints(raster, spatial = TRUE)
  raster.df <- data.frame(raster.pts)
  rm(raster.pts)
  return(raster.df)
}


tau.kappa.to.sp <- function(tau, kappa, log = TRUE){
  if(log){
    kap = exp(kappa)
    range = sqrt(8)/kap
    sigma = 1/(sqrt(4*pi)*kap*exp(tau))
    return(data.frame(range = range, sigma = sigma))
  } else{
    range = sqrt(8)/kappa
    sigma = 1/(sqrt(4*pi)*kappa*tau) 
    return(data.frame(range = range, sigma = sigma))
  }
}

sim.and.jitter.data <- function(obs.data.list, obs.locs.list, pred.data.list, true.locs.list,
                                prior.list){
  for(i in 1:num.datasets){
    obs.locs.temp <- list()
    inputs.temp <- list()
    
    # Sample EAs uniformly
    EA.idxs <- sample.int(nrow(master.frame), sample.size, replace = FALSE)
    true.locs.lat.lon <- master.frame[EA.idxs, c("x", "y")]
    true.locs.km <- convert.deg.to.km(true.locs.lat.lon)
    colnames(true.locs.km)[1:2] <- c("east", "north") 
    true.locs.list[[i]] <- true.locs.km
    
    # Simulate responses
    sim.locs <-rbind(true.locs.list[[i]], pred.locs)
    pred.idxs <- nrow(rbind(true.locs.list[[i]]) + 1):(nrow(rbind(true.locs.list[[i]])) + nrow(pred.locs))
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
      if(scheme[1] == "d"){donut.info <- list(num.pts.per.ring = 10,
                                              num.rings = 5,
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
  return(list(obs.data.list = obs.data.list, 
              obs.locs.list = obs.locs.list, 
              pred.data.list = pred.data.list,
              true.locs.list = true.locs.list))
}

geomask.pdf <- function(s.true, s.obs, R2){
  d <- sqrt((s.true[1] - s.obs[1])^2 + (s.true[2] - s.obs[2])^2)
  return(1/(2*pi*R2*d))
}

calc.pdf <- function(A, X, y, z.samples, nug.std, log = FALSE){
  mc.samples <- ncol(z.samples)
  num.betas <- ncol(X)
  pdfs <- matrix(NA, nrow = nrow(X), ncol = mc.samples)
  for(j in 1:mc.samples){
    beta <- z.samples[1:num.betas, j]; w <- z.samples[(num.betas + 1):nrow(z.samples), j]
    mu <- X%*%beta + A%*%w
    pdfs[,j] = dnorm(y, mean = as.numeric(mu), sd = nug.std)
  }
  if(log){return(log(pdfs))}
  return(pdfs)
}

calc.integrals <- function(z.samples, a, x, weights, y, nug.std, log = FALSE){
  K <- ncol(z.samples)
  num.betas = length(x)
  integrals <- rep(NA, K)
  for(k in 1:K){
    beta <- z.samples[1:num.betas, k]; w <- z.samples[(num.betas + 1):nrow(z.samples), k]
    fe <- x %*% beta; proj.w <- a %*% w
    tmp.pdf <- dnorm(y, fe + proj.w, nug.std)
    integrals[k] = robust.mix(tmp.pdf, weights)
  }
  if(log){return(log(integrals))}
  return(integrals)
}

robust.mix <- function(tmp.pdf, weights){
  log.pdf <- log(tmp.pdf)
  max.val <- max(log.pdf + log(weights))
  tmp.val = 0.0
  for(k in 1:length(tmp.pdf)){
    tmp.val = tmp.val + exp(log.pdf[k] +log(weights[k]) - max.val)
  }
  tmp.val = exp(max.val) + tmp.val
  return(tmp.val)
}

calculate.privacy.protection <- function(tmb.res, obs.locs, data.jitter, cov.raster.list, mesh, true.locs, scheme = c("rp", 0, 2), grid.res = 30){
  N <- nrow(obs.locs)
  R1 <- as.numeric(scheme[2])
  R2 <- as.numeric(scheme[3])
  dist.pri <- dist.post <-  rep(NA, N)
  areas.pri <- areas.post <-  rep(NA, N)
  w.vec <- rep(NA, N)
  post.df.list <- list()
  exp.locs <- matrix(NA, nrow = N, ncol = 2)
  mc.samples <- 1000 # Number of monte carlo samples
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = N,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)
  print("Calculating privacy protection!")
  for(idx in 1:N){
    nug.std <- sqrt(tmb.res$eval.jitter$summ.fixed[3, 1])
    Q.z <- tmb.res$eval.jitter$Q
    L.z <- Cholesky(Q.z, super = T)
    mu.z <- tmb.res$eval.jitter$mu
    z.samples <- rmvnorm.prec(as.vector(mu.z), L.z, mc.samples)

    # Make grid
    if(scheme[1] == "rp"){
      center.coord <- data.frame("east" = obs.locs[idx, 1], "north" = obs.locs[idx, 2])
      x.range <- seq(from = center.coord$east - R2, to = center.coord$east + R2, length.out = grid.res)
      y.range <- seq(from = center.coord$north - R2, to = center.coord$north + R2, length.out = grid.res)
      wx <- 2*R2/(grid.res - 1) # width of cells in x dir
      wy <- 2*R2/(grid.res - 1) # width of cells in y direction
      w.vec[idx] = wx
      grid <- expand.grid(x.range, y.range); colnames(grid) <- c("east", "north")
      cc.sp <- st_as_sf(center.coord * 1000, coords = c("east", "north"), crs = st_crs(26392))
      pts.sp <- st_as_sf(grid * 1000, coords = c("east", "north"), crs = st_crs(26392))
      circle.sf <- st_buffer(cc.sp, dist = R2*1000)
      pts.int <- st_intersection(pts.sp, circle.sf)
      pts.km <- st_coordinates(pts.int)/1000
    } else if(scheme[1] == "d"){
      center.coord <- data.frame("east" = obs.locs[idx, 1], "north" = obs.locs[idx, 2])
      x.range <- seq(from = center.coord$east - R2, to = center.coord$east + R2, length.out = grid.res)
      y.range <- seq(from = center.coord$north - R2, to = center.coord$north + R2, length.out = grid.res)
      wx <- 2*R2/(grid.res - 1) # width of cells in x dir
      wy <- 2*R2/(grid.res - 1) # width of cells in y direction
      w.vec[idx] = wx
      grid <- expand.grid(x.range, y.range); colnames(grid) <- c("east", "north")
      cc.sp <- st_as_sf(center.coord * 1000, coords = c("east", "north"), crs = st_crs(26392))
      pts.sp <- st_as_sf(grid * 1000, coords = c("east", "north"), crs = st_crs(26392))
      circle.sf <- st_buffer(cc.sp, dist = R2*1000)
      circle.outer.sf <- st_buffer(cc.sp, dist = R2*1000)
      circle.inner.sf <- st_buffer(cc.sp, dist = R1*1000)
      donut.sf <- st_difference(circle.outer.sf, circle.inner.sf)
      pts.int <- st_intersection(pts.sp, donut.sf)
      pts.km <- st_coordinates(pts.int)/1000
    }
    
    # Calculate integral in MC-integration sum
    if(scheme[1] == "rp"){
      weights <- data.jitter$wUrban[idx, ]
      a <- data.jitter$AprojUrban[idx, ]
      x <- data.jitter$X_urban[idx, ]
      y <- data.jitter$y_iUrban[idx]
    } else if(scheme[1] == "d"){
      weights <- data.jitter$w[idx, ]
      a <- data.jitter$Aproj[idx, ]
      x <- data.jitter$X_mat[idx, ]
      y <- data.jitter$y_i[idx]
    }
    ints <- calc.integrals(z.samples, a, x, weights, y, nug.std, log = FALSE)
    
    # Calculate numerator in MC-integration
    X <- make.design.mat2(cov.raster.list, locs = pts.km)
    A <- inla.spde.make.A(mesh, pts.km)
    pdf.y <- calc.pdf(A, X, y,  z.samples, nug.std, log = FALSE)
    
    # sum it up
    mc.ints <- sweep(pdf.y, 2, ints, FUN = "/")
    mc.ints <- apply(mc.ints, 1, mean)
    
    
    # Add contributions from locations
    posterior <- rep(NA, nrow(pts.km))
    prior <- rep(NA, nrow(pts.km))
    for(i in 1:nrow(pts.km)){
      pri <- geomask.pdf(pts.km[i,], as.numeric(center.coord), R2 = R2)
      prior[i] = pri
      posterior[i] = pri * mc.ints[i]
    }
   
    area <- sum(posterior*wx*wy)
    areas.post[idx] = area
    area.pri <- sum(prior*wx*wy)
    areas.pri[idx] = area.pri
    prior.scaled <- prior/area.pri
    post.scaled <- posterior/area
    
    if(any(is.na(posterior))){
      browser()
      post.df <- data.frame(post = posterior,
                            pri = prior,
                            x = pts.km[, 1],
                            y = pts.km[, 2],
                            mc.ints = mc.ints,
                            na.ints = any(is.na(ints))
      )
    } else{
      post.df <- data.frame(post = posterior,
                            pri = prior,
                            x = pts.km[, 1],
                            y = pts.km[, 2]
                            )
    }
    post.df.list[[idx]] = post.df
    expected.loc <- apply(post.scaled*wx*wy*pts.km, 2, sum)
    exp.locs[idx, ] = expected.loc
    true.loc <- true.locs[idx, ]
    dist.post[idx] <- as.numeric(sqrt((expected.loc[1] - true.loc[1])^2 + (expected.loc[2] - true.loc[2] )^2))
    dist.pri[idx] <- as.numeric(sqrt((center.coord[1] - true.loc[1])^2 + (center.coord[2] - true.loc[2] )^2))
    pb$tick()
  }
  return(list(dfs = post.df.list, dist.post = dist.post, dist.pri = dist.pri, exp.locs = exp.locs,
              areas.post = areas.post, areas.pri = areas.pri, w.vec = w.vec))
}


