# Integration scheme for donut geomasking ----


# Get integration points and weights
# R1 is inner radius, while R2 is outer radius
# num.rings: number of rings in integration donut
# num.pts.per.ring: number of points per ring in donut
get.pts.and.weights <- function(coords, 
                                R1, R2, 
                                num.rings, 
                                num.pts.per.ring,
                                verbose = TRUE){
  ms <- rep(num.pts.per.ring, num.rings) # Number of points in each ring
  r1.vec <- rev(cumsum(ms))
  r2.vec <- c(0, head(cumsum(ms), -1))
  rs <- (R1*r1.vec + R2*r2.vec)/sum(ms)
  rs <- c(tail(rs, -1), R2) # radii of rings (outer)
  
  # Areas for each ring segment
  As <- pi*rs^2
  As <- As - c(pi*R1^2, As[1:(num.rings-1)])
  As/ms
  
  # Helper function for calculating values of pi(s_i | s_ijk^*) ?
  densityFun = function(d, R1 = R1, R2 = R2) {
    out = 1 / (R2 * 2 * pi * d)
    out[d > R2] = 0
    out[d < R1] = 0
    out[d < 0] = NaN
    return(out)
  }
  
  rsIntegrationPointsMidpoint = c(R1, rs[-num.rings] + diff(rs) / 2) # midpoint solution
  
  # Calculate int point position
  aDiff = 2 * pi / ms
  shrinkFactor = sqrt(2 * (1 - cos(aDiff))) / aDiff # a scaling factor smaller than one
  shrinkFactor[1] = 1
  rsIntegrationPoints = rsIntegrationPointsMidpoint * shrinkFactor
  tooShrunk = rsIntegrationPoints < c(R1, rs[-num.rings])
  if(any(tooShrunk)) {
    warning(paste0("Center of mass is outside integration area for rings ", 
                   paste(tooShrunk, collapse=", "), 
                   ". Setting integration point to closest within the integration area"))
    rsIntegrationPoints[toShrunk] = c(0, rs[-num.rings])[tooShrunk]
  }
  
  # Weights
  ws = rep(1/sum(ms), num.rings) #(rs - c(R1, rs[-num.rings])) / rs[num.rings] * 1 / ms
  if(verbose){
    print(paste("This should sum to one:", sum(ms*ws)))
    #print(paste("This should be 1 (the sum of the weights), and should all be equal:", sum(ms*ws)))
    #print(matrix(ws, nrow=1))
  }
  
  # Collect angles and int. pts. in lists
  pts = list()
  as = list()
  for(j in 1:num.rings) {
    thisas = seq(0, 2*pi, l=ms[j]+1)[-(ms[j]+1)]
    if(j %% 2 == 1) {
      thisas = thisas + pi / ms[j]
    }
    as = c(as, as=list(thisas))
    thisPts = rsIntegrationPoints[j] * cbind(cos(thisas), sin(thisas))
    pts = c(pts, pts=list(thisPts))
  }
  
  names(as) = paste0("as", 1:num.rings)
  names(pts) = paste0("pts", 1:num.rings)
  
  return(list(rs=rs, ms=ms, As=As, ws=ws, 
              pts=pts, as=as, ptRs=rsIntegrationPoints, 
              densityFun=densityFun))
}

# Make the integration donut
# R1 is inner radius, while R2 is outer radius
# num.rings: number of rings in integration donut
# num.pts.per.ring: number of points per ring in donut
make.integration.donut <- function(coords, 
                                    R1, R2, 
                                    num.rings, 
                                    num.pts.per.ring,
                                    admin.map,
                                    admin.level = 2,
                                    verbose = TRUE){
  # Get weights and points
  out <- get.pts.and.weights(coords, R1, R2, num.rings, num.pts.per.ring)
  
  xs.vec <- unlist(lapply(out$pts, function(x){x[,1]}))
  ys.vec <- unlist(lapply(out$pts, function(x){x[,2]}))
  x.mat <- outer(coords[,1], xs.vec, "+")
  y.mat <- outer(coords[,2], ys.vec, "+")
  
  ws.vec <- rep(out$ws, out$ms)
  ws.mat <- matrix(ws.vec, ncol = length(ws.vec), nrow = nrow(coords))
  
  if(!is.null(admin.map)){
    coordsLonLat = convert_km_to_deg(coords)
    spCoordsLonLat = SpatialPoints(coordsLonLat, admin.map@proj4string)
    temp = over(spCoordsLonLat, admin.map, returnList=FALSE)
    if(admin.level == 2){
      adminID = temp$GID_2
    } else{
      adminID = temp$NAME_1
    }
    # Make dictionary so that area names can be used as OBJECTID
    if(admin.level == 2){
      counties <- admin.map$GID_2
    } else{
      counties <- admin.map$NAME_1
    }
    dict <- c()
    for(i in 1:length(counties)){
      dict[counties[i]] = i
    }
    
    # calculate distances to admin boundaries
    adminMapPolygons = as.SpatialPolygons.PolygonsList(admin.map@polygons, admin.map@proj4string)
    require(geosphere)
    dists = sapply(1:nrow(coords), function(ind) {dist2Line(spCoordsLonLat[ind], adminMapPolygons[dict[adminID[ind]]])[1]}) * (1/1000)
    
    #browser()
    # set whether or not to update weights based on distance to admin boundaries
    updateI = dists < R2
    
    # calculate updated weights for the integration points near the borders
    tempCoords = coords[updateI,]
    require(fields)
    
    tempNewWs <- update.weights(coords=tempCoords,
                                         admin.map=admin.map, 
                                         integrationPoints = out, R1 = R1,
                                         nSubAPerPoint=nSubAPerPoint, 
                                         nSubRPerPoint=nSubRPerPoint,
                                         admin.level = admin.level)
    
    # Update the weights with the new values
    #browser()
    ws.mat[updateI, ] = tempNewWs$ws.adj
  }
  
  return(list(x = x.mat, y = y.mat, w = ws.mat))
}

update.weights <- function(coords, admin.map, # Number of observations
                           integrationPoints, R1, 
                           nSubAPerPoint=10, nSubRPerPoint=10, 
                           testMode=FALSE, admin.level = 2) {
  
  adminMapPoly = as.SpatialPolygons.PolygonsList(admin.map@polygons, admin.map@proj4string)
  
  # calculate set of typical sub-integration points for urban and rural clusters
  subIntegrationPoints = get.sub.integration.pts(integrationPoints=integrationPoints, R1 = R1, 
                                                 nSubAPerPoint=nSubAPerPoint, 
                                                 nSubRPerPoint=nSubRPerPoint)
  
  
  # get admin areas associated with coordinates
  coordsLonLat = convert_km_to_deg(coords)
  spCooordsLonLat = SpatialPoints(coordsLonLat, proj4string=admin.map@proj4string, bbox = NULL)
  out = over(spCooordsLonLat, admin.map, returnList = FALSE)
  adminNames = out$NAME_1
  #adminIDs = out$OBJECTID
  if(admin.level == 2){
    adminIDs = out$GID_2
  } else{
    adminIDs = out$NAME_1
  }
  
  # Make dictionary so that area names can be used as OBJECTID
  if(admin.level == 2){
    counties <- admin.map$GID_2
  } else{
    counties <- admin.map$NAME_1
  }
  
  dict <- c()
  for(i in 1:length(counties)){
    dict[counties[i]] = i
  }
  
  # for each jittered coordinate:
  #   for each integration point:
  #     get associated sub-integration points
  #     get proportion in correct admin area
  #     update integration point weight
  
  ws.adj = matrix(nrow = nrow(coords), ncol = sum(sapply(integrationPoints$pts, function(x) {nrow(x)})))
  for(i in 1:nrow(coords)){
    # time1 = proc.time()[3]
    theseCoords = matrix(coords[i,], nrow=1)
    thisArea = adminNames[i]
    thisAreaID = dict[adminIDs[i]] #adminIDs[i]
    thisPoly = adminMapPoly[thisAreaID]
    
    # get sub-integration points
    thisSubWs = subIntegrationPoints$subWs
    thisSubPts = subIntegrationPoints$subPts
    thisSubPts = lapply(thisSubPts, function(x) {sweep(x, 2, unlist(theseCoords), "+")})
    thisSubPtsSP = lapply(thisSubPts, function(x) {
      SpatialPoints(x*1000, proj4string=CRS("EPSG:26392"))
    })
    
    # project subPts to correct projection
    thisSubPtsSPLonLat = lapply(thisSubPtsSP, function(x) {spTransform(x, admin.map@proj4string)})
    
    # determine if each sub-integration point is in correct admin area
    # the following code is commented out because it takes far too long
    
    # subAreas = lapply(thisSubPtsSPLonLat, function(x) {over(x, adminMap, returnList=FALSE)$NAME_1})
    # goodAreas = lapply(subAreas, function(x) {x == thisArea})
    goodAreas <- lapply(thisSubPtsSPLonLat, function(x) {!is.na(over(x, thisPoly, returnList=FALSE))})
    # require(fields)
    # system.time(goodAreas2 <- lapply(thisSubPtsSPLonLat, function(x) {in.poly(attr(x, "coords"), attr(thisPoly@polygons[[1]]@Polygons[[1]], "coords"))}))
    # require(ptinpoly)
    # system.time(goodAreas3 <- lapply(thisSubPtsSPLonLat, function(x) {pip2d(attr(thisPoly@polygons[[1]]@Polygons[[1]], "coords"), attr(x, "coords"))}))
    # require(Rcpp)
    # system.time(goodAreas4 <- lapply(thisSubPtsSPLonLat, function(x) {pnt.in.poly2(attr(x, "coords"), attr(thisPoly@polygons[[1]]@Polygons[[1]], "coords"))}))
    
    # update weights for sub-integration points
    updatedSubWs = thisSubWs
    updatedSubWs = lapply(1:length(updatedSubWs), function(x) {
      temp = updatedSubWs[[x]]
      temp[!goodAreas[[x]]] = 0
      temp
    })
    
    # sum sub-integration weights to get new (unnormalized) integration weights
    nSubPts = nSubRPerPoint * nSubAPerPoint
    tempWs = lapply(updatedSubWs, function(x) {
      nIntPts = length(x) / nSubPts
      aggIDs = rep(1:nIntPts, each=nSubPts)
      aggregate(x, by=list(aggIDs), FUN=sum)$x
    })
    tempWs = unlist(tempWs)
    
    # normalize new integration weights to sum to 1
    finalWs = tempWs/sum(tempWs)
    
    # update weights matrix
    ws.adj[i, ] = finalWs
    
    # time2 = proc.time()[3]
    # print(paste0("Iteration ", i, "/", nrow(coords), " took ", round(time2-time1, 2), " seconds"))
  }
  return(list(ws.adj = ws.adj, subPts=thisSubPts, goodPts=goodAreas, updatedSubWs=updatedSubWs))
}


get.sub.integration.pts <- function(integrationPoints, R1, 
                                    nSubAPerPoint=10, nSubRPerPoint=10) {
  rs = integrationPoints$rs
  ms = integrationPoints$ms
  As = integrationPoints$As
  ws = integrationPoints$ws
  pts = integrationPoints$pts
  as = integrationPoints$as
  ptRs = integrationPoints$ptRs
  densityFun = integrationPoints$densityFun
  
  #centerCoords = pts[[1]]
  
  # generates sub-integration points for any point
  getSubIntegrationPointsForOnePoint = function(minR, maxR, minA, maxA){
    widthR = (maxR - minR)/(nSubRPerPoint)
    widthA = (maxA - minA)/(nSubAPerPoint)
    
    # calculate radial coordinates of the sub-integration points
    
    # get angular coordinates of the centers of mass
    if(minA <= maxA) {
      theseAs = seq(minA + widthA/2, maxA - widthA/2, by=widthA)
    } else {
      theseAs = seq(minA + widthA/2 - 2*pi, maxA - widthA/2, by=widthA)
      theseAs[theseAs < 0] = theseAs[theseAs < 0] + 2*pi
    }
    
    # now get radial centers of mass
    theseRs = seq(minR + widthR/2, maxR - widthR/2, by=widthR)
    rsIntegrationPointsMidpoint = theseRs # midpoint solution
    
    aDiff = widthA
    shrinkFactor = sqrt(2 * (1 - cos(aDiff))) / aDiff # a scaling factor smaller than one
    shrinkFactor[1] = 1
    rsIntegrationPoints = rsIntegrationPointsMidpoint * shrinkFactor
    tooShrunk = rsIntegrationPoints < (rsIntegrationPointsMidpoint - widthR/2)
    if(any(tooShrunk)) {
      warning(paste0("Center of mass is outside integration area for rings ", 
                     paste(tooShrunk, collapse=", "), 
                     ". Setting integration point to closest within the integration area"))
      rsIntegrationPoints[toShrunk] = rsIntegrationPointsMidpoint[toShrunk] - widthR/2
    }
    theseRs = rsIntegrationPoints
    thesePointsRadial = make.surface.grid(list(rs=theseRs, as=theseAs))
    
    # convert to Euclidean coordinates
    thesePointsEuclidean = cbind(thesePointsRadial[,1]*cos(thesePointsRadial[,2]), 
                                 thesePointsRadial[,1]*sin(thesePointsRadial[,2]))
    
    # translate coordinates based on the jittered observation coordinates
    return(thesePointsEuclidean)
  }
  
  # for every ring:
  #   for every point:
  #     get sub-integration points
  subWs = list()
  subPts = list()
  for(i in 1:length(rs)) {
    thisIntW = ws[i]
    theseSubWs = rep(thisIntW/(nSubRPerPoint*nSubAPerPoint), 
                     each=nSubRPerPoint*nSubAPerPoint*nrow(pts[[i]]))
    
    theseas = as[[i]]
    theseMinR = ifelse(i==1, R1, rs[i-1])
    theseMaxR = rs[i]
    if(length(theseas) != 1) {
      aWidth = theseas[2] - theseas[1]
    } else {
      aWidth = 2*pi
    }
    theseSubPts = c()
    
    #browser()
    for(j in 1:length(theseas)) {
      # determine boundaries of this integration area
      thisMinA = theseas[j] - aWidth/2
      thisMaxA = theseas[j] + aWidth/2
      thisMinR = theseMinR
      thisMaxR = theseMaxR
      
      # obtain sub-integration points
      thisSubPts = getSubIntegrationPointsForOnePoint(minR=thisMinR, maxR=thisMaxR, 
                                                      minA=thisMinA, maxA=thisMaxA)
      theseSubPts = rbind(theseSubPts, thisSubPts)
    }
    
    subPts = c(subPts, list(theseSubPts))
    subWs = c(subWs, list(theseSubWs))
  }
  
  return(list(subPts=subPts, subWs=subWs))
}
  
# Testing ----
# res <- make.integration.donut(coords = matrix(c(-3328, 994), nrow = 4, ncol = 2, byrow = TRUE),
#                                1, 2, 5, 10, nigeria.adm2) 
# 
# res$w
# rbPal <- colorRampPalette(c('red','blue'))
# pal <-  rbPal(10)[as.numeric(cut(res$w[1, ],breaks = 10))]
# df <- data.frame(x = res$x[1,], y = res$y[1, ], col = res$w[1, ])
# ggplot(df) + geom_point(aes(x = x, y = y, color = col))
# 
# plot(res$pts$pts5)
# points(res$pts$pts4)
# points(res$pts$pts3)
# points(res$pts$pts2)
# points(res$pts$pts1)
# 
# set.seed(2)
# adm2 <- readOGR("~/Github/master/gadm41_NGA_shp/gadm41_NGA_2.shp")
# adm2.sf <- st_read("~/Github/master/gadm41_NGA_shp/gadm41_NGA_2.shp")
# 
# locs <- (simulate.locations(100, p.urban = 0.2, uniform = TRUE))$obs.locs
# colnames(locs)[1:2] <- c("east", "north") 
# true.lat.lon <- convert.km.to.deg(locs[,c(1, 2)])
# true.lat.lon <- SpatialPoints(true.lat.lon, proj4string = adm2@proj4string, bbox = NULL)
# check1 <- over(true.lat.lon, adm2, returnList = FALSE)
# obs.locs <- jitter.locs(scale = c(1),
#                         locKM = locs[, c(1, 2)],
#                         urbanRural = locs$urban,
#                         shapeFile = adm2,
#                         check1 = check1,
#                         boundary = TRUE,
#                         adm.level = 2)
# obs.locs <- obs.locs[[1]][, c(1, 2)]
# 
# res2 <- make.intergration.donut(coords = obs.locs,
#                                 1, 2, 5, 10, nigeria.adm2) 
# 
# 
# adminMapPolygons = as.SpatialPolygons.PolygonsList(nigeria.adm2@polygons, nigeria.adm2@proj4string)
# adminMapPolygonsEN = spTransform(adminMapPolygons, CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"))
# adm2.fort <- fortify(adminMapPolygonsEN, region = "ID")
# id = 8 # 8 good
# df2 <- data.frame(x = res2$x[id,], y = res2$y[id, ], col = res2$w[id, ])
# b = 0 # buffer
# ggplot(df2) + geom_point(aes(x = x, y = y, color = col)) + 
#   geom_polygon(data = adm2.fort, aes(x = long, y = lat, group = group), colour='black', fill=NA, size = 0.1) +
#   coord_fixed(ratio = 1, xlim = c(min(res2$x[id,]) - b, max(res2$x[id, ]) + b), ylim = c(min(res2$y[id, ]) - b, max(res2$y[id, ]) + b))
# 
# pts.ll <- convert.km.to.deg(obs.locs)
# sp.pt <- SpatialPoints(pts.ll, nigeria.adm2@proj4string)
# tmp <- over(sp.pt, nigeria.adm2, returnList = FALSE)
# admin.id <- tmp$GID_2
# 
# counties <- nigeria.adm2$GID_2
# dict <- c()
# for(i in 1:length(counties)){
#   dict[counties[i]] = i
# }
# dict[admin.id[2]]
# sum(res2$w[id, ])
# 
# adminMapPolygons = as.SpatialPolygons.PolygonsList(nigeria.adm2@polygons, nigeria.adm2@proj4string)
# plot(adminMapPolygons[495])
# 
# # test if old function adjusts 2. coordinate
# test2 <- makeAllIntegrationPoints(obs.locs, c(rep(TRUE, 50), rep(FALSE, 50)), adminMap = nigeria.adm2)

# Update weights for points that have are outside administrative area

  

