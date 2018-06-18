library(rgdal)
library(sp)
library(maptools)
library(rgeos)
library(maps)
library(raster)
library(geosphere)

# need hucs in global environment; needs network in global environment
if (Erin){
    hucdir <- "~/Documents/Avian_Influenza/Data/hydrologic_units/"
} else {
    hucdir <- "~/Honors Thesis/Project/hydrologic_units/"
}
huc4 <- shapefile(paste(hucdir, "huc4.shp", sep = ""))

# HP - delete me - if read in works ok                  
# setwd("~/Honors Thesis/Project/hydrologic_units")
# huc4 <- shapefile("huc4.shp")
projection(huc4) <- CRS("+proj=longlat +ellps=WGS84")
sel <- !(huc4$STATES %in% c("AK", "AS", "AK,CN", "HI", "PR", "GU", "MP", "VI"))
huc4 <- huc4[sel, ]
nsites <- length(huc4@polygons)  # n = 202
names <- huc4@data$HUC4


#needs ai data with origin and destination watersheds as columns
#t <- readRDS("AVHS_Band_v4.rds") # consider if needs made from just band data
#t <- t[!is.na(t$o.huc4), ]
#t <- t[!is.na(t$d.huc4), ]

# update based on bird movement data
# setwd("~/HP/Data")
if (!Erin){
setwd("~/Github")
}
t <- readRDS("birdmovementnets.rds")

# function to define shape of autocorrelation in space
############################################
############################################
# Functions to think about autocorreltion in space
############################################
############################################
scaleW <- function(W){
    B <- matrix(0, nrow = dim(W)[1], ncol = dim(W)[1],
        dimnames = list(rownames(W), colnames(W)))
    for(i in 1:dim(W)[1]){
        if (sum(W[i, ]) > 0){
            B[i, ] <- W[i, ]/ sum(W[i, ])
        }
    }
    return(B)
}

define.neighborhood <- function(method){
    # function returns matrix defining neighborhood
    # row = target watershed: column = is this a neighbor of target watershed
    # requires huc4 defined globally; and t
    nsites <- length(huc4@polygons)
    names <- huc4@data$HUC4
    W <- matrix(0, nrow = nsites, ncol = nsites, 
                dimnames = list(names, names))
    if (method == "queens") {
        # defines neighbors as wheter or not the polygons are touching 
        for (i in 1:nsites){
            W[i, ] <-  gIntersects(spgeom1 = huc4, spgeom2 = huc4[i,], byid = T)
            W[i, i] <- 0
        }
    }
    
    if (method == "weightedqueens"){
        # weight = length of shared boundaries between the watersheds
        touches <- matrix(0, nrow = nsites, ncol = nsites, 
                          dimnames = list(names, names))
        for (i in 1:nsites){
            touches[i, ] <- gIntersects(spgeom1 = huc4, spgeom2 = huc4[i,], byid = T)
            touches[i, i] <- 0
        }
        for (i in 1:nsites){
            for (j in 1:nsites){
                if (touches[i, j] == 0){
                    W[i, j] <- 0
                }
                if (touches[i, j] == 1){
                    line <- gIntersection(huc4[i,], huc4[j,])
                    ifelse(class(line) == "SpatialLines", 
                           W[i, j] <- SpatialLinesLengths(line, longlat = FALSE),
                           W[i, j] <- sum(SpatialLinesLengths(line@lineobj, 
                                longlat = FALSE)) )
                }
            }
        }
    }
    
    if (method == "knn3"){
        # calculate matrix as 3 nearest neighbors
        dist <- matrix(0, nrow = nsites, ncol = nsites, 
                       dimnames = list(names, names))
        centroids <- gCentroid(huc4, byid = TRUE)
        for (i in 1:nsites){
            for (j in 1:nsites){
                dist[i,j] <- distVincentyEllipsoid(p1 = centroids@coords[i, 1:2], 
                                                   p2 = centroids@coords[j, 1:2])
            }
            minvals <- sort(dist[i, ])[2:4]
            W[i, which(dist[i, ] %in% minvals)] <- 1
        }
    }
    if (method == "network") {
        for (i in 1:length(colnames(W))) {
            for (j in 1:length(rownames(W))) {
                temp <- t[t$o.huc4 == colnames(W)[i] & 
                    t$d.huc4 == rownames(W)[j], ]
                temp2 <- t[t$d.huc4 == colnames(W)[i] & 
                    t$o.huc4 == rownames(W)[j], ]
                temp <- length(temp[,1]) + length(temp2[,1])
                ifelse(temp < 1, W[i, j] <- 0, W[i, j] <- 1)
               # ifelse(length(t$sampleid[
               #     t$o.huc4 %in% c(rownames(W)[i], rownames(W)[j]) &
               #     t$d.huc4 %in% c(rownames(W)[i], rownames(W)[j]) ]) < 1,
               #     W[i, j] <- 0,
               #     W[i, j] <- 1)
                W[i, i] <- 0
            }
        }
    }
    if (method == "weightednetwork") {
        for (i in 1:length(colnames(W))) {
            for (j in 1:length(rownames(W))) {
                temp <- t[t$o.huc4 == colnames(W)[i] & 
                    t$d.huc4 == rownames(W)[j], ]
                temp2 <- t[t$d.huc4 == colnames(W)[i] & 
                    t$o.huc4 == rownames(W)[j], ]
                W[i, j] <- length(temp[,1]) + length(temp2[,1])
                W[i, i] <- 0
            }
        }
    }
    return(W)
}



