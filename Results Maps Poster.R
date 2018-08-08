# Figure 1
#4 panels - apparent prevalence, median prevalence, max prevalence, median variance
#need to wait for winning model to finalize
#currently for CAR model with a queens W matrix

library(raster)
#library(tmap)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(rgdal)
library(rgeos)
library(gridExtra)

##########################################################################################
#prepare spatial data
#########################################################################################
ai <- readRDS("~/Github/location_samplingevent_n_y_speciesgroup.rds")
ai <- ai[!(ai$collection.year %in% c("2011", "2014", "2016", "2017")), ]

#add coluns to tell which watershed
setwd("~/HP/hydrologic_units")
huc4 <- shapefile("huc4.shp")
projection(huc4) <- CRS("+proj=longlat +ellps=WGS84")
pt <- ai[, c("lat", "long")]
coordinates(pt) <- ~ long + lat
proj4string(pt) <- CRS("+proj=longlat +ellps=WGS84")

#select only wanted hucs
sel <- !(huc4$STATES %in% c("AK", "AS", "AK,CN", "HI", "PR", "GU", "MP", "VI"))
huc4 <- huc4[sel, ]

#fill columns
ai$huc4 <- extract(huc4, pt)$HUC4
ai$huc4name <- extract(huc4, pt)$NAME

#number of sites
sites <- unique(ai$huc4)
sites <- sites[-17]

######################################################################################
#read in model output (CAR/queens)
#will need to change to reflect winning model
######################################################################################

#from workflow.R
nmonths <- 12
#set nyears = 5 if not using 2016 data
nyears <- 5
nsites <- length(sites) #should be base=136, spatial=202
#set niter and thin the match data set
niter <- 10000
thin <- 10
niter <- niter/thin
niter <- niter/2
ncol <- 2+(nmonths*nyears*nsites)

#read in results to be used
setwd("~/Honors Thesis/Thesis Work/Model Selection")
data <- readRDS("icar_wtnetwork_coda_10000iterations.rds")
data.1 <- as.matrix(data[[1]][1:niter, ])
data.2 <- as.matrix(data[[2]][1:niter, ])
data.3 <- as.matrix(data[[3]][1:niter, ])
test <- rbind(data.1, data.2, data.3)
#will need to change dimensions for models with 202 hucs
test <- test[,2:ncol]
# rm(data, data.1, data.2, data.3)

######################################################################################################
#quantile calculations
###############################################################################################
#create an array that holds all interations for each month/year combination, across all hucs
c <- array(0, dim=c((3*niter), nyears*nmonths, nsites))
for (i in 1:nsites) {
  for (j in 1:(nyears*nmonths)){
    k <- ((i-1)*60)+j #should be 2+((i-1)*60)+j if mapping without 2016 data
    hold <- test[,k]
    c[ ,j,i] <- hold
    #rm(hold)
  }
}

#produce an array to hold the median/95th/maximum for each iteration
#commented out maximum
med <- array(0, dim=c(1, (3*niter), nsites))
ninetyfive <- array(0, dim=c(1, (3*niter), nsites))
# max <- array(0, dim=c(1, (3*(niter/thin)), nsites))
for (i in 1:nsites) {
  for (j in 1:(3*niter)) {
    hold <- quantile(c[j, ,i], probs=0.5)
    hold.n <- quantile(c[j, ,i], probs=0.95)
    # hold.max <- quantile(c[j, ,i], probs=1)
    med[ ,j,i] <- hold
    ninetyfive[ ,j,i] <- hold.n
    # max[ ,j,i] <- hold.max
    #rm(hold, hold.n)
  }
}

#find the median value of across all interations for each level of measurement
#store the values in a column vector
pi.med <- matrix(0, nsites, 1)
pi.ninetyfive <- matrix(0, nsites, 1)
# pi.max <- matrix(0, nsites, 1)
var.med <- matrix(0, nsites, 1)
for (i in 1:nsites) {
  pi.med[i, ] <- median(med[1, ,i])
  var.med[i, ] <- var(med[1, ,i])
  pi.ninetyfive[i, ] <- median(ninetyfive[1, ,i])
  # pi.max[i, ] <- median(max[1, ,i])
}

#rm(med, ninetyfive, )
#create data frame that holds location information and true prevalence data
prevalence <- data.frame(hucname = sites, hucnumber = seq(1, 181),
                   median = pi.med, ninetyfive = pi.ninetyfive,
                   variance = var.med)   # max = pi.max,)

##########################################################################################
#matching prevalence to huc4
##########################################################################################

#select years in of ai data used in fitting the model
ai <- readRDS("~/Honors Thesis/Thesis Work/HP_dataset.rds")
aired <- ai[!(ai$collection.year %in% c("2011", "2014", "2017")), ]
ai <- SpatialPointsDataFrame(coords = cbind(aired$long, aired$lat),
                             proj4string = huc4@proj4string, data=aired)

#read in watersheds as an OGR
setwd("~/Honors Thesis/Project/hydrologic_units")
huc4<-readOGR(dsn=".", layer="huc4")
projection(huc4) <- CRS("+proj=longlat +ellps=WGS84")

#fill columns
huc4$median <- prevalence$median[match(huc4$HUC4, prevalence$hucname)]
huc4$var.med <- prevalence$variance[match(huc4$HUC4, prevalence$hucname)]
huc4$ninetyfive <- prevalence$ninetyfive[match(huc4$HUC4, prevalence$hucname)]
# huc4$maximum <- prevalence$max[match(huc4$HUC4, prevalence$hucname)]

#calculate apparent prevalence and fill columns
test <- over(ai, huc4)
allsamples <- as.data.frame(table(test$HUC4))
allsamples$Var1 <- as.character(allsamples$Var1)
huc4$num.samples <- allsamples$Freq[match(huc4$HUC4, allsamples$Var1)]
neg <- ai$AIpcr_susneg == "negative"
test <- over(ai[!neg, ], huc4)
possamples <- as.data.frame(table(test$HUC4))
possamples$Var1 <- as.character(possamples$Var1)
huc4$num.pos.samples <- possamples$Freq[match(huc4$HUC4, possamples$Var1)]
huc4$num.pos.samples[is.na(huc4$num.pos.samples)] <- 0
huc4$apparent.prev <- huc4$num.pos.samples/huc4$num.samples

##########################################################################################
#make some maps
##########################################################################################

#simplify data for faster plotting
temp <- gSimplify(huc4, tol=0.01, topologyPreserve = TRUE)
plot.huc <- SpatialPolygonsDataFrame(temp, data=huc4@data)
plot.huc <- plot.huc[sel, ]
plot <- fortify(plot.huc)
plot.huc$id <- row.names(plot.huc)
plot.huc <- left_join(plot, plot.huc@data)

pal.2 <- c(brewer.pal(3, "Blues"), "white")
pal <- brewer.pal(3, "Blues")

#apparent prev
plot.huc$apparent.prev[is.na(plot.huc$apparent.prev)] <- 1
qa <- quantile(plot.huc$apparent.prev[!(plot.huc$apparent.prev == "1")], c(0, 0.33, 0.66, 1))
qa.a <- quantile(plot.huc$apparent.prev, c(0, 0.33, 0.66, 1))
qa
plot.huc$quantile1 <- cut(plot.huc$apparent.prev, 
                           c(qa[1], qa[2], qa[3], qa[4], qa.a[4]),
                           labels = c("0-0.077", "0.0077-0.223", "0.223-0.857", "No Samples"),
                           include.lowest = TRUE)

#Median
qa <- quantile(plot.huc$median, c(0, 0.33, 0.66, 1), na.rm=TRUE)
qa
plot.huc$quantile2 <- cut(plot.huc$median, 
                           c(qa[1], qa[2], qa[3], qa[4]),
                           labels = c("0.009-0.060", "0.060-0.119", "0.119-0.510"),
                           include.lowest = TRUE)

#95%
qa <- quantile(plot.huc$ninetyfive, c(0, 0.33, 0.66, 1), na.rm=TRUE)
qa
plot.huc$quantile3 <- cut(plot.huc$ninetyfive, c(qa[1], qa[2], qa[3], qa[4]),
                          labels = c("0.114-0.421", "0.421-0.612", "0.612-0.922"), 
                          include.lowest = TRUE)
#variance
qa <- quantile(plot.huc$var.med, c(0, 0.33, 0.66, 1), na.rm=TRUE)
qa
plot.huc$quantile4 <- cut(plot.huc$var.med, c(qa[1], qa[2], qa[3], qa[4]), 
                          labels = c("0.00002-0.0003", "0.0003-0.0007", "0.0007-0.01"),
                          include.lowest = TRUE)
#apparent prevalence map
apparent.plot <- ggplot() +
  geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, 
                                  fill=quantile1),
               colour = "black") +
  scale_fill_manual(values = pal.2) +
  coord_fixed(1.3) +
  theme(panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), 
        legend.title = element_blank()) +
  ggtitle("Apparent Prevalence")

jpeg("~/Honors Thesis/Thesis Work/Paper Figures/Figure 1a.jpeg")
apparent.plot
dev.off()

#median map
median.plot <- ggplot() +
  geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, fill = quantile2), 
               colour = "black") +
  scale_fill_manual(values = pal) +
  coord_fixed(1.3) +
  theme(panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), 
        legend.title = element_blank()) +
  ggtitle("Median True Prevalence")

jpeg("~/Honors Thesis/Thesis Work/Paper Figures/Figure 1b.jpeg")
median.plot
dev.off()

#95% map
ninetyfive.plot <- ggplot() +
  geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, fill = quantile3), 
               colour = "black") +
  scale_fill_manual(values = pal) +
  coord_fixed(1.3) +
  theme(panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size=10)) +
  ggtitle("95th-Percentile True Prevalence")

jpeg("~/Honors Thesis/Thesis Work/Paper Figures/Figure 1c.jpeg")
ninetyfive.plot
dev.off()

#variance map
variance.plot <- ggplot() +
  geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, 
                                  fill = quantile4), colour = "black") +
  scale_fill_manual(values = pal) +
  coord_fixed(1.3) +
  theme(panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), 
        legend.title = element_blank()) +
  ggtitle("Variance")

jpeg("~/Honors Thesis/Thesis Work/Paper Figures/Figure 1d.jpeg")
variance.plot
dev.off()

#composite
jpeg("~/Honors Thesis/Thesis Work/Paper Figures/Figure 1.jpeg")
grid.arrange(apparent.plot, median.plot, ninetyfive.plot, variance.plot, ncol=2)
dev.off()
