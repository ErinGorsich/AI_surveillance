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
ai <- readRDS("~/Github/AVHS_samplingevent_locationname.rds")
ai <- ai[!(ai$collection.year %in% c("2011", "2014", "2016", "2017")), ]
ai <- ai[ai$species.group %in% c("1", "2", "4"), ]
ai$species.group[ai$species.group == "4"] <- 3

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
# sites <- sites[-17]

######################################################################################
#read in model output (CAR/queens)
#will need to change to reflect winning model
######################################################################################

#from workflow.R
nmonths <- 12
#set nyears = 5 if not using 2016 data
nyears <- 5
nsites <- length(sites)-1 #194
#set niter and thin the match data set
niter <- 80000
thin <- 10
niter <- niter/thin
ncol <- 6+(nmonths*nyears*nsites)

#read in results to be used
setwd("~/Github")
data <- readRDS("base_fit.rds")
data.11 <- as.matrix(data[[1]][1:niter, 7:11646])
data.21 <- as.matrix(data[[2]][1:niter, 7:11646])
data.31 <- as.matrix(data[[3]][1:niter, 7:11646])
data.12 <- as.matrix(data[[1]][1:niter, 11647:23286])
data.22 <- as.matrix(data[[2]][1:niter, 11647:23286])
data.32 <- as.matrix(data[[3]][1:niter, 11647:23286])
data.13 <- as.matrix(data[[1]][1:niter, 23287:34926])
data.23 <- as.matrix(data[[2]][1:niter, 23287:34926])
data.33 <- as.matrix(data[[3]][1:niter, 23287:34926])

test.1 <- rbind(data.11, data.21, data.31)
test.2 <- rbind(data.12, data.22, data.32)
test.3 <- rbind(data.13, data.23, data.33)

rm(data, data.11, data.12, data.13, data.21, data.22, data.23, data.31, data.32, data.33)

######################################################################################################
#quantile calculations - dabbling
###############################################################################################
#create an array that holds all interations for each month/year combination, across all hucs
c <- array(0, dim=c((3*niter), nyears*nmonths, nsites))
for (i in 1:nsites) {
  for (j in 1:(nyears*nmonths)){
    k <- ((i-1)*60)+j #should be 2+((i-1)*60)+j if mapping without 2016 data
    hold <- test.1[,k]
    c[ ,j,i] <- hold
    #rm(hold)
  }
}

#produce an array to hold the median/95th/maximum for each iteration
#commented out maximum
med <- array(0, dim=c(1, (3*niter), nsites))
for (i in 1:nsites) {
  for (j in 1:(3*niter)) {
    hold <- quantile(c[j, ,i], probs=0.5)
    med[ ,j,i] <- hold
    #rm(hold)
  }
}

#find the median value of across all interations for each level of measurement
#store the values in a column vector
pi.med <- matrix(0, nsites, 1)
for (i in 1:nsites) {
  pi.med[i, ] <- median(med[1, ,i])
}

#rm(med)
#create data frame that holds location information and true prevalence data
prevalence.dabbling <- data.frame(hucname = sites[1:194], hucnumber = seq(1, 194),
                                  median = pi.med)

######################################################################################################
#quantile calculations - diving
###############################################################################################
#create an array that holds all interations for each month/year combination, across all hucs
c <- array(0, dim=c((3*niter), nyears*nmonths, nsites))
for (i in 1:nsites) {
  for (j in 1:(nyears*nmonths)){
    k <- ((i-1)*60)+j #should be 2+((i-1)*60)+j if mapping without 2016 data
    hold <- test.2[,k]
    c[ ,j,i] <- hold
    #rm(hold)
  }
}

#produce an array to hold the median/95th/maximum for each iteration
#commented out maximum
med <- array(0, dim=c(1, (3*niter), nsites))
for (i in 1:nsites) {
  for (j in 1:(3*niter)) {
    hold <- quantile(c[j, ,i], probs=0.5)
    med[ ,j,i] <- hold
    #rm(hold)
  }
}

#find the median value of across all interations for each level of measurement
#store the values in a column vector
pi.med <- matrix(0, nsites, 1)
for (i in 1:nsites) {
  pi.med[i, ] <- median(med[1, ,i])
}

#rm(med)
#create data frame that holds location information and true prevalence data
prevalence.diving <- data.frame(hucname = sites[1:194], hucnumber = seq(1, 194),
                                  median = pi.med)

######################################################################################################
#quantile calculations - geese
###############################################################################################
#create an array that holds all interations for each month/year combination, across all hucs
c <- array(0, dim=c((3*niter), nyears*nmonths, nsites))
for (i in 1:nsites) {
  for (j in 1:(nyears*nmonths)){
    k <- ((i-1)*60)+j #should be 2+((i-1)*60)+j if mapping without 2016 data
    hold <- test.3[,k]
    c[ ,j,i] <- hold
    #rm(hold)
  }
}

#produce an array to hold the median/95th/maximum for each iteration
#commented out maximum
med <- array(0, dim=c(1, (3*niter), nsites))
for (i in 1:nsites) {
  for (j in 1:(3*niter)) {
    hold <- quantile(c[j, ,i], probs=0.5)
    med[ ,j,i] <- hold
    #rm(hold)
  }
}

#find the median value of across all interations for each level of measurement
#store the values in a column vector
pi.med <- matrix(0, nsites, 1)
for (i in 1:nsites) {
  pi.med[i, ] <- median(med[1, ,i])
}

#rm(med)
#create data frame that holds location information and true prevalence data
prevalence.geese <- data.frame(hucname = sites[1:194], hucnumber = seq(1, 194),
                                  median = pi.med)

##########################################################################################
#matching prevalence to huc4
##########################################################################################

#select years in of ai data used in fitting the model
ai <- readRDS("~/Github/AVHS_samplingevent_locationname.rds")
aired <- ai[!(ai$collection.year %in% c("2011", "2014", "2016", "2017")), ]
aired <- aired[aired$species.group %in% c("1", "2", "4"), ]
aired$species.group[aired$species.group == "4"] <- 3
ai <- SpatialPointsDataFrame(coords = cbind(aired$long, aired$lat),
                             proj4string = huc4@proj4string, data=aired)

#read in watersheds as an OGR
setwd("~/HP/hydrologic_units")
huc4<-readOGR(dsn=".", layer="huc4")
projection(huc4) <- CRS("+proj=longlat +ellps=WGS84")

#fill columns
huc4$median.dabbling <- prevalence.dabbling$median[match(huc4$HUC4, prevalence.dabbling$hucname)]
huc4$median.diving <- prevalence.diving$median[match(huc4$HUC4, prevalence.diving$hucname)]
huc4$median.geese <- prevalence.geese$median[match(huc4$HUC4, prevalence.geese$hucname)]

#dabbling apparent prevalence
ai.dabbling <- ai[ai$species.group == "1", ]
test <- over(ai.dabbling, huc4)
allsamples <- as.data.frame(table(test$HUC4))
allsamples$Var1 <- as.character(allsamples$Var1)
huc4$num.samples.dabbling <- allsamples$Freq[match(huc4$HUC4, allsamples$Var1)]
neg <- ai.dabbling$AIpcr_susneg == "negative"
test <- over(ai.dabbling[!neg, ], huc4)
possamples <- as.data.frame(table(test$HUC4))
possamples$Var1 <- as.character(possamples$Var1)
huc4$num.pos.samples.dabbling <- possamples$Freq[match(huc4$HUC4, possamples$Var1)]
huc4$num.pos.samplesdabbling[is.na(huc4$num.pos.samples.dabbling)] <- 0
huc4$apparent.prev.dabbling <- huc4$num.pos.samples.dabbling/huc4$num.samples.dabbling

#diving apparent prevalence
ai.diving <- ai[ai$species.group == "2", ]
test <- over(ai.diving, huc4)
allsamples <- as.data.frame(table(test$HUC4))
allsamples$Var1 <- as.character(allsamples$Var1)
huc4$num.samples.diving <- allsamples$Freq[match(huc4$HUC4, allsamples$Var1)]
neg <- ai.diving$AIpcr_susneg == "negative"
test <- over(ai.diving[!neg, ], huc4)
possamples <- as.data.frame(table(test$HUC4))
possamples$Var1 <- as.character(possamples$Var1)
huc4$num.pos.samples.diving <- possamples$Freq[match(huc4$HUC4, possamples$Var1)]
huc4$num.pos.samples.diving[is.na(huc4$num.pos.samples.diving)] <- 0
huc4$apparent.prev.diving <- huc4$num.pos.samples.diving/huc4$num.samples.diving

#geese apparent prevalence
ai.geese <- ai[ai$species.group == "3", ]
test <- over(ai.geese, huc4)
allsamples <- as.data.frame(table(test$HUC4))
allsamples$Var1 <- as.character(allsamples$Var1)
huc4$num.samples.geese <- allsamples$Freq[match(huc4$HUC4, allsamples$Var1)]
neg <- ai.geese$AIpcr_susneg == "negative"
test <- over(ai.geese[!neg, ], huc4)
possamples <- as.data.frame(table(test$HUC4))
possamples$Var1 <- as.character(possamples$Var1)
huc4$num.pos.samples.geese <- possamples$Freq[match(huc4$HUC4, possamples$Var1)]
huc4$num.pos.samples.geese[is.na(huc4$num.pos.samples.geese)] <- 0
huc4$apparent.prev.geese <- huc4$num.pos.samples.geese/huc4$num.samples.geese

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

#color pallete
pal.2 <- c(brewer.pal(6, "Blues")[2:6], "white")
pal <- brewer.pal(3, "Blues")

#set a common scale for easy comparison
plot.huc$scale.apdabbling <- cut(plot.huc$apparent.prev.dabbling, c(0, 0.1, 0.2, 0.35, 0.5, 1),
                      labels = c("0-0.1", "0.1-0.2", "0.2-0.35", "0.35-0.5", "0.5-0.1"), 
                      include.lowest = TRUE)
plot.huc$scale.apdiving <- cut(plot.huc$apparent.prev.diving, c(0, 0.1, 0.2, 0.35, 0.5, 1),
                               labels = c("0-0.1", "0.1-0.2", "0.2-0.35", "0.35-0.5", "0.5-0.1"), 
                               include.lowest = TRUE)
plot.huc$scale.apgeese <- cut(plot.huc$apparent.prev.geese, c(0, 0.1, 0.2, 0.35, 0.5, 1),
                               labels = c("0-0.1", "0.1-0.2", "0.2-0.35", "0.35-0.5", "0.5-0.1"), 
                               include.lowest = TRUE)
plot.huc$scale.meddabbling <- cut(plot.huc$median.dabbling, c(0, 0.1, 0.2, 0.35, 0.5, 1),
                               labels = c("0-0.1", "0.1-0.2", "0.2-0.35", "0.35-0.5", "0.5-0.1"), 
                               include.lowest = TRUE)
plot.huc$scale.meddiving <- cut(plot.huc$median.diving, c(0, 0.1, 0.2, 0.35, 0.5, 1),
                               labels = c("0-0.1", "0.1-0.2", "0.2-0.35", "0.35-0.5", "0.5-0.1"), 
                               include.lowest = TRUE)
plot.huc$scale.medgeese <- cut(plot.huc$median.geese, c(0, 0.1, 0.2, 0.35, 0.5, 1),
                              labels = c("0-0.1", "0.1-0.2", "0.2-0.35", "0.35-0.5", "0.5-0.1"), 
                              include.lowest = TRUE)
# #apparent prev
# plot.huc$apparent.prev[is.na(plot.huc$apparent.prev)] <- 1
# qa <- quantile(plot.huc$apparent.prev[!(plot.huc$apparent.prev == "1")], c(0, 0.33, 0.66, 1))
# qa.a <- quantile(plot.huc$apparent.prev, c(0, 0.33, 0.66, 1))
# qa
# plot.huc$quantile1 <- cut(plot.huc$apparent.prev, 
#                            c(qa[1], qa[2], qa[3], qa[4], qa.a[4]),
#                            labels = c("0-0.077", "0.0077-0.223", "0.223-0.857", "No Samples"),
#                            include.lowest = TRUE)
# 
# #Median
# qa <- quantile(plot.huc$median, c(0, 0.33, 0.66, 1), na.rm=TRUE)
# qa
# plot.huc$quantile2 <- cut(plot.huc$median, 
#                            c(qa[1], qa[2], qa[3], qa[4]),
#                            labels = c("0.009-0.060", "0.060-0.119", "0.119-0.510"),
#                            include.lowest = TRUE)
# 
# #95%
# qa <- quantile(plot.huc$ninetyfive, c(0, 0.33, 0.66, 1), na.rm=TRUE)
# qa
# plot.huc$quantile3 <- cut(plot.huc$ninetyfive, c(qa[1], qa[2], qa[3], qa[4]),
#                           labels = c("0.114-0.421", "0.421-0.612", "0.612-0.922"), 
#                           include.lowest = TRUE)
# #variance
# qa <- quantile(plot.huc$var.med, c(0, 0.33, 0.66, 1), na.rm=TRUE)
# qa
# plot.huc$quantile4 <- cut(plot.huc$var.med, c(qa[1], qa[2], qa[3], qa[4]), 
#                           labels = c("0.00002-0.0003", "0.0003-0.0007", "0.0007-0.01"),
#                           include.lowest = TRUE)

#dabbling apparent prevalence map
dabbling.apparent.plot <- ggplot() +
  geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, 
                                  fill=scale.apdabbling),
               colour = "black") +
  scale_fill_manual(values = pal.2) +
  coord_fixed(1.3) +
  theme(panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), 
        legend.title = element_blank()) +
  ggtitle("Apparent Prevalence - Species Group 1")

#diving apparent prevalence map
diving.apparent.plot <- ggplot() +
  geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, 
                                  fill=scale.apdiving),
               colour = "black") +
  scale_fill_manual(values = pal.2) +
  coord_fixed(1.3) +
  theme(panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), 
        legend.title = element_blank()) +
  ggtitle("Apparent Prevalence - Species Group 2")

#geese apparent prevalence map
geese.apparent.plot <- ggplot() +
  geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, 
                                  fill=scale.apgeese),
               colour = "black") +
  scale_fill_manual(values = pal.2) +
  coord_fixed(1.3) +
  theme(panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), 
        legend.title = element_blank()) +
  ggtitle("Apparent Prevalence - Species Group 3")

#dabbling median map
dabbling.median.plot <- ggplot() +
  geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, 
                                  fill=scale.meddabbling),
               colour = "black") +
  scale_fill_manual(values = pal.2) +
  coord_fixed(1.3) +
  theme(panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), 
        legend.title = element_blank()) +
  ggtitle("Median True Prevalence - Species Group 1")

#diving median map
diving.median.plot <- ggplot() +
  geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, 
                                  fill=scale.meddiving),
               colour = "black") +
  scale_fill_manual(values = pal.2) +
  coord_fixed(1.3) +
  theme(panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), 
        legend.title = element_blank()) +
  ggtitle("Median True Prevalence - Species Group 2")

#geese median map
geese.median.plot <- ggplot() +
  geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, 
                                  fill=scale.medgeese),
               colour = "black") +
  scale_fill_manual(values = pal.2) +
  coord_fixed(1.3) +
  theme(panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), 
        legend.title = element_blank()) +
  ggtitle("Median True Prevalence - Species Group 3")

setwd("~/HP/Paper Figures")
jpeg("dabbling_apparent_map_base.jpeg")
dabbling.apparent.plot
dev.off()

jpeg("diving_apparent_map_base.jpeg")
diving.apparent.plot
dev.off()

jpeg("geese_apparent_map_base.jpeg")
geese.apparent.plot
dev.off()

jpeg("dabbling_median_map_base.jpeg")
dabbling.median.plot
dev.off()

jpeg("diving_median_map_base.jpeg")
diving.median.plot
dev.off()

jpeg("geese_median_map_base.jpeg")
geese.median.plot
dev.off()

# jpeg("~/Honors Thesis/Thesis Work/Paper Figures/Figure 1a.jpeg")
# apparent.plot
# dev.off()

# #median map
# median.plot <- ggplot() +
#   geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, fill = quantile2), 
#                colour = "black") +
#   scale_fill_manual(values = pal) +
#   coord_fixed(1.3) +
#   theme(panel.background=element_blank(), panel.grid.major = element_blank(),
#         panel.border = element_blank(), axis.title = element_blank(),
#         axis.ticks = element_blank(), axis.text = element_blank(), 
#         legend.title = element_blank()) +
#   ggtitle("Median True Prevalence")
# 
# jpeg("~/Honors Thesis/Thesis Work/Paper Figures/Figure 1b.jpeg")
# median.plot
# dev.off()
# 
# #95% map
# ninetyfive.plot <- ggplot() +
#   geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, fill = quantile3), 
#                colour = "black") +
#   scale_fill_manual(values = pal) +
#   coord_fixed(1.3) +
#   theme(panel.background=element_blank(), panel.grid.major = element_blank(),
#         panel.border = element_blank(), axis.title = element_blank(),
#         axis.ticks = element_blank(), axis.text = element_blank(), 
#         legend.title = element_blank(),
#         legend.text = element_text(size=10)) +
#   ggtitle("95th-Percentile True Prevalence")
# 
# jpeg("~/Honors Thesis/Thesis Work/Paper Figures/Figure 1c.jpeg")
# ninetyfive.plot
# dev.off()
# 
# #variance map
# variance.plot <- ggplot() +
#   geom_polygon(data=plot.huc, aes(x=long, y=lat, group=group, 
#                                   fill = quantile4), colour = "black") +
#   scale_fill_manual(values = pal) +
#   coord_fixed(1.3) +
#   theme(panel.background=element_blank(), panel.grid.major = element_blank(),
#         panel.border = element_blank(), axis.title = element_blank(),
#         axis.ticks = element_blank(), axis.text = element_blank(), 
#         legend.title = element_blank()) +
#   ggtitle("Variance")
# 
# jpeg("~/Honors Thesis/Thesis Work/Paper Figures/Figure 1d.jpeg")
# variance.plot
# dev.off()
# 
# #composite
# jpeg("~/Honors Thesis/Thesis Work/Paper Figures/Figure 1.jpeg")
# grid.arrange(apparent.plot, median.plot, ninetyfive.plot, variance.plot, ncol=2)
# dev.off()
