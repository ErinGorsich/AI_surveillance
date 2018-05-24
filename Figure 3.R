#figure 3
#four paneled figure to map the different ways of making the W matrix

library(maps)
library(maptools)
library(rgdal)
library(ggplot2)
library(rgeos)
library(RColorBrewer)
library(dplyr)

setwd("~/GitHub/HP_AI_Surveillance")

source("define.neighborhood.R")  # note was called define_neighborhood?
source("multiplot.R")

#prepare huc spatial data
huc <- readOGR(dsn = "hydrologic_units", layer="huc4")

sel <- c(!(huc$STATES %in% c("AK", "AS", "AK,CN", "HI", 
    "PR", "GU", "MP", "VI")))
huc <- huc[sel, ]
simple <- gSimplify(huc, tol=0.01, topologyPreserve=TRUE)

plot <- fortify(huc)
huc$id <- row.names(huc)
plot <- left_join(plot, huc@data)
huc <- plot

neighborcol <- brewer.pal(8, "Blues")[4]
focalcol <-brewer.pal(8, "Blues")[8]

#central <- plot.huc[plot.huc$id == "195", ]
#centralname <- as.character(plot.huc$HUC4[plot.huc$id == 195][1])

central <- plot.huc[plot.huc$id == "200", ] #1706
centralname <- as.character(plot.huc$HUC4[plot.huc$id == 200][1])


#############################################################################
#panel one - queens
#############################################################################
#calculate queens matrix
W <- define.neighborhood("queens")

#isolate data for the central huc
neighbors <- data.frame(huc = rownames(W), 
    weight = W[rownames(W)== centralname, ])

huc$weight <- neighbors$weight[match(huc$HUC4, neighbors$huc)]
plot.huc <- huc

#make a plot
queens <- ggplot() +
    geom_polygon(data = plot.huc, aes(x = long, y = lat, group = group,
        fill = as.factor(weight)), colour = "black") +
    scale_fill_manual(values = c("white", neighborcol), guide = FALSE) +
    geom_polygon(data = central, aes(x = long, y = lat), 
                 fill=focalcol, color="black") +
    coord_fixed(1.3) +
    theme(panel.background=element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank()) +
    ggtitle("Queens") #+

#save plot
#save.location <- "~/Honors Thesis/Thesis Work/Paper Figures/"
save.location <- "~/Documents/Avian_Influenza/HP/"
jpeg(paste(save.location, "Figure 3a.jpeg", sep = ""))
queens
dev.off()

###########################################################################################
#panel two - weighted queens
#need to fix the coloring
###########################################################################################

#calculate weighted queens matrix
W <- define.neighborhood("weightedqueens")

#isolate data for the central huc (0301)
#neighbors <- data.frame(huc = rownames(W), weight = W[195, ])
neighbors <- data.frame(huc = rownames(W), 
    weight = W[rownames(W)== centralname, ])

#reread in huc data in original form
# huc <- readOGR(dsn = "hydrologic_units", layer="huc4")
# sel <- !(huc$STATES %in% c(
#     "AK", "AS", "AK,CN", "HI", "PR", "GU", "MP", "VI"))
# 
# # prep for plotting (subset and simplify)
# temp <- gSimplify(huc, tol = 0.001, topologyPreserve = TRUE)
# plot.huc <- SpatialPolygonsDataFrame(temp, data = huc@data)
# plot.huc <- plot.huc[sel, ]
# 
# # prep for ggplot2 and add the data back
# plot <- fortify(plot.huc)
# plot.huc$id <- row.names(plot.huc)
# plot.huc <- left_join(plot, plot.huc@data) # requires library dplyr

#####################################
# NOTE: This hashed out bit below writes over the weights added,
# should be removed.  I added simplify (lowers the reolution) above. 
#####################################
#convert data for mapping
#huc <- spTransform(huc, CRS("+proj=longlat +datum=WGS84"))
#simple <- gSimplify(huc, tol = 0.01, topologyPreserve = TRUE)
#huc <- fortify(simple)

# save them for yourself
pal <- brewer.pal(8, "Blues")


# Use this one: specify your own limits
################################################
# This way is nice if you want consistent colors ranges between plots-
# so blue type 1 in plot 1 means the same in plot 2, specify the ranges yourself

# Figure out what the 0, 20th, 40th, 60th, 80th and 100th quantiles are. 
# All 0s set to no association
#match weighted data to huc spatial data
plot.huc$weight <- NA
plot.huc$weight <- neighbors$weight[match(plot.huc$HUC4, neighbors$huc)]

plot.huc$scaledwt <- plot.huc$weight/max(plot.huc$weight)
qa <- quantile(plot.huc$scaledwt[plot.huc$scaledwt > 0], c(0, 0.33, 0.66, 1))
qa

plot.huc$nquantile <- cut(plot.huc$scaledwt, c(0, qa), 
    labels = c("None", "0-33%", "33-66%", "66-100%"),
    include.lowest = TRUE, right = FALSE)

# specify palette with 5 blues and white if no weight
pal <-c("white", brewer.pal(8, "Blues")[3:6])
focalcol <-brewer.pal(8, "Blues")[8]

weightedqueens <- ggplot() +
    geom_polygon(data = plot.huc, aes(x = long, y = lat, group = group,
        fill = nquantile), colour = "black") +
    scale_fill_manual(values = pal, guide = FALSE) +#delete guide=FALSE for legend
    geom_polygon(data = central, aes(x = long, y = lat), 
                 fill=focalcol, color="black") +
    coord_fixed(1.3) +
    theme(panel.background=element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank()) +
    ggtitle("Weighted Queens") #+
    #labs(fill = "")  # to get rid of or specify a differnet label


save.location <- "~/Documents/Avian_Influenza/HP/"
jpeg(paste(save.location, "Figure 3b.jpeg", sep = ""))
weightedqueens
dev.off()


# END ERIN'S UPDATES

#########################################################################################
# panel four - network
#########################################################################################

#calculate the network W matrix
W <- define.neighborhood("network")

#isolate data for the central huc
neighbors <- data.frame(huc=rownames(W), weight = W[rownames(W)==centralname, ])

#match weighted data to huc spatial data
huc$weight <- NA
huc$weight <- neighbors$weight[match(huc4$HUC4, neighbors$huc)]

#make a plot
ggplot() +
  geom_polygon(data=usa, aes(x=long, y=lat, group=group), fill=NA, color="white") +
  geom_polygon(data=huc, aes(x=long, y=lat, group=group, color=weight), 
               fill=brewer.pal(8, "Blues"), color="black") +
  #scale_fill_brewer(type="seq", palette="Blues") +
  #geom_polygon(data=central, aes(x=long, y=lat), fill="darkblue", color="black")+
  coord_fixed(1.3) +
  theme(panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank()) +
  ggtitle("Network")

#save plot
jpeg("~/Honors Thesis/Thesis Work/Paper Figures/Figure 3d.jpeg")
network
dev.off()





############################################################################################
#panel three - knn3 (randomly picks three neighbors)
############################################################################################
#will make four small plots to show different combinations

neighbors <- matrix(data=c("34", "197", "36", "180", "182", "202"), nrow=1, ncol=6)

#quadrant 1 data
sample1 <- sample(neighbors, 3, replace=F)
neighbor1 <- huc[huc$id == sample1[1], ]
neighbor2 <- huc[huc$id == sample1[2], ]
neighbor3 <- huc[huc$id == sample1[3], ]
#quadrant 1 plot
version1 <- ggplot()+
    geom_polygon(data=usa, aes(x=long, y=lat, group=group), fill=NA, color="white") +
    geom_polygon(data=huc, aes(x=long, y=lat, group=group), fill=NA, color="black") +
    geom_polygon(data=central, aes(x=long, y=lat, group=group), fill="darkblue", color="black")+
    geom_polygon(data=neighbor1, aes(x=long, y=lat), fill="skyblue", color="black") +
    geom_polygon(data=neighbor2, aes(x=long, y=lat), fill="skyblue", color="black") +
    geom_polygon(data=neighbor3, aes(x=long, y=lat), fill="skyblue", color="black") +
    coord_fixed(1.3) +
    theme(panel.background=element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank()) +
    ggtitle("Random Variation 1")

#quadrant 2 data
sample2 <- sample(neighbors, 3, replace=F)
neighbor1 <- huc[huc$id == sample2[1], ]
neighbor2 <- huc[huc$id == sample2[2], ]
neighbor3 <- huc[huc$id == sample2[3], ]
#quadrant 2 plot
version2 <- ggplot()+
    geom_polygon(data=usa, aes(x=long, y=lat, group=group), fill=NA, color="white") +
    geom_polygon(data=huc, aes(x=long, y=lat, group=group), fill=NA, color="black") +
    geom_polygon(data=central, aes(x=long, y=lat, group=group), fill="darkblue", color="black")+
    geom_polygon(data=neighbor1, aes(x=long, y=lat), fill="skyblue", color="black") +
    geom_polygon(data=neighbor2, aes(x=long, y=lat), fill="skyblue", color="black") +
    geom_polygon(data=neighbor3, aes(x=long, y=lat), fill="skyblue", color="black") +
    coord_fixed(1.3) +
    theme(panel.background=element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank()) +
    ggtitle("Random Variation 2")

#quadrant 3 data
sample3 <- sample(neigbhors, 3, replace = TRUE) #? 
neighbor1 <- huc[huc$id == sample3[1], ]
neighbor2 <- huc[huc$id == sample3[2], ]
neighbor3 <- huc[huc$id == sample3[3], ]
#quadrant 3 plot
version3 <- ggplot()+
    geom_polygon(data=usa, aes(x=long, y=lat, group=group), fill=NA, color="white") +
    geom_polygon(data=huc, aes(x=long, y=lat, group=group), fill=NA, color="black") +
    geom_polygon(data=central, aes(x=long, y=lat, group=group), fill="darkblue", color="black")+
    geom_polygon(data=neighbor1, aes(x=long, y=lat), fill="skyblue", color="black") +
    geom_polygon(data=neighbor2, aes(x=long, y=lat), fill="skyblue", color="black") +
    geom_polygon(data=neighbor3, aes(x=long, y=lat), fill="skyblue", color="black") +
    coord_fixed(1.3) +
    theme(panel.background=element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank()) +
    ggtitle("Random Variation 3")

#quadrant 4 data
sample4 <- sample(neighbors, 3, replace=F)
neighbor1 <- huc[huc$id == sample4[1], ]
neighbor2 <- huc[huc$id == sample4[2], ]
neighbor3 <- huc[huc$id == sample4[3], ]
#quadrant 4 plot
version4 <- ggplot()+
    geom_polygon(data=usa, aes(x=long, y=lat, group=group), fill=NA, color="white") +
    geom_polygon(data=huc, aes(x=long, y=lat, group=group), fill=NA, color="black") +
    geom_polygon(data=central, aes(x=long, y=lat, group=group), fill="darkblue", color="black")+
    geom_polygon(data=neighbor1, aes(x=long, y=lat), fill="skyblue", color="black") +
    geom_polygon(data=neighbor2, aes(x=long, y=lat), fill="skyblue", color="black") +
    geom_polygon(data=neighbor3, aes(x=long, y=lat), fill="skyblue", color="black") +
    coord_fixed(1.3) +
    theme(panel.background=element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank()) +
    ggtitle("Random Variation 4")

#make and save a single plot with all four variations
jpeg("~/Honors Thesis/Thesis Work/Paper Figures/Figure 3c.jpeg")
knn3 <- multiplot(version1, version2, version3, version4, cols=2)
dev.off()
