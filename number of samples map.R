#makes a map about number of samples per huc
#dabbling ducks, diving ducks, geese and swans

#load packages
library(rgdal)
library(rgeos)
library(maptools)
library(dplyr)
library(sp)
library(ggmap)
library(ggplot2)
library(plyr)
library(RColorBrewer)

#read in spatial data
setwd("~/Honors Thesis/Project/hydrologic_units")
huc4 <- readOGR(dsn=".", layer = "huc4")

#read in data
data <- readRDS("~/Github/AVHS_samplingevent_locationname.rds")
data <- SpatialPointsDataFrame(coords = cbind(data$long, data$lat),
                               proj4string = huc4@proj4string, data = data)
data <- data[data$species.group %in% c("1", "2", "4"), ]

#select data in wanted hucs
sel <- !(huc4$STATES %in% c("AK", "AS", "AK,CN", "HI", "PR", "GU", "MP", "VI"))
huc4 <- huc4[sel, ]

#calculate number of samples per watershed
test <- over(data, huc4)
allsamples <- as.data.frame(table(test$HUC4))
allsamples$Var1 <- as.character(allsamples$Var1)
huc4$HUC4 <- as.character(huc4$HUC4)
huc4$NUM_SAMPLES <- allsamples$Freq[match(huc4$HUC4, allsamples$Var1)]

#simplify data
temp <- gSimplify(huc4, tol=0.01, topologyPreserve = TRUE)
huc4.plot <- SpatialPolygonsDataFrame(temp, data=huc4@data)
huc4.plot <- huc4.plot[sel, ]
plot <- fortify(huc4.plot)
huc4.plot$id <- row.names(huc4.plot)
huc4.plot <- left_join(plot, huc4.plot@data)

#group data into bins
qa <- quantile(huc4.plot$NUM_SAMPLES, c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE)
qa
huc4.plot$quantile <- cut(huc4.plot$NUM_SAMPLES, c(qa[1], 1, qa[2], qa[3], qa[4], qa[5], qa[6]),
                          labels = c("No Samples", "1-108", "109-389", "390-934", "935-1588", "> 1588"), 
                          include.lowest = TRUE)
pal <- brewer.pal(7, "Blues")
better <- c("white", pal[2:6])
#make a map
jpeg("~/Honors Thesis/Thesis Work/Figures/numberofsamples_esa.jpeg")
ggplot() +
  geom_polygon(data=huc4.plot, aes(x=long, y=lat, group=group, fill=quantile), colour="black") +
  scale_fill_manual(values = better) +
  coord_fixed(1.3) +
  theme(panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size=10)) +
  ggtitle("Number of Samples by Watershed")
dev.off()
