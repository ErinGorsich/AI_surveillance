# Plot priors
se <- rbeta(10000, shape1 = 20.833, shape2 = 4.148)
sp <- rbeta(10000, shape1 = 8.403, shape2 = 1.001)
par(mfrow = c(1, 2))
hist(se); hist(sp)
quantile(se, c(0.025, 0.5, 0.975))
quantile(sp, c(0.025, 0.5, 0.975))

# if I were choosing priors...
# Se = median at 86.6
# 95% sure the median it is > 60
# Spe = median at 99; 95% sure median is > 60


# Plot a maps with the networks shown
x <- c("rgdal", "rgeos", "maptools", "dplyr", "sp", "maps", 
       "ggmap", "ggplot2", "plyr", "tmap", "RColorBrewer") 
lapply(x, library, character.only = TRUE)

source('define.neighborhood.R', chdir = TRUE)
W <- define.neighborhood(method = "network")

#prepare huc spatial data
huc <- readOGR(dsn = "hydrologic_units", layer="huc4")
sel <- c(!(huc$STATES %in% c("AK", "AS", "AK,CN", "HI", 
                             "PR", "GU", "MP", "VI")))
huc <- huc[sel, ]

plot <- fortify(huc)
huc$id <- row.names(huc)
plot <- left_join(plot, huc@data)

neighborcol <- brewer.pal(8, "Blues")[4]
focalcol <-brewer.pal(8, "Blues")[8]

central <- huc[huc$id == "200", ] #1706
centralname <- as.character(huc$HUC4[huc$id == 200][1])

#isolate data for the central huc
neighbors <- data.frame(huc = rownames(W), 
    weight = W[rownames(W) == centralname, ])

huc <- plot
huc$weight <- neighbors$weight[match(huc$HUC4, neighbors$huc)]


#make a plot
plot.huc <- huc
net <- ggplot() +
    geom_polygon(data = plot.huc, aes(x = long, y = lat, group = group,
                                      fill = as.factor(weight)), colour = "black") +
    scale_fill_manual(values = c("white", neighborcol), guide = FALSE) +
    geom_polygon(data = central, aes(x = long, y = lat), 
                 fill=focalcol, color="black") +
    coord_fixed(1.3) +
    theme(panel.background=element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank()) +
    ggtitle("Network") #+

#save plot
#save.location <- "~/Honors Thesis/Thesis Work/Paper Figures/"
save.location <- "~/Documents/Avian_Influenza/HP/"
jpeg(paste(save.location, "Figure 3 network.jpeg", sep = ""))
net
dev.off()

# Plot a maps for seasonal data 
t <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_n_mall_all.rds")
t2 <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_y_mall_all.rds")
setwd("~/GitHub/HP_AI_Surveillance")
huc <- readOGR(dsn = "hydrologic_units", layer="huc4")
temp <- gSimplify(huc, tol = 0.01, topologyPreserve = TRUE)
temp <- SpatialPolygonsDataFrame(temp, data = huc@data)
sel <- !(temp$STATES %in% c("AK", "AS", "AK,CN", "HI", "PR", "GU", "MP", "VI"))
huc <- temp[sel, ]

subsites <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/subset_sites.rds")
save.location <- "~/Documents/Avian_Influenza/HP/"

make_maps = function (t, month, year, huc) {
    m <- t[month, year, ]
    huc$weight <- NA
    huc$weight[match(huc$HUC4, subsites, nomatch = 0)] <- m
    huc$weight[is.na(huc$weight)] <- 0
    #if (max(huc$weight) > 0 & max(huc$weight) < 6) {
    #    huc$weightcat <- cut(huc$weight, c(0, 1, max(huc$weight)+1), 
    #        include.lowest = TRUE, right = FALSE)
    #    col <- c("white", neighborcol)
    #} else if (max(huc$weight) == 0) {
    #    huc$weightcat <- rep(0, length(huc$weight))
    #    col <- c("white")
    #} else {
    #    huc$weightcat <- cut(huc$weight, c(0, 1, 6, max(huc$weight)+1), 
    #        include.lowest = TRUE, right = FALSE)
    #    col <- c("white", neighborcol, focalcol)
    #}
    
    #plot <- fortify(huc)
    #huc$id <- row.names(huc)
    #plot <- left_join(plot, huc@data)
    jpeg(paste(save.location, "numsamples_", month, "_", year, sep = ""))
    qtm(huc, "weight")
    #ggplot() +
    #    geom_polygon(data = plot, aes(x = long, y = lat, group = group,
    #        fill = weightcat), colour = "black") +
    #    #scale_fill_manual(values = col, guides = FALSE)
    #    coord_fixed(1.3) +
    #    theme(panel.background=element_blank(), panel.grid.major = element_blank(),
    #        panel.border = element_blank(), axis.title = element_blank(),
    #       axis.ticks = element_blank(), axis.text = element_blank())
    dev.off()
}

for (i in 1:12){
    m <- t[i, 1, ]
    huc$weight <- NA
    huc$weight[match(huc$HUC4, subsites, nomatch = 0)] <- m
    huc$weight[is.na(huc$weight)] <- 0
    huc$weightcat <- 2
    huc$weightcat[huc$weight < 6] <- 1
    huc$weightcat[huc$weight == 0] <- 0
    
    jpeg(paste(save.location, "numsamples_", i, "_", 1, sep = ""))
    qtm(huc, "weightcat") + tm_legend(show=FALSE)
    dev.off()
}

# plot aggregate month accross years
i <- 1
for (i in 1:12){
    m <- apply(t[i, , ], 2, sum)
    huc$weight <- NA
    huc$weight[match(huc$HUC4, subsites, nomatch = 0)] <- m
    huc$weight[is.na(huc$weight)] <- 0
    huc$weightcat <- 3
    huc$weightcat[huc$weight <50] <- 2
    huc$weightcat[huc$weight < 6] <- 1
    huc$weightcat[huc$weight == 0] <- 0
    jpeg(paste(save.location, "multiyear_aggregates_", i, "_", 1, sep = ""))
    qtm(huc, "weightcat") + tm_legend(show=FALSE)
    dev.off()
}


# same for number of positive sampels aggregated accross years
i <- 1
for (i in 1:12){
    m <-  apply(t2[i, , ], 2, sum)/apply(t[i, , ], 2, sum)
    huc$weight <- NA
    huc$weight[match(huc$HUC4, subsites, nomatch = 0)] <- m
    huc$weight[is.na(huc$weight)] <- 0
    huc$weightcat <- 3
    huc$weightcat[huc$weight < 0.5] <- 2
    huc$weightcat[huc$weight < 0.25] <- 1
    huc$weightcat[huc$weight == 0] <- 0
    jpeg(paste(save.location, "multiyear_aggregates_prev", i, "_", 1, sep = ""))
    qtm(huc, "weightcat") + tm_legend(show=FALSE)
    dev.off()
}

