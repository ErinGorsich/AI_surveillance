library('rgdal')
library('maps')
library('raster') 
library('plyr')
library('rgeos')
library('rjags')

# Data organization, HP's model
######################################
#ai <- readRDS("~/GitHub/AvianInfluenza/AVHS_v2.rds")

# for full dataset: 
ai <- readRDS("/Users/gorsich/Documents/Avian_Influenza/Data/HP_dataset.rds")
ai <- ai[!(ai$state %in% c("AK", "HI")), ]

######################################
# Add columns to tell us the watershed
######################################
setwd("~/Documents/Avian_Influenza/Data/hydrologic_units")
huc4 <- shapefile("huc4.shp")
projection(huc4) <- CRS("+proj=longlat +ellps=WGS84")

pt <- ai[ , c("lat", "long")]
coordinates(pt) <- ~ long + lat
proj4string(pt) <- CRS("+proj=longlat +ellps=WGS84")

# really crazy slow, maybe save output to avoid this step more than once
ai$huc4 <- extract(huc4, pt)$HUC4
#ai$huc4name <- extract(huc4, pt)$NAME
table(ai$huc4)

# remove poorly sampled years (now only 2007-2010; 2015-2017)
aired <- ai[!(ai$collection.year %in% c("2011", "2014", "2016", "2017")), ]
ai <- aired
ai <- ai[!is.na(ai$huc4), ]

# remove strange hucs
sel <- !(huc4$STATES %in% c("AK", "AS", "AK,CN", "HI", "PR", "GU", "MP", "VI"))
huc4 <- huc4[sel, ]

# order we want the hucs in!
hucs <- huc4@data$HUC4

######################################
# define seasonal data
######################################
nsites <- length(unique(ai$huc4)) # 181
nseasons <- 3
nyears <- length(unique(ai$collection.year))  # 5

# define an empty array
y <- array(NA, dim = c(nseasons, nyears, nsites))
n <- array(NA, dim = c(nseasons, nyears, nsites))

# order months/years chronologically; order huc4 by hucs
sites <- unique(ai$huc4)
sites <- hucs[hucs %in% sites]
seasons <- c("fall", "winter", "springsummer")
seasonlist <- list(
    fall = c("08", "09", "10"), winter = c("11", "12", "01"), 
    springsummer = c("02", "03", "04", "05", "06", "07"))
years <- unique(ai$collection.year)

# save common watersheds
saveRDS(sites, "~/Documents/Avian_Influenza/HP/statistical model/subset_sites.rds")

# fill array
######################################
# consider saving as well 
for (i in 1:nsites){
    for (j in 1:nseasons){
        for (k in 1:nyears){
            n[j, k, i] <- length(ai[ai$huc4 == sites[i] & 
                ai$collection.month %in% seasonlist[[j]] & 
                ai$collection.year == years[k], 1])
            y[j, k, i] <- length(ai[ai$huc4 == sites[i] & 
                ai$collection.month %in% seasonlist[[j]] & 
                ai$collection.year == years[k] & 
                ai$AIpcr_susneg == "positive", 1])
        }
    }
}

# to see it: 
n[1:3, 1:5, 1]
y[1:3, 1:5, 1]

saveRDS(n, "~/Documents/Avian_Influenza/HP/statistical model/data_n_mall3season.rds")
saveRDS(y, "~/Documents/Avian_Influenza/HP/statistical model/data_y_mall3season.rds")

nred <- n[1:2, , ]
yred <- y[1:2, , ]
saveRDS(n, "~/Documents/Avian_Influenza/HP/statistical model/data_n_mall2season.rds")
saveRDS(y, "~/Documents/Avian_Influenza/HP/statistical model/data_y_mall2season.rds")


######################################
# define data
######################################
nsites <- length(unique(ai$huc4)) # 181
nmonths <- 12
nquarters <- 4
nyears <- length(unique(ai$collection.year))  # 5

# define an empty array
y <- array(NA, dim = c(nmonths, nyears, nsites))
n <- array(NA, dim = c(nmonths, nyears, nsites))

# order months/years chronologically; order huc4 by hucs
sites <- unique(ai$huc4)
sites <- hucs[hucs %in% sites]
months <- as.character(unique(ai$collection.month))
months <- months[c(4, 9, 10, 12, 11, 8, 6, 2, 7, 5, 1, 3)]
months
years <- unique(ai$collection.year); years
#years <- years[c(2,1, 3, 4, 5)]

# save common watersheds
saveRDS(sites, "~/Documents/Avian_Influenza/HP/statistical model/subset_sites.rds")

# fill array
######################################
# consider saving as well 
for (i in 1:nsites){
	for (j in 1:nmonths){
		for (k in 1:nyears){
			n[j, k, i] <- length(ai[ai$huc4 == sites[i] & 
			     ai$collection.month == months[j] & 
			     ai$collection.year == years[k], 1])
			y[j, k, i] <- length(ai[ai$huc4 == sites[i] & 
			     ai$collection.month == months[j] & 
			     ai$collection.year == years[k] & 
			     ai$AIpcr_susneg == "positive", 1])
		}
	}
}

# to see it: 
n[1:12, 1:5, 1:5]

# if larger, mallard only dataset, only hucs with data
saveRDS(n, "~/Documents/Avian_Influenza/HP/statistical model/data_n_mall_all.rds")
saveRDS(y, "~/Documents/Avian_Influenza/HP/statistical model/data_y_mall_all.rds")

######################################
# define data - if want all watersheds!
######################################
# nsites <- length(huc4@polygons)  # n = 202 now not 136
# sites <- huc4@data$HUC4
# nmonths <- 12
# nyears <- length(unique(ai$collection.year))  # 6
# saveRDS(sites, "~/Documents/Avian_Influenza/HP/statistical model/sites.rds")
# 
# # define an empty array
# y <- array(NA, dim = c(nmonths, nyears, nsites))
# n <- array(NA, dim = c(nmonths, nyears, nsites))
# 
# months <- as.character(unique(ai$collection.month))
# months <- months[c(4, 9, 10, 12, 11, 8, 6, 2, 7, 5, 1, 3)]
# months
# years <- unique(ai$collection.year)
# years
# 
# # fill array
# ######################################
# # consider saving as well
# for (i in 1:nsites){
#     for (j in 1:nmonths){
#         for (k in 1:nyears){
#             n[j, k, i] <- length(ai[ai$huc4 == sites[i] &
#                 ai$collection.month == months[j] &
#                 ai$collection.year == years[k], 1])
#             y[j, k, i] <- length(ai[ai$huc4 == sites[i] &
#                 ai$collection.month == months[j] &
#                 ai$collection.year == years[k] &
#                 ai$AIpcr_susneg == "positive", 1])
#         }
#     }
# }
# 
# # to see it:
# n[1:12, 1:5, 1:5]
# 
# # if larger, mallard only dataset with all hucs
# saveRDS(n, "~/Documents/Avian_Influenza/HP/statistical model/data_n_mall_all_fullhuc.rds")
# saveRDS(y, "~/Documents/Avian_Influenza/HP/statistical model/data_y_mall_all_fullhuc.rds")
# 


