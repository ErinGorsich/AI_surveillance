# Bird movement network preparation

setwd("~/Documents/Avian_Influenza/Data/BandData/Band_Encounters_withdates")
mall <- read.csv("bands_10.csv")
backup <- mall

# groom banding data for only hunter harvest recoveries within 1 year of banding
mall <- mall[!is.na(mall$GISRLong), ]
mall <- mall[!is.na(mall$GISRLat), ]
mall <- mall[!is.na(mall$r.date), ]
mall <- mall[!is.na(mall$d.date), ]
mall$r.date <- as.character(mall$r.date)
mall$b.date <- as.character(mall$b.date)
mall$r.date <- as.Date(mall$r.date, format = "%Y-%m-%d")
mall$b.date <- as.Date(mall$b.date, format = "%Y-%m-%d")
mall <- mall[!is.na(mall$r.date), ]
mall <- mall[!is.na(mall$b.date), ]
mall$diff <- mall$r.date - mall$b.date

# if recaptured in a year or less
backup <- mall[mall$diff < 365 & mall$diff > 60, ]
mall <- backup

# remove if recovery accuracy is too large
mall <- mall[!(mall$R_Coord_Location_Accuracy_Desc %in% c("Country", 
    "Northwest of this point", "State")), ]

# dates within sampling window
mall <- mall[mall$r.year > 2006 & mall$r.year < 2016, ] # 66689

# only hunter harvested - no bird banders
mall <- mall[mall$Who_VWho == "Finder", ]

# add o.huc and d.huc to the gps locations
##########################################################
ai <- mall
setwd("~/Documents/Avian_Influenza/Data/hydrologic_units")
huc4 <- shapefile("huc4.shp")
sel <- !(huc4$STATES %in% c("AK", "AS", "AK,CN", "HI", "PR", "GU", "MP", "VI"))
huc4 <- huc4[sel, ]
projection(huc4) <- CRS("+proj=longlat +ellps=WGS84")

bpt <- ai[ , c("GISBLat", "GISBLong")]
coordinates(bpt) <- ~ GISBLong + GISBLat
proj4string(bpt) <- CRS("+proj=longlat +ellps=WGS84")

rpt <- ai[ , c("GISRLat", "GISRLong")]
coordinates(rpt) <- ~ GISRLong + GISRLat
proj4string(rpt) <- CRS("+proj=longlat +ellps=WGS84")

# really crazy slow, maybe save output to avoid this step more than once
ai$o.huc4 <- extract(huc4, bpt)$HUC4
table(ai$o.huc4)

ai$d.huc4 <- extract(huc4, rpt)$HUC4
table(ai$d.huc4)

# remove poorly sampled years (now only 2007-2010; 2015-2017)
backup <- ai; length(ai[,1])
ai <- ai[!is.na(ai$o.huc4), ]
ai <- ai[!is.na(ai$d.huc4), ]

birdmovement <- ai[ , c('o.huc4', 'd.huc4')]
# ? 
#ai <- ai[!(ai$state %in% c("AK", "HI")), ]
saveRDS(birdmovement, "birdmovementnets.rds")
