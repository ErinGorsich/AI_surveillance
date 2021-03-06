#cleans data and assigns sampling events and species groups

library(lubridate)
library(stringr)
library(dplyr)
library(rgdal)
library(raster)

##################################################################################
#data read in
#################################################################################
setwd("~/HP/Data")
data <- read.csv("FluA_FY07toFY16_Webb.csv")

colnames(data) <- c("subjectID", "barcode", "band", "species.code.full", 
                   "sex", "age", "collection.group", "lat", "long", "state", "county", "sample.type", 
                   "collection.date.init", "agency", "laboratory", "collection.method", 
                   "livedead", "submitted.date", "tested.date", "regional_AIpcr", "regional_H5pcr",
                   "regional_H7pcr",	"NVSL_AIpcr", "NVSL_H5pcr", "NVSL_H7pcr",
                   "virus.isolation", "subtype", "highlow", "ivpi", "biological", "watershed" )

#####################################################################################
#age and sex
#####################################################################################

#make entries in the sex and age columns characters
data$sex <- as.character(data$sex)
data$age <- as.character(data$age)

#rename variables in the sex column
data$sex[data$sex == "Male"] <- "male"
data$sex[data$sex == "Female"] <- "female"
data$sex[data$sex == "Undetermined"] <- "unkown"

#rename variables in the age column
data$age[data$age == "After Hatch Year"] <- "AHY"
data$age[data$age == "Hatch Year"] <- "HY"
data$age[data$age == "Undetermined"] <- "U"

#shortens species name to only the code
data$species.code.full <- as.character(data$species.code.full)
data$species.code.full <- substr(data$species.code.full, 1,4)
  #substr subsets portions of a character vector
  #this command asks for the first four characters of each character vector
  #leaves only the species code

#####################################################################################
#Dates: change formatting to be consistent
#####################################################################################

data$collection.date.init <- as.character(data$collection.date.init)
data$collection.date.format <- NA
data$collection.month <- NA
data$collection.year <- NA
data$collection.day <- NA

#assigns the date as being formatted as month/day/year or day/month/year
dmy <- grep('[A-Z]', data$collection.date.init, value = TRUE)
data$collection.date.format[which(data$collection.date.init %in% dmy)] <- "dmy"
data$collection.date.format[is.na(data$collection.date.format)] <- "mdy"

#seperates dates into columns for day, month, and year
for (i in 1:length(data$collection.date.init)) {
  if (data$collection.date.format[i] == "mdy"){
    temp <- strsplit(data$collection.date.init[i], "/")[[1]]
    ifelse(nchar(temp[1]) == 1, 
           data$collection.month[i] <- paste("0", temp[1], sep = ""),
           data$collection.month[i] <- temp[1])
    data$collection.year[i] <- temp[3]
    ifelse(nchar(temp[2]) == 1, 
           data$collection.day[i] <- paste("0", temp[2], sep = ""), 
           data$collection.day[i] <- temp[2])
  } 
  if (data$collection.date.format[i] == "dmy"){
    temp <- strsplit(data$collection.date.init[i], "-")[[1]]
    ifelse(nchar(temp[1]) == 1, 
           data$collection.day[i] <- paste("0", temp[1], sep = ""),
           data$collection.day[i] <- temp[1] )
    data$collection.year[i] <- paste("20", temp[3], sep = "")
    ifelse(nchar(temp[2]) == 1, 
           data$collection.month[i] <- paste("0", temp[2], sep = ""),
           data$collection.month[i] <- temp[2]
    )
  }
}

#converts months from words to numbers
data$collection.month[data$collection.month == "Aug"] <- "08"
data$collection.month[data$collection.month == "Dec"] <- "12"
data$collection.month[data$collection.month == "Feb"] <- "02"
data$collection.month[data$collection.month == "Jan"] <- "01"
data$collection.month[data$collection.month == "Jul"] <- "07"
data$collection.month[data$collection.month == "Jun"] <- "06"
data$collection.month[data$collection.month == "May"] <- "05"
data$collection.month[data$collection.month == "Nov"] <- "11"
data$collection.month[data$collection.month == "Oct"] <- "10"
data$collection.month[data$collection.month == "Sep"] <- "09"

#creates a column with the complete date in the number form
data$collection.date.char <- paste(data$collection.month, data$collection.day, data$collection.year, sep = "-")
#creates a column to put the dates in year, month, day form
data$collection.date <- mdy(data$collection.date.char)

backup.dates <- data
######################################################################################
#Sample Type
######################################################################################

# remove sentinal species
OC <- data[data$sample.type == "Oral + Cloacal", ]
CS <- data[data$sample.type == "Cloacal Swab", ]
OP <- data[data$sampl.type == "Oropharyngeal", ]
TS <- data[data$sample.type == "Tracheal Swab", ]

#remove samples that are not of Oral+Cloacal type
data <- data[data$sample.type == "Oral + Cloacal", ]

backup.type <- data

#####################################################################################
#Disease Columns
#####################################################################################

#ultimately, test results follow the national lab
data$AIpcr <- "NoResult"
data$AIpcr[data$regional_AIpcr == "ND"] <- "NEG"
data$AIpcr[data$regional_AIpcr == "SUS"] <- "SUS"
data$AIpcr[data$regional_AIpcr == "DET" & data$NVSL_AIpcr == ""] <- "POS"
data$AIpcr[data$regional_AIpcr == "DET" & data$NVSL_AIpcr == "DET"] <- "POS"
data$AIpcr[data$regional_AIpcr == "DET" & data$NVSL_AIpcr == "ND"] <- "NEG"
data$AIpcr[data$regional_AIpcr == "" & data$NVSL_AIpcr == "ND"] <- "NEG"
data$AIpcr[data$regional_AIpcr == "" & data$NVSL_AIpcr == "DET"] <- "POS"

#suspected results are considered to be negative
data$AIpcr_susneg <- "NoResult"
data$AIpcr_susneg[data$AIpcr == "NEG"] <- "negative"
data$AIpcr_susneg[data$AIpcr == "SUS"] <- "negative"
data$AIpcr_susneg[data$AIpcr == "POS"] <- "positive"

#H5 results are ultimately determined by the national lab, but default to regional lab
data$H5pcr <- "No Result"
data$H5pcr[data$regional_H5pcr == "ND" & data$NVSL_H5pcr == ""] <- "NEG"
data$H5pcr[data$regional_H5pcr == "ND" & data$NVSL_H5pcr == "DET"] <- "POS"
data$H5pcr[data$regional_H5pcr == "ND" & data$NVSL_H5pcr == "ND"] <- "NEG"
data$H5pcr[data$regional_H5pcr == "SUS" & data$NVSL_H5pcr == "DET"] <- "POS"
data$H5pcr[data$regional_H5pcr == "SUS" & data$NVSL_H5pcr == "ND"] <- "NEG"
data$H5pcr[data$regional_H5pcr == "DET" & data$NVSL_H5pcr == ""] <- "POS"
data$H5pcr[data$regional_H5pcr == "DET" & data$NVSL_H5pcr == "DET"] <- "POS"
data$H5pcr[data$regional_H5pcr == "DET" & data$NVSL_H5pcr == "ND"] <- "NEG"
data$H5pcr[data$regional_H5pcr == "" & data$NVSL_H5pcr == "DET"] <- "POS"
data$H5pcr[data$regional_H5pcr == "" & data$NVSL_H5pcr == "ND"] <- "NEG"
data$H5pcr[data$AIpcr_susneg == "negative" & data$H5pcr == "NoResult"] <- "NEG"

#H7 results are ultimately determined by the national lab, but default to the regional lab
data$H7pcr <- "No Result"
data$H7pcr[data$regional_H7pcr == "ND" & data$NVSL_H7pcr == ""] <- "NEG"
data$H7pcr[data$regional_H7pcr == "ND" & data$NVSL_H7pcr == "DET"] <- "POS"
data$H7pcr[data$regional_H7pcr == "ND" & data$NVSL_H7pcr == "ND"] <- "NEG"
data$H7pcr[data$regional_H7pcr == "SUS" & data$NVSL_H7pcr == "DET"] <- "POS"
data$H7pcr[data$regional_H7pcr == "SUS" & data$NVSL_H7pcr == "ND"] <- "NEG"
data$H7pcr[data$regional_H7pcr == "DET" & data$NVSL_H7pcr == ""] <- "POS"
data$H7pcr[data$regional_H7pcr == "DET" & data$NVSL_H7pcr == "DET"] <- "POS"
data$H7pcr[data$regional_H7pcr == "DET" & data$NVSL_H7pcr == "ND"] <- "NEG"
data$H7pcr[data$regional_H7pcr == "" & data$NVSL_H7pcr == "DET"] <- "POS"
data$H7pcr[data$regional_H7pcr == "" & data$NVSL_H7pcr == "ND"] <- "NEG"
data$H7pcr[data$AIpcr_susneg == "negative" & data$H7pcr == "NoResult"] <- "NEG"

#change no result to unknown
data$AIpcr_susneg[data$AIpcr_susneg == "NoResult"] <- "UNKN"
data$H5pcr[data$AIpcr == "NoResult"] <- "UNKN"
data$H7pcr[data$H7pcr == "NoResult"] <- "UNKN"

#change characters into factors
data$AIpcr_susneg <- as.factor(data$AIpcr_susneg)
data$H5pcr <- as.factor(data$H5pcr)
data$H7pcr <- as.factor(data$H7pcr)

backup.test <- data

######################################################################################
#remove years of data with low sampling numbers and 2016
#####################################################################################
data <- data[!(data$collection.year %in% c("2011", "2014", "2016", "2017")), ]

backup.year <- data

# setwd("~/Github")
setwd("~/HP/Data")
saveRDS(data, "AVHS_samplingevent.rds")


###########################################################################################################################
#Add Species Group
###########################################################################################################################
#assigns a species group number to each sample based on species
  #1 <- dabbling ducks (including perching ducks)
  #2 <- diving ducks
  #3 <- sea ducks
  #4 <- aserines (geese and swans)
  #5 <- shore birds
  #6 <- sea birds
  #7 <- land birds
#removes samples with species codes UNDU, OTHR, OHDU 

dabbling <- c("MALL", "WODU", "AMWI", "NOPI", "GADW","BWTE", "AGWT", "CITE", "NSHO",
              "ABDU", "MBDH", "MODU", "MEDU", "EUWI")
dabbling.df <- data.frame(species.code = dabbling, species.group = rep(1, each = length(dabbling)))
dabbling.matrix <- as.matrix(dabbling.df)

diving <- c("REDH", "HOME", "RNDU", "RUDU", "LESC", "GRSC", "COME", "CANV")
diving.df <- data.frame(species.code = diving, species.group = rep(2, each = length(diving)))
diving.matrix <- as.matrix(diving.df)

sea.duck <- c("BUFF", "HARD", "COGO", "BAGO", "COEI", "SUSC", "WWSC", "BLSC", 'LTDU')
sea.duck.df <- data.frame(species.code = sea.duck, species.group = rep(3, each = length(sea.duck)))
sea.duck.matrix <- as.matrix(sea.duck.df)

aserines <- c("HAGO", "CAGO", "LSGO", "ATBR", "CACG", "ROGO", "ACGO", "GSGO",
              "MUSW", "TUSW", "TRUS")
aserines.df <- data.frame(species.code = aserines, species.group = rep(4, each = length(aserines)))
aserines.matrix <- as.matrix(aserines.df)

shore <- c("AMOY", "RUTU", "SBDO", "REKN", "DOWI", "WESA", "PAGP", "BBPL", "DUNL", 
           "WATA", "SAND", "SEPL", "CLRA", "SESA", "LESA", "PUSA", "KILL", "SOSA",
           "SPSA", "HAST", "BTCU", "HAMO")
shore.df <- data.frame(species.code = shore, species.group = rep(5, each = length(shore)))
shore.matrix <- as.matrix(shore.df)

sea.bird <- c("COTE", "ROST", "ARTE", "RBGU", "GBBG", "GWGU", "HERG", "WTSH", "LAAL", 
              "DCCO", "HACO", "AMCO", "COLO", "AWPE")
sea.bird.df <- data.frame(species.code = sea.bird, species.group = rep(6, each = length(sea.bird)))
sea.bird.matrix <- as.matrix(sea.bird.df)

land <- c('BAEA', 'RTHA', 'RSHA', 'SNOW', 'WITU', "COHA")
land.df <- data.frame(species.code = land, species.group = rep(7, each = length(land)))
land.matrix <- as.matrix(land.df)

overall.matrix <- rbind(dabbling.matrix, diving.matrix, sea.duck.matrix, aserines.matrix, shore.matrix, 
                 sea.bird.matrix, land.matrix)
overall.df <- as.data.frame(overall.matrix)

data <- inner_join(data, overall.df, by = c("species.code.full"="species.code"))

# setwd("~/Github")
# setwd("~/HP/Data")
saveRDS(data, "AVHS_samplingevent_speciesgroup.rds")

######################################################################################################
#add a column to tell which watershed
######################################################################################################

# setwd("~/Honors Thesis/Project/hydrologic_units")
setwd("~/HP/hydrologic_units")
huc4 <- shapefile("huc4.shp")
projection(huc4) <- CRS("+proj=longlat +ellps=WGS84")

pt <- data[, c("lat", "long")]
coordinates(pt) <- ~ long + lat
proj4string(pt) <- CRS("+proj=longlat +ellps=WGS84")
sel <- !(huc4$STATES %in% c("AK", "AS", "AK,CN", "HI", "PR", "GU", "MP", "VI"))
huc4 <- huc4[sel, ]

data$huc4 <- extract(huc4, pt)$HUC4
data <- data[!(is.na(data$huc4)), ]

###################################################################################
#Add Sampling Events
###################################################################################
#create a data frame to match all possible sampling locations to all possible weeks and years
#weeks defined by the epidemiological calendar

# data <- readRDS("~/HP/Data/AVHS_samplingevent.rds")

# locations <- unique(cbind(data$long, data$lat))
# total.locations <- length(locations[,1])
# total.years <- length(unique(data$collection.year))
# total.weeks <- 53 * total.years
# location.x <- rep(locations[ ,1], total.weeks)
# location.y <- rep(locations[, 2], total.weeks)
# week <- rep(seq(1, 53, 1), each = total.locations, times = total.years)
# year <- rep(c("2007", "2008", "2009", "2010", "2015"), each = total.locations*53)
# assign.event <- data.frame(location.x = location.x, location.y=location.y, week = week,
#                            year=year)
# assign.event$prelim.number <- seq(1, length(location.x), 1)
# backup.event <- assign.event
# 
# #find location, month, year combinations with data and assign a sample event number to each
# data$collection.date <- mdy(data$collection.date.char)
# data$week <- epiweek(data$collection.date)
# data$collection.year <- as.factor(data$collection.year)
# test <- inner_join(data, assign.event, by = c("long"="location.x", "lat"="location.y", "week"="week", "collection.year"="year"))
# unique <- unique(test$prelim.number)
# temp <- filter(assign.event, prelim.number %in% unique)
# temp$event.number.week <- seq(1, length(temp$prelim.number))
# test <- inner_join(test, temp, by = c("long"="location.x", "lat"="location.y", "week"="week", "collection.year"="year",
#                                       'prelim.number'="prelim.number"))
# data <- test
# setwd("~/HP/Data")
# saveRDS(data, "AVHS_samplingevent.rds")

#assign sampling event numbers based on collection.group name, not lat and long
locations <- unique(data$collection.group)
total.locations <- length(locations)
total.years <- length(unique(data$collection.year))
total.weeks <- 53 * total.years
week <- rep(seq(1, 53, 1), each = total.locations, times = total.years)
year <- rep(c("2007", "2008", "2009", "2010", "2015"), each = total.locations*53)
assign.event <- data.frame(location = locations, week = week,
                           year=year)
assign.event$prelim.number <- seq(1, length(assign.event$location), 1)
backup.event <- assign.event

#find location, month, year combinations with data and assign a sample event number to each
data$collection.date <- mdy(data$collection.date.char)
data$week <- epiweek(data$collection.date)
data$collection.year <- as.factor(data$collection.year)
test <- inner_join(data, assign.event, by = c("collection.group"="location", "week"="week", "collection.year"="year"))
unique <- unique(test$prelim.number)
temp <- filter(assign.event, prelim.number %in% unique)
temp$event.number.week <- seq(1, length(temp$prelim.number))
test <- inner_join(test, temp, by = c("collection.group"="location", "week"="week", "collection.year"="year",
                                      'prelim.number'="prelim.number"))
data <- test
setwd("~/HP/Data")
saveRDS(data, "AVHS_samplingevent_locationname.rds")

##########################################################################################################################
#extract number of samples, number of positive samples for each sampling event and by species group/sampling event
##########################################################################################################################

#isolate unique sampling events including the month, year, and watershed of sampling
events <- dplyr::select(data, event.number.week, collection.month, collection.year, huc4)
events <- events[!duplicated(events[1:1]), ]

#create a dataframe that includes the number of (positive) samples per species group per sampling event
event <- data.frame(sample.event = rep(events$event.number.week, 7), month = rep(events$collection.month, 7),
                    year = rep(events$collection.year, 7), watershed = rep(events$huc4, 7), 
                    species.group = rep(seq(1, 7), each = length(events$event.number.week)), n = NA, y = NA)

for (i in 1:length(event$sample.event)) {
  hold <- filter(data, event.number.week == event[i, 1])
  holdy <- filter(hold, species.group == event[i, 5])
  holdz <- filter(holdy, AIpcr_susneg == "positive")
  event[i, 6] <- length(holdy$subjectID)
  event[i, 7] <- length(holdz$subjectID)
}

setwd("~/Github")
saveRDS(event, "locationsamplingevent_n_y_speciesgroup.rds")

