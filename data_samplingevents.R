#cleans data and assigns sampling events

library(lubridate)
library(stringr)

##################################################################################
#data read in
#################################################################################
setwd("~/Honors Thesis/Project")
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

setwd("~/Github")
saveRDS(data, "AVHS_samplingevent.rds")

###################################################################################
#Add Sampling Events
###################################################################################
#create a data frame to match all possible sampling locations to all possible weeks, months, years
  #month = 28 day period (13/year)
  #weeks belong to the month of the first day of the week
# data <- readRDS("~/HP/Data/AVHS_samplingevent.rds")

data$week <- epiweek(data$collection.date)
data$day <- yday(data$collection.date)

#assign.event
locations <- unique(cbind(data$long, data$lat))
total.locations <- length(locations[,1])
total.years <- length(unique(data$collection.year))
total.days <- 366 * total.years

location.x <- rep(locations[ ,1], total.days)
location.y <- rep(locations[, 2], total.days)
day <- rep(seq(1, 366, 1), each = total.locations, times = total.years)
year <- rep(c("2007", "2008", "2009", "2010", "2015"), each = total.locations*366)

assign.event <- data.frame(location.x = location.x, location.y=location.y, day = day,
                           year=year)
assign.event$prelim.number <- seq(1, length(location.x), 1)
backup.event <- assign.event

test <- inner_join(data, assign.event, by = c("long"="location.x", "lat"="location.y", "day"="day", "collection.year"="year"))

#start here for testing
unique <- unique(data$prelim.number)
temp <- assign.event  %>%
  filter(prelim.number %in% unique)

#find location, month, year combinations with data and assign a sample event number to each
  #add a week number column to the sampling data
# data$collection.date.transformed <- mdy(data$collection.date.char)
# data$week <- week(data$collection.date.transformed)
# data$year <- NA
# for (i in 1:length(data$collection.year)) {
#   if (data[i, ]$collection.year == "2007") {
#     data[i, ]$year <- 0
#   } else if (data[i, ]$collection.year == "2008") {
#     data[i, ]$year <- 1
#   } else if (data[i, ]$collection.year == "2009") {
#     data[i, ]$year <- 2
#   } else if (data[i, ]$collection.year == "2010") {
#     data[i, ]$year <- 3
#   } else if (data[i, ]$collection.year == "2015") {
#     data[i, ]$year <- 4
#   }
# }
# data$week.intermediate <- NA
# for (i in 1:10) {
#   if (data[i, ]$week == "53") {
#     data[i, ]$week <- "52"
#   } else {
#     data[i, ]$week <- data[i, ]$week
#   }
# } 
# data$week.final <- NA
# for (i in 1:length(data$week)) {
#   data[i, ]$week.final <- (data[i,]$year*52) + data[i, ]$week
# }

#assign a week (1-52) for each sample taken
data$week <- epiweek(data$collection.date)
temp <- data[data$week == "53", ]
temp <- na.omit(temp)
# for (i in 1:length(data$week)) {
#   if (data[i, ]$week == "53" & data[i, ]$collection.month == "12") {
#     data[i, ]$week <- 52
#   } else if (data[i, ]$week == "53" & data[i, ]$collection.month == "01") {
#     data[i, ]$week <- 1
#   } else {
#     data[i, ]$week <- data[i, ]$week
#   }
# }
data$week <- epiweek(data$collection.date)
temp <- data[data$week == "53", ]
temp <- na.omit(temp)
for (i in 1:length(temp$week)) {
  if (temp[i, ]$week == "53" & temp[i, ]$collection.month == "12") {
    temp[i, ]$week <- 52
  } else if (temp[i, ]$week == "53" & temp[i, ]$collection.month == "01") {
    temp[i, ]$week <- 1
  } else {
    temp[i, ]$week <- temp[i, ]$week
  }
}

data$collection.week <- temp$week[match(data$subjectID, temp$subjectID)]
for (i in 1:length(data$collection.week)) {
  if (data[i, ]$collection.week == "NA") {
    data[i, ]$collection.week <- data[i, ]$week
  } else {
    data[i, ]$collection.week <- data[i, ]$collection.week
  }
}
#add a sample event column to the main dataset that assigns a sample event number to each sample



