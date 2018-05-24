#cleans data and assigns sampling events

library(lubridate)
library(stringr)

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
data$collection.month[ai2$collection.month == "Aug"] <- "08"
data$collection.month[ai2$collection.month == "Dec"] <- "12"
data$collection.month[ai2$collection.month == "Feb"] <- "02"
data$collection.month[ai2$collection.month == "Jan"] <- "01"
data$collection.month[ai2$collection.month == "Jul"] <- "07"
data$collection.month[ai2$collection.month == "Jun"] <- "06"
data$collection.month[ai2$collection.month == "May"] <- "05"
data$collection.month[ai2$collection.month == "Nov"] <- "11"
data$collection.month[ai2$collection.month == "Oct"] <- "10"
data$collection.month[ai2$collection.month == "Sep"] <- "09"
