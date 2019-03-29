######################################################################
######################################################################
#Create an N and Y matrix from real data with no sampling events
#dabbling only
######################################################################
######################################################################

#read in sorted data
data <- readRDS("~/Github/locationsamplingevent_n_y_speciesgroup.rds")

#remove sampling events and isolate species group 1 only
data <- data[, 2:7]
data <- data[data$species.group == 1, ]

#convert years to numbers
data$year.number <- NA
for (i in 1:length(data$year)) {
  if data$year[i] == "2007" {
    data$year.number[i] <- 1
  } else if data$year[i] == "2008" {
    data$year.number [i] <- 2
  } else if data$year[i] == "2009" {
    data$year.number[i] <- 3
  } else if data$year[i] == "2010" {
    data$year.number[i] <- 4
  } else if data$year[i] == "2015" {
    data$year.number[i] <- 5
  }
}

n <- array(NA, dim = c(12, 5, 222))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:222) {
      n[i, j, k] <- chart$total[chart$species.group == 1 & chart$month == i & chart$year == j & chart$huc == k]
    }
  }
}

y <- array(NA, dim = c(12, 5, 222))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:222) {
      y[i, j, k] <- chart$positive[chart$species.group == 1 & chart$month == i & chart$year == j & chart$huc == k]
    }
  }
}