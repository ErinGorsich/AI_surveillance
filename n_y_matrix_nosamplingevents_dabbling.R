######################################################################
######################################################################
#Create an N and Y matrix from real data with no sampling events
#dabbling only
######################################################################
######################################################################

#read in sorted data
data <- readRDS("/home/webblab/Documents/HP/locationsamplingevent_n_y_speciesgroup.rds")

#remove sampling events and isolate species group 1 only
data <- data[data$species.group == 1, ]

#convert years to numbers
data$year.number <- NA
for (i in 1:10) {
  if (data$year[i] == "2007") {
    data$year.number[i] <- 1
  } else if (data$year[i] == "2008") {
    data$year.number [i] <- 2
  } else if (data$year[i] == "2009") {
    data$year.number[i] <- 3
  } else if (data$year[i] == "2010") {
    data$year.number[i] <- 4
  } else if (data$year[i] == "2015") {
    data$year.number[i] <- 5
  }
}

chart$month <- NA
chart$year <- NA
chart$watershed <- NA
chart$n <- NA #total number of samples
chart$y <- NA #total number of positive 

for (i in 1:length(unique(data$month))) {
  for (j in 1:length(unique(data$year))) {
    for (k in 1:length(unique(data$watershed))) {
      hold <- data.frame(sample.event = NA, month = i, year = j, watershed = k)
      hold$sample.event <- data$sample.event[data$month == i & data$year == j & data$watershed == k]
    }
  }
}

n <- array(NA, dim = c(12, 5, length(unique(data$watershed))))
for (i in 1:length(unique(data$month))) {
  for (j in 1:5) {
    for (k  in 1:length(unique(data$watershed))) {
      n[i, j, k] <- data$n[data$month == i & data$year == j & data$watershed == k]
    }
  }
}

y <- array(NA, dim = c(12, 5, length(unique(data$watershed))))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:length(unique(data$watershed))) {
      y[i, j, k] <- data$y[data$month == i & data$year == j & data$watershed == k]
    }
  }
}

saveRDS(n, "home/webblab/Documents/HP/n_matrix_nosamplingevents_dabbling.rds")
saveRDS(n, "home/webblab/Documents/HP/y_matrix_nosamplingevents_dabbling.rds")