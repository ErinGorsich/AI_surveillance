library(simstudy)
library(dplyr)

see <- readRDS("/home/webblab/Documents/HP/AVHS_samplingevent_speciesgroup.rds")
see$huc <- NA
see$year <- NA
see$month <- as.numeric(see$collection.month)
for (i in 1:length(see$collection.year)) {
  if (see$collection.year[i] == "2007") {
    see$year[i] <- 1
  } else if (see$collection.year[i] == "2008") {
    see$year[i] <- 2
  } else if (see$collection.year[i] == "2009") {
    see$year[i] <- 3
  } else if (see$collection.year[i] == "2010") {
    see$year[i] <- 4
  } else {
    see$year[i] <- 5
  }
}

huc <- data.frame(watershed = unique(see$watershed), huc = NA)
huc$huc = seq(1:length(huc$watershed))

for (i in 1:length(see$watershed)) {
  y <- see$watershed[i]
  see$huc[i] <- huc$huc[huc$watershed == y]
}

total <- count(see, species.group, month, year, huc)
total$total.samples <- total$n
total$n <- NULL
results <- count(see, species.group, month, year, huc, AIpcr_susneg)
results <- results[results$AIpcr_susneg == "positive", ]
results$y <- results$n
results$n <- NULL

long <- data.frame(species.group = as.factor(rep(seq(1,3), each = 202*3*5)),
                   huc = rep(seq(1,202), times = 60*3),
                   month = rep(seq(1,12), each = 202, times = 5*3),
                   year = rep(seq(1,5), each = 202*12, times = 3),
                   app.prev = NA)

long <- left_join(long, total)
long <- left_join(long, results)
long$AIpcr_susneg <- NULL

for (i in 1:length(long$total.samples)) {
  if (is.na(long$total.samples[i]) == TRUE) {
    long$total.samples[i] <- 0
  } else {
    long$total.samples[i] <- long$total.samples[i]
  }
}

for (i in 1:length(long$y)) {
  if (is.na(long$y[i]) == TRUE) {
    long$y[i]<- 0
  } else {
    long$y[i] <- long$y[i]
  }
} 

for (i in 1:length(long$total.samples)) {
  long$app.prev[i] <- long$y[i]/long$total.samples[i]
}

##################################################################################
#simulate a positive number of samples for each species group, huc, month, year combination
#use a binomial distribution (n = # of samples, p = app.prev) to simulate random 
#number of positive samples (sim.pos), then recalculate app.prev (sim.app.prev)
###################################################################################

long$sim.pos <- NA
long$sim.app.prev <- NA
for (i in 1:length(long$species.group)) {
  hold <- rbinom(long$total.samples[i], 1, p = long$app.prev[i])
  long$sim.pos[i] <- length(hold[hold == 1])
  long$sim.app.prev[i] <- long$sim.pos[i]/long$total.samples[i]
}

saveRDS(long, "/home/webblab/Documents/HP/sim_data_2.rds")

#########################################################################################
#make simulated n and y matrices for model fitting
#species.group == 1
#########################################################################################

n <- array(NA, dim = c(12, 5, 202))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:202) {
      n[i, j, k] <- long$total.samples[long$species.group == 1 & long$month == i & long$year == j & long$huc == k]
    }
  }
}

y <- array(NA, dim = c(12, 5, 202))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:202) {
      y[i, j, k] <- long$sim.pos[long$species.group == 1 & long$month == i & long$year == j & long$huc == k]
    }
  }
}

saveRDS(n, "/home/webblab/Documents/HP/sim_data_2_n.rds")
saveRDS(y, "/home/webblab/Documents/HP/sim_data_2_y.rds")

#########################################################################################
#make actual n and y matrices for model fitting (test mixing issues)
#species.group == 1
#########################################################################################

n <- array(NA, dim = c(12, 5, 222))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:222) {
      n[i, j, k] <- long$total.samples[long$species.group == 1 & long$month == i & long$year == j & long$huc == k]
    }
  }
}

y <- array(NA, dim = c(12, 5, 222))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:222) {
      y[i, j, k] <- long$y[long$species.group == 1 & long$month == i & long$year == j & long$huc == k]
    }
  }
}

saveRDS(n, "/home/webblab/Documents/HP/real_data_n_mall.rds")
saveRDS(y, "/home/webblab/Documents/HP/sim_data_y_mall.rds")