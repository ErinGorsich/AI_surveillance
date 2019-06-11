############################################################################################
############################################################################################
#Simulated Data Attempt 6
#create a data set of apparent prevalence based on true prevalence given to the model
############################################################################################
############################################################################################

data <- readRDS("~/Github/locationsamplingevent_n_y_speciesgroup.rds")

data <- data[data$species.group %in% c(1,2,4), ]
for (i in 1:length(data$species.group)) {
  if (data$species.group[i] == 4) {
    data$species.group[i] <- 3
  } else {data$species.group[i] <- data$species.group[i]}
}
data$YEAR <- NA
for (i in 1:length(data$year)) {
  if (data$year[i] == "2007") {
    data$YEAR[i] <- 1
  } else if (data$year[i] == "2008") {
    data$YEAR[i] <- 2
  } else if (data$year[i] == "2009") {
    data$YEAR[i] <- 3
  } else if (data$year[i] == "2010") {
    data$YEAR[i] <- 4
  } else {
    data$YEAR[i] <- 5
  }
}
data$month <- as.numeric(data$month)

huc <- data.frame(watershed = unique(data$watershed), huc = seq(1,length(unique(data$watershed))))

for (i in 1:length(data$watershed)) {
  data$HUC[i] <- huc$huc[huc$watershed == data$watershed[i]]
}

nmonths <- length(unique(data$month))
nyears <- length(unique(data$YEAR))
nhucs <- length(unique(data$HUC))

tally <- data.frame(species.group = rep(seq(1,3), each = nyears*nmonths*nhucs),
                    month = rep(seq(1,12), times = nyears*nhucs*3),
                    year = rep(seq(1,5), each = 12, times = nhucs*3),
                    huc = rep(seq(1, nhucs), each = nyears*nmonths, rep = 3), 
                    n = NA, y = NA, app.prev = NA)
for (h in 1:3) {
  for (i in 1:nmonths) {
    for (j in 1:nyears) {
      for (k in 1:nhucs) {
        hold <- data$n[data$species.group == h & data$month == i & data$YEAR == j & data$HUC == k]
        tally$n[tally$species.group == h & tally$month == i & tally$year == j & tally$huc == k] <-
          sum(hold)
        hold <- data$y[data$species.group == h & data$month == i & data$YEAR == j & data$HUC ==k]
        tally$y[tally$species.group == h & tally$month == i & tally$year == j & tally$huc ==k] <- 
          sum(hold)
      }
    }
  }
}
for (i in 1:length(tally$species.group)) {
  tally$app.prev[i] <- tally$y[i]/tally$n[i]
  if (is.na(tally$app.prev[i] == TRUE)) {
    tally$app.prev[i] <- 0
  }
}

tally$true.prev <- NA
spec = 0.999
sens = 0.83
for (i in 1:length(tally$species.group)){
  tally$true.prev[i] <- (tally$app.prev[i]+(spec-1))/(spec+(sens-1))
  if (tally$true.prev[i] > 1) {
    tally$true.prev[i] <- 1
  } else if (tally$true.prev[i] < 0) {
    tally$true.prev[i] <- abs(tally$true.prev[i])
  } else {
    tally$true.prev[i] <- tally$true.prev[i]
  }
}
for (i in 1:length(tally$true.prev)) {
  if (tally$app.prev[i]<1 & tally$true.prev[i]==1) {
    tally$true.prev[i] <- tally$true.prev[i] - (1-tally$app.prev[i])
  } else {
    tally$true.prev[i] <- tally$true.prev[i]
  }
}

############################################################################################################
#simulate data using model function
#species.group==1
############################################################################################################
nmonths <- length(unique(tally$month))
nyears <- length(unique(tally$year))
nhucs <- length(unique(tally$huc))

n <- array(dim=c(nmonths, nyears, nhucs), NA)
lambda <- array(dim=c(nmonths, nyears, nhucs), NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      n[i,j,k] <- tally$n[tally$species.group==1 & tally$month==i & tally$year==j & tally$huc==k]
      lambda[i,j,k] <- tally$true.prev[tally$species.group==1 & tally$month==i & tally$year==j & tally$huc==k]
    }
  }
}
Se <- rbeta(nmonths*nyears*nhucs, 20.833, 4.148)
Sp <- rbeta(nmonths*nyears*nhucs, 8.403, 1.001)
Se <- mean(Se)
Sp <- mean(Sp)
#Sp=0.8952586, Se=0.8329987

simulated <- base.model(n=n, lambda=lambda, Se=Se, Sp=Sp)

saveRDS(simulated, "/home/webblab/Documents/HP/sim_data_attempt6.rds")
saveRDS(n, "/home/webblab/Documents/HP/sim_data_6_n.rds")

#######################################################################################################
#create a y matrix for model fittiing based on the generated apparent prevalences
#######################################################################################################
y <- array(dim=c(nmonths, nyears, nhucs), NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      positive <- simulated[i,j,k]*n[i,j,k]
      y[i,j,k] <- round(positive, digits = 0)
    }
  }
}

saveRDS(y, "/home/webblab/Documents/HP/sim_data_6_y.rds")