####################################################################################
####################################################################################
#workflow to run diagnostic tests on the base model to improve speed
####################################################################################
####################################################################################

#load packages
library(rjags)
library(coda)

setwd("~/Github/AI_surveillance")
source("define_models.r")

#read in data
setwd("~/Github")
data <- readRDS('samplingevent_n_y_speciesgroup.rds')

#change month and year columns to numerics
data$month <- as.numeric(data$month)
data$year <- as.numeric(as.character(data$year))

#add a hucid column
temp <- data.frame(watershed = (unique(data$watershed)))
temp$hucid <- seq(1, length(temp[,1]))
data$hucid <- temp$hucid[match(data$watershed, temp$watershed)]

#add a yearid column
temp <- data.frame(year = unique(data$year))
temp$yearid <- seq(1, length(temp[,1]))
data$yearid <- temp$yearid[match(data$year, temp$year)]

#change month column to integers
data$month <- as.integer(data$month)

#define variables
nsamplingevents <- length(unique(data$sample.event)) #12787
nspecies <- length(unique(data$species.group)) #7
nmonths <- length(seq(1, 12)) #12
nyears <- length(unique(data$year)) #5
nhucs <- length(unique(data$watershed)) #195

y <- matrix(NA, nrow = nsamplingevents, ncol = nspecies)
for (s in 1:nsamplingevents) {
  for (l in 1:nspecies){
    if(length( data$y[data$sample.event == s & data$species.group == l]) >0 ){
      y[s, l] <- data$y[data$sample.event == s & data$species.group == l]}
  }
}

n <- matrix(NA, nrow = nsamplingevents, ncol = nspecies)
for (s in 1:nsamplingevents) {
  for (l in 1:nspecies){
    if(length(data$n[data$sample.event == s & data$species.group == l]) >0 ){
      n[s, l] <- data$n[data$sample.event == s & data$species.group == l]}
  }
}

nadapt <- 1000
niter <- 100
thin <- 2

##########################################################################################
#test one
#run model with sampling events but without species groups
##########################################################################################
ai <- readRDS("~/HP/Data/AVHS_samplingevent.rds")
sampling.events <- data.frame(sampling.event = unique(ai$event.number.week),
                              watershed = NA, n = NA, y = NA)

for (i in 1:length(sampling.events$sampling.event)) {
  hold <- ai[ai$event.number.week == i, ]
  sampling.events[sampling.events$sampling.event == i, ]$watershed <- hold$huc4[1]
  hold$tally <- 1
  sampling.events[sampling.events$sampling.event == i, ]$n <- sum(hold$tally)
  hold.y <- hold[hold$AIpcr_susneg == "positive", ]
  sampling.events[sampling.events$sampling.event == i, ]$y <- sum(hold.y$tally)
}

jags.data <- list(nsamplingevents = nsamplingevents, nmonths = nmonths, nyears = nyears, 
                  nhucs = nhucs, n = sampling.events$n, y = sampling.events$y, 
                  month = data$month[1:nsamplingevents], year = data$yearid[1:nsamplingevents], 
                  huc = data$hucid[1:nsamplingevents])
variable.names = c("Se", "Sp", "pi")

setwd("~/Github/AI_surveillance")
base.nospecies.mod <- jags.model(file = "base_sampling_test_nospecies.txt", data = jags.data,
                                 n.chains = 3, n.adapt=nadapt)
setwd("~/Github/AI_surveillance/model runs")
saveRDS(base.nospecies.mod, "base_test_nospecies_adapt.rds")
update(base.nospecies.mod, nadapt)
base.nospecies.mod.fit <- coda.samples(model=base.nospecies.mod, variable.names = variable.names,
                                       n.iter = niter, thin=thin)
saveRDS(base.nospecies.mod.fit, "base_test_nospecies_fit.rds")

#########################################################################################
#test 2
#run the base model seperately for each species
#########################################################################################

#run only dabbling ducks
data.dabbling <- data[data$species.group == "1", ]

y.dabbling <- as.data.frame(y)
y.dabbling <- y.dabbling$V1

n.dabbling <- as.data.frame(n)
n.dabbling <- n.dabbling$V1

jags.data <- list(nsamplingevents = nsamplingevents, nmonths = nmonths, nyears = nyears, 
                  nhucs = nhucs, n = n.dabbling, y = y.dabbling, month = data.dabbling$month,
                  year = data.dabbling$yearid, huc = data.dabbling$hucid)
variable.names = c("Se", "Sp", "pi")

setwd("~/Github/AI_surveillance")
base.dabbling.mod <- jags.model(file = "base_sampling_test_nospecies.txt", data = jags.data,
                                n.chains = 3, n.adapt=nadapt)
setwd("~/Github/AI_surveillance/modelruns")
saveRDS(base.dabbling.mod, "base_dabbling_only_adapt.rds")
base.dabbling.mod <- update(base.dabbling.mod, nadapt)
base.dabbling.mod.fit <- coda.samples(model = base.dabbling.mod, variable.names = variable.names,
                                      n.iter = niter, thin=thin)
saveRDS(base.dabbling.mod.fit, "base_dabbling_only_fit.rds")

#run only diving ducks
data.diving <- data[data$species.group == "2", ]

y.diving <- as.data.frame(y)
y.diving <- y.diving$v1

n.diving <- as.data.frame(n)
n.diving <- n.diving$V1

jags.data <- list(nsamplingevents = nsamplingevents, nmonths = nmonths, nyears = nyears, 
                  nhucs = nhucs, n = n.diving, y = y.diving, month = data.diving$month,
                  year = data.diving$yearid, huc = data.diving$hucid)
variable.names = c("Se", "Sp", "pi")

setwd("~/Github/AI_surveillance")
base.diving.mod <- jags.model(file = "base_sampling_test_nospecies.txt", data = jags.data,
                                n.chains = 3, n.adapt=nadapt)
setwd("~/Github/AI_surveillance/modelruns")
saveRDS(base.diving.mod, "base_diving_only_adapt.rds")
base.diving.mod <- update(base.diving.mod, nadapt)
base.diving.mod.fit <- coda.samples(model = base.diving.mod, variable.names = variable.names,
                                      n.iter = niter, thin=thin)
saveRDS(base.diving.mod.fit, "base_diving_only_fit.rds")

################################################################################################
#test 3
#run base model as is, but with only dabbling and diving ducks
################################################################################################

nspecies = 2

jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, 
                  nmonths = nmonths, nyears = nyears, nhucs = nhucs, n = n, y = y, 
                  month = data$month[1:nsamplingevents],
                  year = data$yearid[1:nsamplingevents], huc = data$hucid[1:nsamplingevents])

variable.names = c("Se", "Sp", "pi")

setwd("~/Github/AI_surveillance")
base.dandd.mod <- jags.model(file = "base_sampling_events.txt", data = jags.data, 
                             n.chains = 3, n.adapt=nadapt)
setwd("~/Github/AI_surveillance/modelruns")
saveRDS(base.dandd.mod, "base_dandd_adapt.rds")
base.dandd.mod <- update(base.dabbling.mod, nadapt)
base.dandd.mod.fit <- coda.samples(model = base.dandd.mod, variable.names = variable.names, 
                                   n.iter=niter, thin=thin)
saveRDS(base.dandd.mod.fit, "base_dand_fit.rds")