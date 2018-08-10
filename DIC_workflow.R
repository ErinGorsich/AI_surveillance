#######################################################################################
#######################################################################################
#run DIC calculations for all models
#######################################################################################
#######################################################################################

#load packages
library(rjags)
library(coda)
library(dclone)

#source files
setwd("~/Github/AI_surveillance")
source("define_models.r")
source("define_neighborhood.r")

#read in data
setwd("~/Github")
data <- readRDS("locationsamplingevent_n_y_speciesgroup.rds")

#read in sampling data
data$month <- as.numeric(data$month)
data$year <- as.numeric(as.character(data$year))

#limit data to only the species groups wanted and change group 4 to group 3
data <- data[data$species.group %in% c("1", "2", "4"), ]
data$species.group[data$species.group == "4"] <- 3

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
nsamplingevents <- length(unique(data$sample.event)) #12349
nspecies <- length(unique(data$species.group)) #3
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

#set number of iterations
nadapt <- 10000
nburn <- 50000
niter <- 20000
thin <- 10

#############################################################################################
#base model 
#############################################################################################
jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, 
                  nmonths = nmonths, nyears = nyears, nhucs = nhucs, n = n, y = y, 
                  month = data$month[1:nsamplingevents],
                  year = data$yearid[1:nsamplingevents], huc = data$hucid[1:nsamplingevents])
variable.names = c('Se', 'Sp', 'pi')

base.adapt <- jags.model(file = "~/Github/AI_surveillance/base_sampling_events.txt", data = jags.data,
                         n.chains = 3, n.adapt = nadapt)
update(base.adapt, nburn)
base.dic <- dic.samples(model = base.adapt, n.iter = niter, thin=thin, type = "pD")
saveRDS(base.dic, "~/Github/base_dic.rds")

rm(base.adapt, base.dic)

###########################################################################################
#queens
###########################################################################################

W <- define.neighborhood(method="queens")
B <- scaleW(W)
D <- diag(rowSums(W))
I <- diag(1, nrow = dim(W)[1], ncol = dim(W)[1])

jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, nmonths=nmonths,
                  nyears=nyears, nhucs = nhucs, n = n, y=y, month = data$month[1:nsamplingevents],
                  year = data$yearid[1:nsamplingevents], huc = data$hucid[1:nsamplingevents], W=W, I=I)
variable.names = c('Se', 'Sp', 'pi')

queens.adapt <- jags.model(file = "~/Github/AI_surveillance/icar_sampling_events.txt", data = jags.data,
                         n.chains = 3, n.adapt = nadapt)
update(queens.adapt, nburn)
queens.dic <- dic.samples(model = queens.adapt, n.iter = niter, thin=thin, type = "pD")
saveRDS(queens.dic, "~/Github/queens_dic.rds")

rm(queens.adapt, queens.dic)

###############################################################################################
#weighted queens
###############################################################################################

W <- define.neighborhood(method="weightedqueens")
B <- scaleW(W)
D <- diag(rowSums(W))
I <- diag(1, nrow = dim(W)[1], ncol = dim(W)[1])

jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, nmonths=nmonths,
                  nyears=nyears, nhucs = nhucs, n = n, y=y, month = data$month[1:nsamplingevents],
                  year = data$yearid[1:nsamplingevents], huc = data$hucid[1:nsamplingevents], W=W, I=I)
variable.names = c('Se', 'Sp', 'pi')

wtqueens.adapt <- jags.model(file = "~/Github/AI_surveillance/icar_sampling_events.txt", data = jags.data,
                           n.chains = 3, n.adapt = nadapt)
update(wtqueens.adapt, nburn)
wtqueens.dic <- dic.samples(model = wtqueens.adapt, n.iter = niter, thin=thin, type = "pD")
saveRDS(wtqueens.dic, "~/Github/wtqueens_dic.rds")

rm(wtqueens.adapt, wtqueens.dic)

###############################################################################################
#network
###############################################################################################

W <- define.neighborhood(method="network")
B <- scaleW(W)
D <- diag(rowSums(W))
I <- diag(1, nrow = dim(W)[1], ncol = dim(W)[1])

jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, nmonths=nmonths,
                  nyears=nyears, nhucs = nhucs, n = n, y=y, month = data$month[1:nsamplingevents],
                  year = data$yearid[1:nsamplingevents], huc = data$hucid[1:nsamplingevents], W=W, I=I)
variable.names = c('Se', 'Sp', 'pi')

network.adapt <- jags.model(file = "~/Github/AI_surveillance/icar_sampling_events.txt", data = jags.data,
                           n.chains = 3, n.adapt = nadapt)
update(network.adapt, nburn)
network.dic <- dic.samples(model = network.adapt, n.iter = niter, thin=thin, type = "pD")
saveRDS(network.dic, "~/Github/network_dic.rds")

rm(network.adapt, network.dic)

####################################################################################################
#weighted network
###################################################################################################

W <- define.neighborhood(method="weightednetwork")
B <- scaleW(W)
D <- diag(rowSums(W))
I <- diag(1, nrow = dim(W)[1], ncol = dim(W)[1])

jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, nmonths=nmonths,
                  nyears=nyears, nhucs = nhucs, n = n, y=y, month = data$month[1:nsamplingevents],
                  year = data$yearid[1:nsamplingevents], huc = data$hucid[1:nsamplingevents], W=W, I=I)
variable.names = c('Se', 'Sp', 'pi')

wtnetwork.adapt <- jags.model(file = "~/Github/AI_surveillance/icar_sampling_events.txt", data = jags.data,
                            n.chains = 3, n.adapt = nadapt)
update(wtnetwork.adapt, nburn)
wtnetwork.dic <- dic.samples(model = wtnetwork.adapt, n.iter = niter, thin=thin, type = "pD")
saveRDS(wtnetwork.dic, "~/Github/wtnetwork_dic.rds")

rm(wtnetwork.adapt, wtnetwork.dic)

#######################################################################################################
#temporal
#######################################################################################################

jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, nmonths=nmonths,
                  nyears=nyears, nhucs = nhucs, n = n, y=y, month = data$month[1:nsamplingevents],
                  year = data$yearid[1:nsamplingevents], huc = data$hucid[1:nsamplingevents])
variable.names = c('Se', 'Sp', 'pi')

temporal.adapt <- jags.model(file = "~/Github/AI_surveillance/ar1_sampling_events.txt", data = jags.data,
                             n.chains = 3, n.adapt = nadapt)
update(temporal.adapt, nburn)
temporal.dic <- dic.samples(model = temporal.adapt, n.iter=niter, thin=thin, type="pD")
saveRDS(temporal.dic, "~/Github/temporal_dic.rds")

