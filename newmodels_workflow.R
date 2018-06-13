####################################################################################
#workflow for updated models
###################################################################################

#load packages
library(rjags)
library(coda)

#load models
setwd("~/Github/AI_surveillance")
source("~/Github/AI_surveillance/define_models.r")
source("~/Github/AI_surveillance/define_neighborhood.r")

#read in sampling data
data <- readRDS('~/Github/samplingevent_n_y_speciesgroup.rds')
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
nsamplingevents <- length(unique(data$sample.event)) #13589
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
#######################################################################################
#base model
#########################################################################################

jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, nmonths=nmonths,
                  nyears=nyears, nhucs = nhucs, n = n, y=y, month = data$month[1:nsamplingevents],
                  year = data$yearid[1:nsamplingevents], huc = data$hucid[1:nsamplingevents])
#need initial values for lambda
# jags.inits <- function(){
#   list("Se" = runif(1, 0.6, 1), 'Sp' = runif(1, 0.6, 1), 
#        "pi" = array(runif(nsamplingevents*nmonths*nyears), dim = c(nmonths, nyears, nsites)))
# }
variable.names = c('Se', 'Sp', 'lambda', 'pi') #apparent prevalence?

nadapt <- 1000
niter <- 10
thin <- 1

setwd("~/Github/AI_surveillance")
base.mod <- jags.model(file = "base_sampling_events.txt", data=jags.data,
                      n.chains=3, n.adapt=nadapt)
saveRDS(base.mod, "modelruns/base_sampling_events_adapt.rds")
base.mod.fit <- coda.samples(model=base.mod, variable.names=variable.names, n.iter=niter, 
                             thin=thin)
saveRDS(base.mod.fit, "modelruns/base_sampling_event_fit.rds")

###################################################################################################
#spatial model
###################################################################################################
W <- define.neighborhood(method="queens")
B <- scaleW(W)
D <- diag(rowSums(W))
I <- diag(1, nrow = dim(W)[1], ncol = dim(W)[1])

jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, nmonths=nmonths,
                  nyears=nyears, nhucs = nhucs, n = n, y=y, month = data$month[1:nsamplingevents],
                  year = data$yearid[1:nsamplingevents], huc = data$hucid[1:nsamplingevents], W=W, I=I)
#need initial values for lambda
# jags.inits <- function(){
#   list("Se" = runif(1, 0.6, 1), 'Sp' = runif(1, 0.6, 1), 
#        "pi" = array(runif(nsamplingevents*nmonths*nyears), dim = c(nmonths, nyears, nsites)))
# }
variable.names = c('Se', 'Sp', 'lambda', 'pi') #apparent prevalence?

nadapt <- 1000
niter <- 10
thin <- 1

setwd("~/Github/AI_surveillance")
spatial.mod <- jags.model(file = "icar_sampling_events.txt", data=jags.data, n.chains = 3,
                          n.adapt=nadapt)
saveRDS(spatial.mod, "modelruns/icar_sampling_events_adapt.txt")
spatial.mod.fit <- coda.samples(model=spatial.mod, variable.names=variable.names, n.iter=niter,
                                thin=thin)
saveRDS(spatial.mod.fit, "modelruns/icar_sample_events_fit.txt")
