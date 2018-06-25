####################################################################################
#workflow for updated models
###################################################################################

#load packages
library(rjags)
library(coda)

Erin <- FALSE

#load models
if (Erin){
    setwd("~/Github/AI_surveillance")
    data <- readRDS('samplingevent_n_y_speciesgroup.rds')
} else {
    setwd("~/Github/AI_surveillance")
    data <- readRDS('~/Github/samplingevent_n_y_speciesgroup.rds')
}
source("define_models.r")
source("define_neighborhood.r")

#read in sampling data
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
nadapt <- 1000
niter <- 10
thin <- 1

# test run no species
# jags.data <- list(nsamplingevents = nsamplingevents,
#     nmonths = nmonths, nyears = nyears, nhucs = nhucs, n = n[ ,1], y = y[ ,1], 
#     month = data$month[1:nsamplingevents],
#     year = data$yearid[1:nsamplingevents], 
#     huc = data$hucid[1:nsamplingevents])
# base.mod <- jags.model(file = "base_sampling_test.txt", data=jags.data,
#     n.chains=3, n.adapt=nadapt)

# model with species and sampling events
jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, 
    nmonths = nmonths, nyears = nyears, nhucs = nhucs, n = n, y = y, 
    month = data$month[1:nsamplingevents],
    year = data$yearid[1:nsamplingevents], huc = data$hucid[1:nsamplingevents])
#need initial values for lambda
# jags.inits <- function(){
#   list("Se" = runif(1, 0.6, 1), 'Sp' = runif(1, 0.6, 1), 
#        "pi" = array(runif(nsamplingevents*nmonths*nyears), dim = c(nmonths, nyears, nsites)))
# }
variable.names = c('Se', 'Sp', 'lambda', 'pi') #apparent prevalence?



#setwd("~/Github/AI_surveillance")
base.mod <- jags.model(file = "base_sampling_events.txt", data=jags.data,
    n.chains=3, n.adapt=nadapt)
saveRDS(base.mod, "modelruns/base_sampling_events_adapt.rds")
update(base.mod, nadapt)
base.mod.fit <- coda.samples(model=base.mod, variable.names=variable.names, n.iter=niter, 
                             thin=thin)
saveRDS(base.mod.fit, "modelruns/base_sampling_event_fit.rds")

#plots
setwd("~/HP/Plots")
pdf("base_trace_density_sensitivity_specificity.pdf")
plot(base.mod.fit[,1:2])
dev.off()

pdf("base_trace_density_sensitivity_specificity.pdf")
plot(base.mod.fit[,1:2])
dev.off()

pdf("base_sensitivity_autocorrelation.pdf")
autocorr.plot(base.mod.fit[,'Se'], main="Sensitivity Autocorrelation")
dev.off()

pdf("base_specificity_autocorrelation.pdf")
autocorr.plot(base.mod.fit[, 'Sp'], main="Specificity Autocorrelation")
dev.off()

rm(base.mod, base.mod.fit)

###################################################################################################
#spatial model - queens
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
queens.mod <- jags.model(file = "icar_sampling_events.txt", data=jags.data, n.chains = 3,
                          n.adapt=nadapt)
saveRDS(queens.mod, "modelruns/icar_sampling_events_queens_adapt.rds")
update(queens.mod, nadapt)
queens.mod.fit <- coda.samples(model=queens.mod, variable.names=variable.names, n.iter=niter,
                                thin=thin)
saveRDS(queens.mod.fit, "modelruns/icar_sample_events_queens_fit.rds")

#plots
setwd("~/HP/Plots")
pdf("iCAR_queens_trace_density_sensitivity_specificity.pdf")
plot(queens.mod.fit[,1:2])
dev.off()

pdf("iCAR_queens_trace_density_sensitivity_specificity.pdf")
plot(queens.mod.fit[,1:2])
dev.off()

pdf("iCAR_queens_sensitivity_autocorrelation.pdf")
autocorr.plot(queens.mod.fit[,'Se'], main="Sensitivity Autocorrelation")
dev.off()

pdf("iCAR_queens_specificity_autocorrelation.pdf")
autocorr.plot(queens.mod.fit[, 'Sp'], main="Specificity Autocorrelation")
dev.off()

###################################################################################################
#spatial model - weighted queens
###################################################################################################

W <- define.neighborhood(method="weightedqueens")
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
wtqueens.mod <- jags.model(file = "icar_sampling_events.txt", data=jags.data, n.chains = 3,
                         n.adapt=nadapt)
saveRDS(wtqueens.mod, "modelruns/icar_sampling_events_wtqueens_adapt.rds")
update(wtqueens.mod, nadapt)
wtqueens.mod.fit <- coda.samples(model=wtqueens.mod, variable.names=variable.names, n.iter=niter,
                               thin=thin)
saveRDS(wtqueens.mod.fit, "modelruns/icar_sample_events_wtqueens_fit.rds")

#plots
setwd("~/HP/Plots")
pdf("iCAR_wtqueens_trace_density_sensitivity_specificity.pdf")
plot(wtqueens.mod.fit[,1:2])
dev.off()

pdf("iCAR_wtqueens_trace_density_sensitivity_specificity.pdf")
plot(wtqueens.mod.fit[,1:2])
dev.off()

pdf("iCAR_wtqueens_sensitivity_autocorrelation.pdf")
autocorr.plot(wtqueens.mod.fit[,'Se'], main="Sensitivity Autocorrelation")
dev.off()

pdf("iCAR_wtqueens_specificity_autocorrelation.pdf")
autocorr.plot(wtqueens.mod.fit[, 'Sp'], main="Specificity Autocorrelation")
dev.off()

###################################################################################################
#spatial model - network
###################################################################################################

W <- define.neighborhood(method="network")
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
network.mod <- jags.model(file = "icar_sampling_events.txt", data=jags.data, n.chains = 3,
                         n.adapt=nadapt)
saveRDS(network.mod, "modelruns/icar_sampling_events_network_adapt.rds")
update(network.mod, nadapt)
network.mod.fit <- coda.samples(model=network.mod, variable.names=variable.names, n.iter=niter,
                               thin=thin)
saveRDS(network.mod.fit, "modelruns/icar_sample_events_network_fit.rds")

#plots
setwd("~/HP/Plots")
pdf("iCAR_network_trace_density_sensitivity_specificity.pdf")
plot(network.mod.fit[,1:2])
dev.off()

pdf("iCAR_network_trace_density_sensitivity_specificity.pdf")
plot(network.mod.fit[,1:2])
dev.off()

pdf("iCAR_network_sensitivity_autocorrelation.pdf")
autocorr.plot(network.mod.fit[,'Se'], main="Sensitivity Autocorrelation")
dev.off()

pdf("iCAR_network_specificity_autocorrelation.pdf")
autocorr.plot(network.mod.fit[, 'Sp'], main="Specificity Autocorrelation")
dev.off()

###################################################################################################
#spatial model - weighted network
###################################################################################################

W <- define.neighborhood(method="weightednetwork")
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
wtnetwork.mod <- jags.model(file = "icar_sampling_events.txt", data=jags.data, n.chains = 3,
                         n.adapt=nadapt)
saveRDS(wtnetwork.mod, "modelruns/icar_sampling_events_wtnetwork_adapt.rds")
update(wtnetwork.mod, nadapt)
wtnetwork.mod.fit <- coda.samples(model=wtnetwork.mod, variable.names=variable.names, n.iter=niter,
                               thin=thin)
saveRDS(wtnetwork.mod.fit, "modelruns/icar_sample_events_wtnetwork_fit.rds")

#plots
setwd("~/HP/Plots")
pdf("iCAR_wtnetwork_trace_density_sensitivity_specificity.pdf")
plot(wtnetwork.mod.fit[,1:2])
dev.off()

pdf("iCAR_wtnetwork_trace_density_sensitivity_specificity.pdf")
plot(wtnetwork.mod.fit[,1:2])
dev.off()

pdf("iCAR_wtnetwork_sensitivity_autocorrelation.pdf")
autocorr.plot(wtnetwork.mod.fit[,'Se'], main="Sensitivity Autocorrelation")
dev.off()

pdf("iCAR_wtnetwork_specificity_autocorrelation.pdf")
autocorr.plot(wtnetwork.mod.fit[, 'Sp'], main="Specificity Autocorrelation")
dev.off()

rm(wtnetwork.mod, wtnetwork.mod.fit)
##########################################################################################
#AR1 Model - Temporal Correlation
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

nadapt <- 500
niter <- 10
thin <- 1

setwd("~/Github/AI_surveillance")
ar1.mod <- jags.model(file = "ar1_sampling_events.txt", data=jags.data, n.chains = 3,
                            n.adapt=nadapt)
saveRDS(ar1.mod, "modelruns/icar_sampling_events_ar1_adapt.rds")
update(ar1.mod, nadapt)
ar1.mod.fit <- coda.samples(model=ar1.mod, variable.names=variable.names, n.iter=niter,
                                  thin=thin)
saveRDS(ar1.mod.fit, "modelruns/icar_sample_events_ar1_fit.rds")

#plots
setwd("~/HP/Plots")
pdf("ar1_trace_density_sensitivity_specificity.pdf")
plot(ar1.mod.fit[,1:2])
dev.off()

pdf("ar1_trace_density_sensitivity_specificity.pdf")
plot(ar1.mod.fit[,1:2])
dev.off()

pdf("ar1_sensitivity_autocorrelation.pdf")
autocorr.plot(ar1.mod.fit[,'Se'], main="Sensitivity Autocorrelation")
dev.off()

pdf("ar1_specificity_autocorrelation.pdf")
autocorr.plot(ar1.mod.fit[, 'Sp'], main="Specificity Autocorrelation")
dev.off()

rm(ar1.mod, ar1.mod.fit)

###################################################################################
####################################################################################
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

