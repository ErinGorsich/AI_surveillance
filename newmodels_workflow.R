####################################################################################
#workflow for updated models
###################################################################################

#load packages
library(rjags)
library(coda)

#load models
setwd("~/Github/AI_surveillance")
source("~/Github/AI_surveillance/define_models.r")

#read in sampling data
data <- readRDS('~/Github/samplingevent_n_y_speciesgroup.rds')

#define variables
nsamplingevents <- length(unique(data$sample.event)) #13638
nspecies <- length(unique(data$species.group)) #7
nmonths <- length(rep(seq(1, 12), times = 5)) #60
nyears <- length(unique(data$year)) #5
nhucs <- length(unique(data$watershed)) #206

#######################################################################################
#base model
#########################################################################################

jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, nmonths=nmonths,
                  nyears=nyears, nhucs = nhucs)
#need initial values for lambda
jags.inits <- function(){
  list("Se" = runif(1, 0.6, 1), 'Sp' = runif(1, 0.6, 1), 
       "pi" = array(runif(nsamplingevents*nmonths*nyears), dim = c(nmonths, nyears, nsites)))
}
variable.names = c('Se', 'Sp', 'lambda', 'pi') #apparent prevalence?

nadapt <- 1000
niter <- 10
thin <- 1

setwd("~/Github/AI_surveillance")
base.mod <- jags.model(file = "base_sampling_events.txt", data=jags.data, inits=jags.inits,
                      n.chains=3, n.adapt=nadapt)
saveRDS(base.mod, "modelruns/base_sampling_events_adapt.rds")
base.mod.fit <- coda.samples(model=base.mod, variable.names=variable.names, n.iter=niter, 
                             thin=thin)
saveRDS(base.mod.fit, "modelruns/base_sampling_event_fit.rds")


#alpha line 17 -  needs iniital value/definition