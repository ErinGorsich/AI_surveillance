###########################################################################################
###########################################################################################
#Attempt to run jags in parallel 
#using this tutorial:
#https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/585b0e4c/
###########################################################################################
###########################################################################################

library(snow)
library(rjags)
library(coda)

# setwd("~/Github/AI_surveillance")
# source("define_models.rds")

##########################################################################################
#data
##########################################################################################

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
nspecies <- 2  #length(unique(data$species.group)) #7
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

##########################################################################################
#run the model
##########################################################################################
month = data$month[1:nsamplingevents]
year = data$yearid[1:nsamplingevents]
huc = data$hucid[1:nsamplingevents]

jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, 
                  nmonths = nmonths, nyears = nyears, nhucs = nhucs, n = n, y = y, 
                  month = month, year = year, huc = huc)

variable.names = c("Se", "Sp", "pi")

nchains = 3

coda.samples.wrapper <- function(j)
{
  jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, 
                    nmonths = nmonths, nyears = nyears, nhucs = nhucs, n = n, y = y, 
                    month = month, year = year, huc = huc)
  variable.names = c("Se", 'Sp', "pi")
  nadapt <- 1000
  niter <- 10
  thin <- 2
  nchains <- 3
  setwd("~/Github/AI_surveillance")
  base.parallel.mod <- jags.model("base_sampling_events.txt", data = jags.data,
                                  n.chains=nchains, n.adapt = nadapt)
  
  coda.samples(base.parallel.mod, variable.names=variable.names, n.iter=niter, thin=thin)
}

parallel.start.time = proc.time()
  chain.1 <- makeCluster(nchains, "SOCK")
    clusterEvalQ(chain.1, library(rjags))
      clusterExport(chain.1,
                    list("nsamplingevents", "nspecies", "nmonths", "nyears", "nhucs", "n", 'y',
                          'month', "year", "huc", "niter", "thin"))
    par.samples <- clusterApply(chain.1, 1:nchains, coda.samples.wrapper)
    
    for(i in 1:length(par.samples)) {
      par.samples[[i]] <- par.samples[[i]][[1]]
    }
    class(par.samples) <- "mcmc.list"
  stopCluster(chain.1)
parallel.end.time = proc.time()
parallel.total.time = parallel.end.time - parallel.start.time


######################################################################################
######################################################################################
#Attempt to run jags in parallel using this tutuorial:
#https://stephendavidgregory.github.io/statistics/Jags-in-parallel
#and the dclone package with jags.parfit
######################################################################################
######################################################################################

#data prepped as before

month = data$month[1:nsamplingevents]
year = data$yearid[1:nsamplingevents]
huc = data$hucid[1:nsamplingevents]

jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, 
                  nmonths = nmonths, nyears = nyears, nhucs = nhucs, n = n, y = y, 
                  month = month, year = year, huc = huc)

variable.names = c("Se", "Sp", "pi")

nchains = 3

timer <- proc.time()
cl <- makePSOCKcluster(nchains)
tmp <- clusterEvalQ(cl, library(dclone))
setwd("~/Github/AI_surveillance")
fit <- jags.parfit(cl = cl, data = jags.data, params = variable.names,
                   model = "~/Github/AI_surveillance/base_sampling_events.txt", 
                   n.chains = 3,
                   n.adapt = 1000, n.update = 100, n.iter = 10, thin = 1)
stopCluster(cl)
time.taken <- proc.time() - timer

######################################################################################
######################################################################################
#Attempt to run jags in parallel using this tutuorial:
#https://stephendavidgregory.github.io/statistics/Jags-in-parallel
#and the dclone package with parJagsModel and parCodasamples
######################################################################################
######################################################################################

#data prepped as before

month = data$month[1:nsamplingevents]
year = data$yearid[1:nsamplingevents]
huc = data$hucid[1:nsamplingevents]

jags.data <- list(nsamplingevents = nsamplingevents, nspecies = nspecies, 
                  nmonths = nmonths, nyears = nyears, nhucs = nhucs, n = n, y = y, 
                  month = month, year = year, huc = huc)

variable.names = c("Se", "Sp", "pi")

nchains = 3

timer <- proc.time()
cl <- makePSOCKcluster(nchains)
parJagsModel(cl = cl, name = 'res', file = "~/Github/AI_surveillance/base_sampling_events.txt",
             data = jags.data, n.chains=nchains, n.adapt = 1000)
parUpdate(cl=cl, object = 'res', n.iter = 100)
fit <- parCodaSamples(cl = cl, model = "res", variable.names = variable.names,
                      n.iter = 10, thin = 1)
time.taken <- proc.time() - timer
