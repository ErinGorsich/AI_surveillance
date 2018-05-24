library(rjags)
library(coda)

setwd("~/GitHub/HP_AI_Surveillance")

# load script holding the text files that define the model
source('define.models.R', chdir = TRUE)

# load script holding the function that defines the neighborhoods
source('define.neighborhood.R', chdir = TRUE)

###############################################
###############################################
# read in data 
###############################################
###############################################
# mallards only, seasonal sampling including nobands, watersheds with data
malln2season <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_n_mall2season.rds")
mally2season <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_n_mall2season.rds")
malln3season <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_n_mall3season.rds")
mally3season <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_n_mall3season.rds")

nadapt <- 5000
niter <- 5000
thin <- 10

nseasons <- dim(malln3season)[1] 
nyears <- dim(malln3season)[2]
nsites <- dim(malln3season)[3]
npi <- nseasons*nyears*nsites
ndates <- nseasons*nyears


###############################################
###############################################
# Run independence model - 3 seasons
###############################################
###############################################
# reduced data - 181 hucs
ndata <- malln3season
ydata <- malln3season

# define data, parameters, initial conditions, iterations, thinning
jags.data <- list(nsites = nsites, nmonths = nseasons, nyears = nyears,
    y = ydata, n = ndata)
jags.inits <- function(){
    list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.9, 1), 
        "pi" = array(runif(npi, 0, 0.4), dim = c(nseasons, nyears, nsites)) ) 
}
jags.parameters <- c('Se', 'Sp', 'pi')

# run model
bayes.mod <- jags.model(file = 'constant_sensitivity.txt', data = jags.data, 
    inits = jags.inits, n.chains = 3, n.adapt = nadapt)
saveRDS(bayes.mod, "indep_initialize_5000adapt_10thin_s1.rds")
bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = niter, thin = thin)
saveRDS(bayes.mod.fit, "indep_coda_5000iterations_10thin_s1.rds")

pdf("s1_param1.pdf")
plot(bayes.mod.fit[, 1])
dev.off()
pdf("s1_param2.pdf")
plot(bayes.mod.fit[, 2])
dev.off()
pdf("s1_param3.pdf")
plot(bayes.mod.fit[, 3])
dev.off()
pdf("s1_gelman1.pdf")
gelman.plot(bayes.mod.fit[ ,1])
dev.off()
pdf("s1_gelman2.pdf")
gelman.plot(bayes.mod.fit[ ,2])
dev.off()

rm(bayes.mod) 
rm(bayes.mod.fit)

# run base default inits
bayes.mod <- jags.model(file = 'constant_sensitivity.txt', data = jags.data, 
    n.chains = 3, n.adapt = nadapt)
saveRDS(bayes.mod, "indep_initialize_5000adapt_10thin_s2.rds")
bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = niter, thin = thin)
saveRDS(bayes.mod.fit, "indep_coda_5000iterations_10thin_s2.rds")

pdf("s2_param1.pdf")
plot(bayes.mod.fit[, 1])
dev.off()
pdf("s2_param2.pdf")
plot(bayes.mod.fit[, 2])
dev.off()
pdf("s2_param3.pdf")
plot(bayes.mod.fit[, 3])
dev.off()
pdf("s2_gelman1.pdf")
gelman.plot(bayes.mod.fit[ ,1])
dev.off()
pdf("s2_gelman2.pdf")
gelman.plot(bayes.mod.fit[ ,2])
dev.off()

rm(bayes.mod) 
rm(bayes.mod.fit)


###############################################
###############################################
# Run independence model - 2 seasons
###############################################
###############################################
# reduced data - 181 hucs
ndata <- malln2season
ydata <- malln2season
nseasons <- 2
npi <- nseasons*nyears*nsites
    
# define data, parameters, initial conditions, iterations, thinning
jags.data <- list(nsites = nsites, nmonths = nseasons, nyears = nyears,
    y = ydata, n = ndata)
jags.inits <- function(){
    list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.9, 1), 
        "pi" = array(runif(npi, 0, 0.4), dim = c(nseasons, nyears, nsites)) ) 
}
jags.parameters <- c('Se', 'Sp', 'pi')

# run model
bayes.mod <- jags.model(file = 'constant_sensitivity.txt', data = jags.data, 
                        inits = jags.inits, n.chains = 3, n.adapt = nadapt)
saveRDS(bayes.mod, "indep_initialize_5000adapt_10thin_s3.rds")
bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = niter, thin = thin)
saveRDS(bayes.mod.fit, "indep_coda_5000iterations_10thin_s3.rds")

pdf("s3_param1.pdf")
plot(bayes.mod.fit[, 1])
dev.off()
pdf("s3_param2.pdf")
plot(bayes.mod.fit[, 2])
dev.off()
pdf("s3_param3.pdf")
plot(bayes.mod.fit[, 3])
dev.off()
pdf("s3_gelman1.pdf")
gelman.plot(bayes.mod.fit[ ,1])
dev.off()
pdf("s3_gelman2.pdf")
gelman.plot(bayes.mod.fit[ ,2])
dev.off()

rm(bayes.mod) 
rm(bayes.mod.fit)

# no inits
# run model
bayes.mod <- jags.model(file = 'constant_sensitivity.txt', data = jags.data, 
                        inits = jags.inits, n.chains = 3, n.adapt = nadapt)
saveRDS(bayes.mod, "indep_initialize_5000adapt_10thin_s3.rds")
bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = niter, thin = thin)
saveRDS(bayes.mod.fit, "indep_coda_5000iterations_10thin_s3.rds")

pdf("s4_param1.pdf")
plot(bayes.mod.fit[, 1])
dev.off()
pdf("s4_param2.pdf")
plot(bayes.mod.fit[, 2])
dev.off()
pdf("s4_param3.pdf")
plot(bayes.mod.fit[, 3])
dev.off()
pdf("s4_gelman1.pdf")
gelman.plot(bayes.mod.fit[ ,1])
dev.off()
pdf("s4_gelman2.pdf")
gelman.plot(bayes.mod.fit[ ,2])
dev.off()

rm(bayes.mod) 
rm(bayes.mod.fit)



