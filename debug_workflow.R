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
# original all bird species, reduced AIV dataset with only birds with bands, watersheds with data
n <- readRDS("data_n_red.rds")
y <- readRDS("data_y_red.rds")

# mallards only, full AIV samples including nobands, watersheds with data
malln <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_n_mall_all.rds")
mally <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_y_mall_all.rds")

# mallards only, full AIV samples including nobands, all watersheds
mallnfull <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_n_mall_all_fullhuc.rds")
mallyfull <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_y_mall_all_fullhuc.rds")

# mallards only, seasonal sampling including nobands, watersheds with data
malln2season <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_n_mall2season.rds")
mally2season <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_n_mall2season.rds")
malln3season <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_n_mall3season.rds")
mally3season <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/data_n_mall3season.rds")

nadapt <- 5000
niter <- 5000
thin <- 10

###############################################
###############################################
# Run iCAR model with queens neighborhood for tests
###############################################
###############################################
# define W & D
W <- define.neighborhood(method = "queens")
B <- scaleW(W)
D <- diag(rowSums(W))
I <- diag(1, nrow = dim(W)[1], ncol = dim(W)[1])

# 1) Full data - no inits 
###############################################
ndata <- mallnfull
ydata <- mallyfull

nmonths <- dim(ndata)[1] 
nyears <- dim(ndata)[2]
nsites <- dim(ndata)[3]
npi <- nmonths*nyears*nsites
ndates <- nmonths*nyears

# define data, parameters, initial conditions, iterations, thinning
jags.data <- list(nsites = nsites, nmonths = nmonths, nyears = nyears,
    y = ydata, n = ndata, W = B, I = I)
jags.parameters <- c('Se', 'Sp', 'pi')

# run model
bayes.mod <- jags.model(file = 'icar_constant_sensitivity.txt', data = jags.data, 
   n.chains = 3, n.adapt = nadapt)
saveRDS(bayes.mod, "icarqueens_initialize_5000adapt_10thin_t1.rds")
bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = niter, thin = thin)
saveRDS(bayes.mod.fit, "icarqueens_coda_5000iterations_10thin_t1.rds")
pdf("test1_param1.pdf")
plot(bayes.mod.fit[, 1])
dev.off()
pdf("test1_param2.pdf")
plot(bayes.mod.fit[, 2])
dev.off()
pdf("test1_param3.pdf")
plot(bayes.mod.fit[, 3])
dev.off()
pdf("test1_gelman1.pdf")
gelman.plot(bayes.mod.fit[ ,1])
dev.off()
pdf("test1_gelman2.pdf")
gelman.plot(bayes.mod.fit[ ,2])
dev.off()

rm(bayes.mod) 
rm(bayes.mod.fit)

# 2) Full data - dispersed inits 
###############################################
# define data, parameters, initial conditions, iterations, thinning
jags.data <- list(nsites = nsites, nmonths = nmonths, nyears = nyears,
    y = ydata, n = ndata, W = B, I = I)
jags.inits <- function(){
    list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.6, 1), 
         "lmupi" = array(runif(ndates, 0, 0.4), dim = c(nmonths, nyears)) ) 
}
jags.parameters <- c('Se', 'Sp', 'pi')

# run model
bayes.mod <- jags.model(file = 'icar_constant_sensitivity.txt', data = jags.data, 
    inits = jags.inits, n.chains = 3, n.adapt = nadapt)
saveRDS(bayes.mod, "icarqueens_initialize_5000adapt_10thin_t2.rds")
bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = niter, thin = thin)
saveRDS(bayes.mod.fit, "icarqueens_coda_5000iterations_10thin_t2.rds")
pdf("test2_param1.pdf")
plot(bayes.mod.fit[, 1])
dev.off()
pdf("test2_param2.pdf")
plot(bayes.mod.fit[, 2])
dev.off()
pdf("test2_param3.pdf")
plot(bayes.mod.fit[, 3])
dev.off()
pdf("test2_gelman1.pdf")
gelman.plot(bayes.mod.fit[ ,1])
dev.off()
pdf("test2_gelman2.pdf")
gelman.plot(bayes.mod.fit[ ,2])
dev.off()

rm(bayes.mod) 
rm(bayes.mod.fit)

# 3) Reduced data - no inits 
###############################################
ndata <- malln
ydata <- mally

subsites <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/subset_sites.rds")
Bred <- B[rownames(B) %in% subsites, rownames(B) %in% subsites]
Ired <- diag(1, nrow = dim(Bred)[1], ncol = dim(Bred)[1])

nmonths <- dim(ndata)[1] 
nyears <- dim(ndata)[2]
nsites <- dim(ndata)[3]
npi <- nmonths*nyears*nsites
ndates <- nmonths*nyears

# define data, parameters, initial conditions, iterations, thinning
jags.data <- list(nsites = nsites, nmonths = nmonths, nyears = nyears,
    y = ydata, n = ndata, W = B, I = I)
jags.parameters <- c('Se', 'Sp', 'pi')

# run model
bayes.mod <- jags.model(file = 'icar_constant_sensitivity.txt', data = jags.data, 
    n.chains = 3, n.adapt = nadapt)
saveRDS(bayes.mod, "icarqueens_initialize_5000adapt_10thin_t3.rds")
bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = niter, thin = thin)
saveRDS(bayes.mod.fit, "icarqueens_coda_5000iterations_10thin_t3.rds")

pdf("test3_param1.pdf")
plot(bayes.mod.fit[, 1])
dev.off()
pdf("test3_param2.pdf")
plot(bayes.mod.fit[, 2])
dev.off()
pdf("test3_param3.pdf")
plot(bayes.mod.fit[, 3])
dev.off
pdf("test3_gelman1.pdf")
gelman.plot(bayes.mod.fit[ ,1])
dev.off()
pdf("test3_gelman2.pdf")
gelman.plot(bayes.mod.fit[ ,2])
dev.off()

rm(bayes.mod) 
rm(bayes.mod.fit)

# 4) Reduced data - dispersed inits 
###############################################

# define data, parameters, initial conditions, iterations, thinning
jags.data <- list(nsites = nsites, nmonths = nmonths, nyears = nyears,
    y = ydata, n = ndata, W = Bred, I = Ired)
jags.inits <- function(){
    list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.6, 1), 
        "lmupi" = array(runif(ndates, 0, 0.4), dim = c(nmonths, nyears)) ) 
}
jags.parameters <- c('Se', 'Sp', 'pi')

# run model
bayes.mod <- jags.model(file = 'icar_constant_sensitivity.txt', data = jags.data, 
    inits = jags.inits, n.chains = 3, n.adapt = nadapt)
saveRDS(bayes.mod, "icarqueens_initialize_5000adapt_10thin_t4.rds")
bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = niter, thin = thin)
saveRDS(bayes.mod.fit, "icarqueens_coda_5000iterations_10thin_t4.rds")

pdf("test4_param1.pdf")
plot(bayes.mod.fit[, 1])
dev.off()
pdf("test4_param2.pdf")
plot(bayes.mod.fit[, 2])
dev.off()
pdf("test4_param3.pdf")
plot(bayes.mod.fit[, 3])
dev.off()
pdf("test4_gelman1.pdf")
gelman.plot(bayes.mod.fit[ ,1])
dev.off()
pdf("test4_gelman2.pdf")
gelman.plot(bayes.mod.fit[ ,2])
dev.off()

rm(bayes.mod) 
rm(bayes.mod.fit)


###############################################
###############################################
# Run independence model
###############################################
###############################################
# reduced data - 181 hucs
ndata <- malln
ydata <- mally
nmonths <- 12
nsites <- dim(ndata)[3]
nyears <- dim(ndata)[2]
npi <- nsites*nyears*nmonths

# define data, parameters, initial conditions, iterations, thinning
jags.data <- list(nsites = nsites, nmonths = nmonths, nyears = nyears,
    y = ydata, n = ndata)
jags.inits <- function(){
    list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.9, 1), 
    "pi" = array(runif(npi, 0, 0.4), dim = c(nmonths, nyears, nsites)) ) 
}
jags.parameters <- c('Se', 'Sp', 'pi')

# run model
bayes.mod <- jags.model(file = 'constant_sensitivity.txt', data = jags.data, 
    inits = jags.inits, n.chains = 3, n.adapt = nadapt)
saveRDS(bayes.mod, "indep_initialize_1000adapt_10thin_t5.rds")
bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = niter, thin = thin)
saveRDS(bayes.mod.fit, "indep_coda_5000iterations_10thin_t5.rds")

pdf("test5_param1.pdf")
plot(bayes.mod.fit[, 1])
dev.off()
pdf("test5_param2.pdf")
plot(bayes.mod.fit[, 2])
dev.off()
pdf("test5_param3.pdf")
plot(bayes.mod.fit[, 3])
dev.off()
pdf("test5_gelman1.pdf")
gelman.plot(bayes.mod.fit[ ,1])
dev.off()
pdf("test5_gelman2.pdf")
gelman.plot(bayes.mod.fit[ ,2])
dev.off()

rm(bayes.mod) 
rm(bayes.mod.fit)
pd <- dic.samples(model = bayes.mod, n.iter = 50000, thin = thin, type="pD")
saveRDS(pd, "base_dic_mallard_monthly.rds", sep="")

###############################################
###############################################
# Run CAR model with weighted queens neighborhood
###############################################
###############################################
# define W & D
W <- define.neighborhood(method = "weightednetwork")
B <- scaleW(W)
D <- diag(rowSums(W))
I <- diag(1, nrow = dim(W)[1], ncol = dim(W)[1])

subsites <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/subset_sites.rds")
Bred <- B[rownames(B) %in% subsites, rownames(B) %in% subsites]
Ired <- diag(1, nrow = dim(Bred)[1], ncol = dim(Bred)[1])

# 1) Full data - no inits 
###############################################
nmonths <- dim(ndata)[1] 
nyears <- dim(ndata)[2]
nsites <- dim(ndata)[3]
npi <- nmonths*nyears*nsites
ndates <- nmonths*nyears

# define data, parameters, initial conditions, iterations, thinning
jags.data <- list(nsites = nsites, nmonths = nmonths, nyears = nyears,
    y = ydata, n = ndata, W = Bred, I = Ired)
jags.inits <- function(){
    list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.9, 1), 
    "lmupi" = array(runif(ndates, 0, 0.4), dim = c(nmonths, nyears)) ) 
}
jags.parameters <- c('Se', 'Sp', 'pi')

# quick run model
bayes.mod <- jags.model(file = 'car_constant_sensitivity.txt', data = jags.data, 
    inits = jags.inits, n.chains = 3, n.adapt = nadapt)
bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = 500, thin = thin)
plot(bayes.mod.fit[, 1:2])

# DIC calculation
nadapt <- 10000
bayes.mod <- jags.model(file = 'car_constant_sensitivity.txt', data = jags.data, 
    inits = jags.inits, n.chains = 3, n.adapt = nadapt)
update(bayes.mod, 10000)
bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = 5000, thin = thin)

pdf("test6_param1.pdf")
plot(bayes.mod.fit[, 1])
dev.off()
pdf("test6_param2.pdf")
plot(bayes.mod.fit[, 2])
dev.off()
pdf("test6_param3.pdf")
plot(bayes.mod.fit[, 3])
dev.off()
pdf("test6_gelman1.pdf")
gelman.plot(bayes.mod.fit[ ,1])
dev.off()
pdf("test6_gelman2.pdf")
gelman.plot(bayes.mod.fit[ ,2])
dev.off()

pd <- dic.samples(model = bayes.mod, n.iter = 10000, thin = thin, type="pD")
saveRDS(pd, "dic_mallard_monthly_wtnetwork.rds")
save(bayes.mod, "burnin_mallard_monthly_wtnetwork.rds")

###############################################
###############################################
# Run iCAR model with network neighborhood for tests
###############################################
###############################################
# reduced data - 181 hucs
ndata <- malln
ydata <- mally

# define W & D
W <- define.neighborhood(method = "network")
B <- scaleW(W)
D <- diag(rowSums(W))
I <- diag(1, nrow = dim(W)[1], ncol = dim(W)[1])

subsites <- readRDS("~/Documents/Avian_Influenza/HP/statistical model/subset_sites.rds")
Bred <- B[rownames(B) %in% subsites, rownames(B) %in% subsites]
Ired <- diag(1, nrow = dim(Bred)[1], ncol = dim(Bred)[1])

nmonths <- dim(ndata)[1] 
nyears <- dim(ndata)[2]
nsites <- dim(ndata)[3]
npi <- nmonths*nyears*nsites
ndates <- nmonths*nyears

# define data, parameters, initial conditions, iterations, thinning
jags.data <- list(nsites = nsites, nmonths = nmonths, nyears = nyears,
    y = ydata, n = ndata, W = Bred, I = Ired)
jags.inits <- function(){
    list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.9, 1), 
    "lmupi" = array(runif(ndates, 0, 0.4), dim = c(nmonths, nyears)) ) 
}
jags.parameters <- c('Se', 'Sp', 'pi')

# run model
bayes.mod <- jags.model(file = 'icar_constant_sensitivity.txt', data = jags.data, 
    inits = jags.inits, n.chains = 3, n.adapt = 50000)
saveRDS(bayes.mod, "icarqueens_initialize_5000adapt_10thin_t7.rds")
bayes.mod.fit <- coda.samples(model = bayes.mod,
    variable.names = jags.parameters, n.iter = 50000, thin = thin)
saveRDS(bayes.mod.fit, "icarqueens_coda_5000iterations_10thin_t7.rds")

pdf("test7_param1.pdf")
plot(bayes.mod.fit[, 1])
dev.off()
pdf("test7_param2.pdf")
plot(bayes.mod.fit[, 2])
dev.off()
pdf("test7_param3.pdf")
plot(bayes.mod.fit[, 3])
dev.off()
pdf("test7_gelman1.pdf")
gelman.plot(bayes.mod.fit[ ,1])
dev.off()
pdf("test7_gelman2.pdf")
gelman.plot(bayes.mod.fit[ ,2])
dev.off()

rm(bayes.mod) 
rm(bayes.mod.fit)
