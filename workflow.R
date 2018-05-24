library(rjags)
library(coda)

setwd("~/GitHub/HP_AI_Surveillance")

# load script holding the text files that define the model
source('define.models.R', chdir = TRUE)

# load script holding the function that defines the neighborhoods
source('define.neighborhood.R', chdir = TRUE)

# load script holding the function to run the model and save results to file
#source('~/Honors Thesis/Thesis Work/Model Selection/run_jags_model.R', chdir = TRUE)

# load other scripts needed to run this file
# ??? diagnostics?


# read in data (base model data- 136 watersheds)
#####################################
n <- readRDS("data_n_red.rds")
y <- readRDS("data_y_red.rds")

nmonths <- dim(n)[1]
nyears <- dim(n)[2]
nsites <- dim(n)[3]
npi <- nmonths*nyears*nsites

nadapt <- 5000
niter <- 50000
thin <- 10


#####################################
#####################################
# Run base model
#####################################
#####################################
# define data, parameters, initial conditions, iterations, thinning
jags.data <- list(nsites = nsites, nmonths = nmonths, nyears = nyears, 
    y = y, n = n)
jags.inits <- function(){
    list("Se" = runif(1, 0.65, 1), "Sp" = runif(1, 0.65, 1), 
         "pi" = array(runif(npi), dim = c(nmonths, nyears, nsites)) )
}
jags.inits2 <- function(){
    list("Se" = runif(1, 0.65, 1), "Sp" = runif(1, 0.65, 1), 
    "pi" = array(runif(npi, 0, 0.5), dim = c(nmonths, nyears, nsites)) )
}
jags.parameters <- c('Se', 'Sp', 'pi')

# run model
bayes.mod <- jags.model(file = 'constant_sensitivity.txt', data = jags.data, 
    inits = jags.inits2, 
    n.chains = 3, n.adapt = nadapt)
#saveRDS(bayes.mod, "nospace_initialize_5000adapt_10thin_redprevinits.rds")

bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = 1000, thin = thin)
saveRDS(bayes.mod.fit, "nospace_coda_50000iterations_10thin_redprevinits.rds")

rm(bayes.mod, bayes.mod.fit)




# read in data (spatial model data- 202 watersheds)
#####################################
n <- readRDS("data_n-2016.rds"); dim(n)
y <- readRDS("data_y-2016.rds"); dim(n) # should be 12,6,202

nmonths <- dim(n)[1]
nyears <- dim(n)[2]
nsites <- dim(n)[3]
npi <- nmonths*nyears*nsites
ndates <- nmonths*nyears

nadapt <- 10000
niter <- 50000
thin <- 10

#####################################
#####################################
# Run CAR model wiht queens neighborhood 
#####################################
#####################################
# define W & D
W <- define.neighborhood(method = "queens")
B <- scaleW(W)
D <- diag(rowSums(W))
I <- diag(1, nrow = dim(W)[1], ncol = dim(W)[1])

# define data, parameters, initial conditions, iterations, thinning
jags.data <- list(nsites = nsites, nmonths = nmonths, nyears = nyears,
    y = y, n = n, W = B, I = I)
jags.inits <- function(){
    list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.6, 1), 
         "lmupi" = array(runif(ndates), dim = c(nmonths, nyears)) )  # changed no more pi0
}

# run model
bayes.mod <- jags.model(file = 'icar_constant_sensitivity.txt', data = jags.data, 
    inits = jags.inits, n.chains = 3, n.adapt = nadapt)
saveRDS(bayes.mod, "carqueens_initialize_10000adapt_10thin.rds")
bayes.mod.fit <- coda.samples(model = bayes.mod, 
    variable.names = jags.parameters, n.iter = niter, thin = thin)
saveRDS(bayes.mod.fit, "carqueens_coda_50000iterations_10thin.rds")


# remove model
rm(bayes.mod, bayes.mod.fit)

#####################################
#####################################
# Run CAR model with weighted queens neighborhood 
#####################################
#####################################
# define W & D
W <- define.neighborhood(method = "weightedqueens")
D <- diag(rowSums(W))

# define data, parameters, initial conditions, iterations, thinning
jags.data <- list(nsites = nsites, nmonths = nmonths, nyears = nyears,
                  y = y, n = n, W = W, D = D)

# inits and parameters should be the same as above

# run model
testrun <- run_jags_model(textfile = 'car_constant_sensitivity.txt', 
    data = jags.data, inits = jags.inits, n.chains = 3, n.adapt = 50, 
    variable.names = jags.parameters, n.iter = 10, thin = 1, 
    name ="car_weighted_queens")

# fill in df
testrun
#df[2, ] <- c("car_weighted_queens", 0, 0, 0)

# check diagnostics later by reading back in the coda files. 
# maybe make another function to make thee figures
jags.mod <- readRDS("car_weighted_queens.rds")

# remove model
rm(jags.mod, testrun)

#####################################
#####################################
# Run CAR model with knn3 neighborhood 
#####################################
#####################################


#####################################
#####################################
# Run CAR model with network-based neighborhood 
#####################################
#####################################





#####################################
#####################################
# weighted network
#####################################
#####################################
#define a W matrix
w <- define.neighborhood("weightednetwork")
B <- scaleW(w); B[1:10, 1:10]
D <- diag(1, nsites, nsites); D[1:10, 1:10]

jags.data <- list(nsites=nsites, nmonths=nmonths, nyears=nyears, y=y, 
                  n=n, W=B, D=D)
jags.inits <- function(){
    list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.6, 1), 
         "pi0" = array(runif(npi), dim = c(nmonths, nyears)), 
         "phi" = runif(0, 1))
}

#run CAR model with weighted queens DIC
net.mod <- jags.model(file="car_constant_sensitivity.txt", data = jags.data, #inits = jags.inits,
                      n.chains=3, n.adapt=nadapt)
saveRDS(net.mod, "CAR_weightednetwork_adaption.rds")
pd.net <- dic.samples(model=net.mod, n.iter=niter, thin=thin, type='pD')
summary(pd.net)
saveRDS(pd.net, "CAR_weightednetwork_DIC.rds")  # updated name save!

rm(net.mod, pd.net)









