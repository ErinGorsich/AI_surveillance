setwd('/home/webblab/Documents/Github/AI_surveillance')

# load models
source('define_models.R')
source('define_neighborhood.R')
# source('model selection function.R')

#packages
library(rjags)
library(coda)

# #base model data read in 
# setwd("~/Honors Thesis/Thesis Work/JAGS")
# n <- readRDS("n matrix.rds")
# y <- readRDS("y matrix.rds")
# nsites <- dim(y)[3]
# nyears <- dim(y)[2]
# nmonths <- dim(y)[1]

#read in n and y matrices with the data from all hucs
n <- readRDS('/home/webblab/Documents/HP/sim_data_6_n_short.rds')
y <- readRDS('/home/webblab/Documents/HP/sim_data_6_y_short.rds')

#define nsites to match new y matrix
nsites <- dim(y)[3]
nyears <- dim(y)[2]
nmonths <- dim(y)[1]
npi <- nsites*nyears*nmonths

#base model, no W matrix needed
jags.data <- list(nsites = nsites, nmonths=nmonths, nyears=nyears, y=y, n=n)
set.seed(100)
jags.inits <- function(){
  list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.6, 1),
       "pi" = array(runif(nsites*nmonths*nyears), dim = c(nmonths, nyears, nsites)))
}
variable.names = c('Se', 'Sp', 'pi')
nadapt = 30000
niter = 10000
thin = 10
setwd("/home/webblab/Documents/Github/AI_surveillance")
base.mod <- jags.model(file = "constant_sensitivity.txt", data=jags.data, inits=jags.inits,
                      n.chains=3, n.adapt = nadapt)
base.mod.fit <- coda.samples(model=base.mod, variable.names = variable.names,
                            n.iter=niter, thin=thin)

setwd("/home/webblab/Documents/HP/simulated_results")
saveRDS(base.mod.fit, 'base_model_fit_sim_6_short20.rds')

pdf("base_trace_density_sensitivity_specificity_mall_attempt6_short20.pdf")
plot(base.mod.fit[,1:2])
dev.off()

pdf("base_sensitivity_autocorrelation_mall_attempt6_short20.pdf")
autocorr.plot(base.mod.fit[,'Se'], main="Sensitivity Autocorrelation")
dev.off() 

pdf("base_specificity_autocorrelation_mall_attempt6_short20.pdf")
autocorr.plot(base.mod.fit[, 'Sp'], main="Specificity Autocorrelation")
dev.off()

#define data for running the AR1 model, no W matrix
# jags.data <- list(nsites=nsites, nmonths=nmonths, nyears=nyears, y=y, n=n)
# jags.inits <- function(){
#   list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.6, 1),
#        "pi0" = array(runif(npi), dim = c(nmonths, nyears)), 
#        "phi" = runif(0, 1))  # tau?
# }
# 
# #run AR1 model
# model.selection(textfile = "ar1_constant_sensitivity.txt", data=jags.data, inits=jags.inits,
#                 variable.names=variable.names, n.chains = 3, n.adapt=nadapt, n.iter=niter,
#                 thin=thin, name="AR1")

#make a queens W matrix
setwd("/home/webblab/Documents/Github/AI_surveillance")
W <- define.neighborhood(method="queens")
B <- scaleW(W)
D <- diag(rowSums(W))
I <- diag(1, nrow = dim(W)[1], ncol = dim(W)[1])

#define data for running JAGS models
jags.data <- list(nsites=nsites, nmonths=nmonths, nyears=nyears, y=y, n=n, W=W, I=I)
jags.inits <- function(){
  list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.6, 1), 
       "pi0" = array(runif(npi), dim = c(nmonths, nyears)),
       "phi" = runif(0, 1))  # tau?
}

#run CAR model
# model.selection(textfile = "car_constant_sensitivity.txt", data=jags.data, inits=jags.inits,
#                 variable.names=variable.names, n.chains = 3, n.adapt = nadapt, n.iter=niter,
#                 thin=thin, name="CAR_queens")
#run CAR model without the model selection function
car.mod <- jags.model(file = "icar_constant_sensitivity.txt", data=jags.data, inits=jags.inits,
                      n.chains=3, n.adapt = nadapt)
car.mod.fit <- coda.samples(model=car.mod, variable.names = variable.names,
                            n.iter=niter, thin=thin)
setwd("/home/webblab/Documents/HP/simulated_results")
# saveRDS(car.mod.fit, "CAR_queens.rds")

# pd  <- dic.samples(model=car.mod, n.iter = n.iter, thin=thin, type="pD")
# saveRDS(pd, "CAR_queens_dic.rds")

pdf("CAR_queens_trace_density_sensitivity_specificity_attempt2.pdf")
plot(car.mod.fit[,1:2])
dev.off()

pdf("CAR_queens_sensitivity_autocorrelation_attempt2.pdf")
autocorr.plot(car.mod.fit[,'Se'], main="Sensitivity Autocorrelation")
dev.off()

pdf("CAR_queens_specificity_autocorrelation_attempt2.pdf")
autocorr.plot(car.mod.fit[, 'Sp'], main="Specificity Autocorrelation")
dev.off()

#run CAR &AR1 model
# model.selection(textfile = "ar1_car_constant_sensitivity.txt", data=jags.data, inits=jags.inits,
#                 variable.names = variable.names, n.chains=3, n.adapt=nadapt, n.iterniter, 
#                 thin=thin, name="CAR&AR1_queens")

#define a weighted queens matrix
w <- define.neighborhood("weightedqueens")
D <- diag(rowSums(w))

#define data for running JAGS models
jags.data <- list(nsites=nsites, nmonths=nmonths, nyears=nyears, y=y, n=n, W=w, D=D)
jags.inits <- function(){
  list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.6, 1), 
       "pi0" = array(runif(npi), dim = c(nmonths, nyears)), 
       "phi" = runif(0, 1))  # tau?
}

#run CAR model
model.selection(textfile = "car_constant_sensitivity.txt", data=jags.data, inits=jags.inits,
                variable.names=variable.names, n.chains = 3, n.adapt = nadapt, n.iter=niter,
                thin=thin, name="CAR_wqueens")

#run CAR &AR1 model
# model.selection(textfile = "ar1_car_constant_sensitivity.txt", data=jags.data, inits=jags.inits,
#                 variable.names = variable.names, n.chains=3, n.adapt=nadapt, n.iterniter, 
#                 thin=thin, name="CAR&AR1_wqueens")

#define a knn3 W matrix
w <- define.neighborhood("knn3")
D <- diag(rowSums(w))

#define data for running JAGS models
jags.data <- list(nsites=nsites, nmonths=nmonths, nyears=nyears, y=y, n=n, W=w, D=D)
jags.inits <- function(){
  list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.6, 1), 
       "pi0" = array(runif(npi), dim = c(nmonths, nyears)), 
       "phi" = runif(0, 1))  # tau?
}

#run CAR model
model.selection(textfile = "car_constant_sensitivity.txt", data=jags.data, inits=jags.inits,
                variable.names=variable.names, n.chains = 3, n.adapt = nadapt, n.iter=niter,
                thin=thin, name="CAR_knn3")

#run CAR &AR1 model
# model.selection(textfile = "ar1_car_constant_sensitivity.txt", data=jags.data, inits=jags.inits,
#                 variable.names = variable.names, n.chains=3, n.adapt=nadapt, n.iterniter, 
#                 thin=thin, name="CAR&AR1_knn3")

#define a network W matrix
w <- define.neighborhood("network")
D <- diag(rowSums(w))

#define data for running JAGS models
jags.data <- list(nsites=nsites, nmonths=nmonths, nyears=nyears, y=y, n=n, W=w, D=D)
jags.inits <- function(){
  list("Se" = runif(1, 0.6, 1), "Sp" = runif(1, 0.6, 1), 
       "pi0" = array(runif(npi), dim = c(nmonths, nyears)), 
       "phi" = runif(0, 1))  # tau?
}

#run CAR model
model.selection(textfile = "car_constant_sensitivity.txt", data=jags.data, inits=jags.inits,
                variable.names=variable.names, n.chains = 3, n.adapt = nadapt, n.iter=niter,
                thin=thin, name="CAR_network")

# #run CAR &AR1 model
# model.selection(textfile = "ar1_car_constant_sensitivity.txt", data=jags.data, inits=jags.inits,
#                 variable.names = variable.names, n.chains=3, n.adapt=nadapt, n.iterniter, 
#                 thin=thin, name="CAR&AR1_network")

#make a dic table
setwd("~/Honors Thesis/Thesis Work/Model Selection/Plots")

dic.1 <- readRDS("base_dic.rds")
#dic.2 <- readRDS("AR1_dic.rds")
dic.3 <- readRDS("CAR_queens_dic.rds")
#dic.4 <- readRDS("CAR&AR1_queens_dic.rds")
dic.5<- readRDS("CAR_wqueens_dic.rds")
#dic.6 <- readRDS("CAR&AR1_wqueens_dic.rds")
dic.7 <- readRDS("CAR_knn3_dic.rds")
#dic.8 <- readRDS("CAR&AR1_knn3_dic.rds")
dic.9 <- readRDS("CAR_network_dic.rds")
#dic.10 <- readRDS("CAR&AR1_network_dic.rds")
print(c(dic.1, dic.3, dic.5, dic.7, dic.9))
#print(c(dic.1, dic.2, dic.3, dic.4, dic.5, dic.6, dic.7, dic.8, dic.9, dic.10))

