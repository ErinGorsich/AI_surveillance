# Define a sequence of JAGS models for model selection

################################################################
################################################################
# Post Thesis models accounting for sampling events
################################################################
################################################################
# data needed:
## month, year, huc = a vector with sampling events as rows; same for all species
## nsamplingevents, nspecies, nmonths, nyears, nhucs

# see example: http://www.columbia.edu/~cjd11/charles_dimaggio/DIRE/styled-4/styled-11/code-8/

# See notes on priors in Chapter 5 of Gelman et al. BDA, and 
# http://andrewgelman.com/2009/10/21/some_practical/
setwd("~/Github/AI_surveillance")
sink("base_sampling_events.txt")
cat("model {
    #likelihood (i = month, j = year, k = huc, l = species, s = sampling events)
    for (l in 1:nspecies){
        for (s in 1:nsamplingevents) {
            p[s, l] <- (Se[l] * lambda[s, l]) + ((1-Sp[l]) * (1-lambda[s, l]))
            y[s, l] ~ dbin(p[s, l], n[s, l])
            lambda[s, l] ~ dbeta(alpha[month[s], year[s], huc[s], l], beta[month[s], year[s], huc[s], l])
        }
    }
    
    # Hierarchial step for lambda 
    for (l in 1:nspecies) {
        for (i in 1:nmonths) {
            for (j in l:nyears) {
                for (k in 1:nhucs) {
                    pi[i, j, k, l] <- alpha[i, j, k, l] / (alpha[i, j, k, l] + beta[i, j, k, l])
                    alpha[i, j, k, l] ~ dunif(1, 5)
                    beta[i, j, k, l] ~ dunif(1, 5)
                    # sdpi[i, j, k, l] <- 1/sqrt(alpha[i, j, k, l] + beta[i, j, k, l])
                    # pi[i, j, k, l] ~ dunif(0, 1)
                    # sdpi[i, j, k, l] ~ dt(0, 1, 1)T(0, )  # half cauchy
                }
            }
        }
    }
    
    #priors for Se/Sp
    for (l in 1:nspecies) {
        Se[l] ~ dbeta(20.833, 4.148)
        Sp[l] ~ dbeta(8.403, 1.001)
    }
    }", fill = TRUE)
sink()

# icar model - with sampling events

setwd("~/Github/AI_surveillance")
sink("icar_sampling_events.txt")
cat("model {
    #likelihood (i = month, j = year, k = huc, l = species, s = sampling events)
    for (l in 1:nspecies){
      for (s in 1:nsamplingevents) {
        p[s, l] <- (Se[l] * lambda[s, l]) + ((1-Sp[l]) * (1-lambda[s, l]))
        y[s, l] ~ dbin(p[s, l], n[s, l])
        lambda[s, l] ~ dbeta(alpha[month[s], year[s], huc[s], l], beta[month[s], year[s], huc[s], l])
      }
    }
    
    # Hierarchial step for lambda
    for (l in 1:nspecies) {
      for (i in 1:nmonths) {
        for (j in l:nyears) {
          for (k in 1:nhucs) {
            pi[i, j, k, l] <- alpha[i, j, k, l] / (alpha[i, j, k, l] + beta[i, j, k, l])
            alpha[i, j, k, l] ~ dunif(1, 5)
            beta[i, j, k, l] ~ dunif(1, 5)
            # sdpi[i, j, k, l] <- 1/sqrt(alpha[i, j, k, l] + beta[i, j, k, l])
    
            logit(pi[i, j, k, l]) <- lmupi[i, j, l] + epsilon[k]
    
            # hierarchial priors
            # sdpi[i, j, k, l] ~ dt(0, 1, 1)T(0, )  # half cauchy
          }
        }
      }
    }
    
    
    # CAR Priors - Assumes spatial pattern is the same in all species...
    for (k in 1:nhucs){
      muepsilon[k] <- 0
    }
    Omega <- (tau * I) * (I - W)
    epsilon[1:nhucs] ~ dmnorm(muepsilon[1:nhucs], Omega[1:nhucs, 1:nhucs])
    tau ~ dgamma(0.1, 0.1)
    
    for (l in 1:nspecies) {
      for (i in 1:nmonths){
        for (j in 1:nyears){
          lmupi[i, j, l] ~ dt(0, 0.4, 1)
        }
      }
    }
    # priors for Se/Sp
    for (l in 1:nspecies) {
      Se[l] ~ dbeta(20.833, 4.148)
      Sp[l] ~ dbeta(8.403, 1.001)
    }
    }", fill = TRUE)
sink()





################################################################
################################################################
# Models used for Thesis - not accounting for sampling events
################################################################
################################################################
# Base model - independence in space and time
sink("constant_sensitivity.txt")
cat("model{
    #likelihood (i=month, t=year, k=site)
    for(i in 1:nmonths){
        for(t in 1:nyears){
            for(k in 1:nsites){
                g[i,t,k] <- (Se*pi[i,t,k])+(1-Sp)*(1-pi[i,t,k])
                y[i,t,k] ~ dbin(g[i,t,k], n[i,t,k])
            }
        }
    }
    #priors
    Se ~ dbemonth[s], year[s], huc[s],ta(20.833,4.148)
    Sp ~ dbeta(8.403,1.001)
    for(i in 1:nmonths) {
        for(t in 1:nyears) {
            for(k in 1:nsites) {
                pi[i,t,k] ~ dbeta(1,1)
            }
        }  
    }
}", fill=TRUE)
sink()


# CAR Model - autocorrelation in space (not used)
sink("car_constant_sensitivity.txt")
cat("model{
    # Priors
    Se~dbeta(20.833,4.148)
    Sp~dbeta(8.403,1.001)
    
    # CAR Priors
    for (k in 1:nsites){
        mualpha[k] <- 0
    }
    Omega <- (tau * I) * (I - phi * W)
    alpha[1:nsites] ~ dmnorm(mualpha[1:nsites], Omega[1:nsites, 1:nsites])
    tau ~ dgamma(0.1, 0.1) # check me
    phi ~ dnorm(0.5, 4)
    #phi ~ dunif(-0.999, 0.999)
    
    for (i in 1:nmonths){
        for (j in 1:nyears){
            lmupi[i, j] ~ dt(0, 0.4, 1)
        }
    }
    
    #likelihood (i = huc, j = month, k = year)
    for(i in 1:nmonths){
        for(j in 1:nyears){
            for(k in 1:nsites){
                g[i,j,k] <- (Se*pi[i,j,k])+(1-Sp)*(1-pi[i,j,k])
                y[i,j,k] ~ dbin(g[i,j,k], n[i,j,k])
                logit(pi[i,j,k]) <- lmupi[i, j] + alpha[k]
            }
        }
    } 
    }", fill=TRUE)
sink()


# icar Model
sink("icar_constant_sensitivity.txt")
cat("model{
    # Priors
    Se~dbeta(20.833,4.148)
    Sp~dbeta(8.403,1.001)
 
    # CAR Priors
    for (k in 1:nsites){
        mualpha[k] <- 0
    }
    Omega <- (tau * I) * (I - W)
    alpha[1:nsites] ~ dmnorm(mualpha[1:nsites], Omega[1:nsites, 1:nsites])
    tau ~ dgamma(0.1, 0.1)

    for (i in 1:nmonths){
        for (j in 1:nyears){
            lmupi[i, j] ~ dt(0, 0.4, 1)
        }
    }
    
    for(i in 1:nmonths){
        for(j in 1:nyears){
            for(k in 1:nsites){
                g[i,j,k] <- (Se*pi[i,j,k])+(1-Sp)*(1-pi[i,j,k])
                y[i,j,k] ~ dbin(g[i,j,k], n[i,j,k])
                logit(pi[i,j,k]) <- lmupi[i, j] + alpha[k]
            }
        }
    }
}", fill=TRUE)
sink()

# AR1 - autocorrelation in time (months)
sink("ar1_constant_sensitivity.txt")
cat("model{
    # Priors
    Se~dbeta(20.833,4.148)
    Sp~dbeta(8.403,1.001)
    
    # AR1 Priors
    for (i in 1:nmonths){
        for (j in 1:nyears){
            lmupi[i, j] <- log(pi0[i, j]/(1 - pi0[i, j]))
            pi0[i, j] ~ dunif(0, 1)
        }
    }
    month_ar1 ~ dunif(-0.01, 0.01)   
    
    #likelihood (i = huc, j = month, k = year)
    for(k in 1:nsites){
        for (j in 1:nyears) {
            g[1,j,k] <- (Se*pi[1,j,k])+(1-Sp)*(1-pi[1,j,k])
            y[1,j,k] ~ dbin(g[1,j,k], n[1,j,k])
            logit(pi[1,j,k]) <- lmupi[1, j, k]
        }
    }
    for(i in 2:nmonths){
        for(j in 2:nyears){
            for(k in 1:nsites){
                g[i,j,k] <- (Se*pi[i,j,k])+(1-Sp)*(1-pi[i,j,k])
                y[i,j,k] ~ dbin(g[i,j,k], n[i,j,k])
                logit(pi[i,j,k]) <- lmupi[i, j, k] + month_ar1*lmupi[i-1, j, k]
            }
        }
    } 
    }", fill=TRUE)
sink()
