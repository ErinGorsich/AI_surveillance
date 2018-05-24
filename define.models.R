# Define a sequence of JAGS models for model selection

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
    Se ~ dbeta(20.833,4.148)
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



# Model 2 - CAR - autocorrelation in space
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


# sink("car_constant_sensitivity.txt")
# cat("model{
#     # Priors
#     Se~dbeta(20.833,4.148)
#     Sp~dbeta(8.403,1.001)
#     
#     # CAR Priors
#     for (k in 1:nsites){
#     mualpha[k] <- 0
#     }
#     Omega <- tau * (D - phi * W)
#     alpha[1:nsites] ~ dmnorm(mualpha[1:nsites], Omega[1:nsites, 1:nsites])
#     tau ~ dgamma(0.1, 0.1) # check me
#     phi ~ dunif(-0.999, 0.999)
#     
#     for (i in 1:nmonths){
#     for (j in 1:nyears){
#     lmupi[i, j] <- log(pi0[i, j]/(1 - pi0[i, j]))
#     pi0[i, j] ~ dunif(0, 1)
#     }
#     }
#     
#     #likelihood (i = huc, j = month, k = year)
#     for(i in 1:nmonths){
#     for(j in 1:nyears){
#     for(k in 1:nsites){
#     g[i,j,k] <- (Se*pi[i,j,k])+(1-Sp)*(1-pi[i,j,k])
#     y[i,j,k] ~ dbin(g[i,j,k], n[i,j,k])
#     logit(pi[i,j,k]) <- lmupi[i, j] + alpha[k]
#     }
#     }
#     } 
#     }", fill=TRUE)
# sink()

# Model 2.5 icar
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









# Model 3 - CAR - autocorrelation in space and time (months)
sink("ar1_car_constant_sensitivity.txt")
cat("model{
    # Priors
    Se~dbeta(20.833,4.148)
    Sp~dbeta(8.403,1.001)
    
    # CAR Priors
    for (k in 1:nsites){
        mualpha[k] <- 0
    }
    Omega <- tau * (D - phi * W)
    alpha[1:nsites] ~ dmnorm(mualpha[1:nsites], Omega[1:nsites, 1:nsites])
    tau ~ dgamma(0.1, 0.1) # check me
    phi ~ dunif(-0.999, 0.999)
    
    for (i in 1:nmonths){
        for (j in 1:nyears){
            lmupi[i, j] ~ dt(0, 0.4, 1)
        }
    }
    
    # AR1 Priors
    month_ar1 ~ dunif(-0.01, 0.01)   

    #likelihood (i = huc, j = month, k = year)
    for(k in 1:nsites){
        for (j in 1:nyears){}
            g[1,j,k] <- (Se*pi[1,j,k])+(1-Sp)*(1-pi[1,j,k])
            y[1,j,k] ~ dbin(g[1,j,k], n[1,j,k])
            logit(pi[1,j,k]) <- lmupi[1, j] + alpha[k]
        }
    }
    for(i in 2:nmonths){
        for(j in 2:nyears){
            for(k in 1:nsites){
                g[i,j,k] <- (Se*pi[i,j,k])+(1-Sp)*(1-pi[i,j,k])
                y[i,j,k] ~ dbin(g[i,j,k], n[i,j,k])
                logit(pi[i,j,k]) <- lmupi[i, j] + alpha[k] + month_ar1*lmupi[i-1, j, k]
            }
        }
    } 
    }", fill=TRUE)
sink()


# Model 3 - AR1 - autocorrelation in time (months)
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
