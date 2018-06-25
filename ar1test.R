setwd("~/GitHub/AI_surveillance")

sink("ar1_sampling_events.txt")
cat("model {
    #likelihood (i = month, j = year, k = huc, l = species, s = sampling events)
    for (l in 1:nspecies){
        for (s in 1:nsamplingevents) {
            y[s, l] ~ dbin(p[s, l], n[s, l])
            p[s, l] <- (Se[l] * lambda[s, l]) + ((1-Sp[l]) * (1-lambda[s, l]))
            lambda[s, l] ~ dbeta(alpha[month[s], year[s], huc[s], l], beta[month[s], year[s], huc[s], l])T(0.001,0.999)
        }
    }
    # Hierarchial step for lambda 
    for (l in 1:nspecies) {
        for (i in 1:nmonths) {
            for (j in 1:nyears) {
                for (k in 1:nhucs) {
                    alpha[i, j, k, l] <- pi[i, j, k, l] * scale[i, j, k, l]
                    beta[i, j, k, l] <- (1 - pi[i, j, k, l]) * scale[i, j, k, l]
                    scale[i, j, k, l] ~ dpar(1.5, 1)
                    pi[i, j, k, l] <- exp(logitpi[i, j, k, l])/ (1 + exp(logitpi[i, j, k, l]))
                    logitpi[i, j, k, l] ~ dt(0, 0.4, 1)
                }
            }
        }
    }
    
    # AR1 prior (residuals at i=1 related to i=2 in the same site/year)
    for (l in 1:nspecies) {
        for (j in 1:nyears) {
            for (k in 1:nhucs) {
                epsilon[1, j, k, l] <- pi[1, j, k, l]/(1-pi[1, j, k, l]) - lmupi[j, k, l]
                lmupi[j, k, l] ~ dt(0, 0.4, 1)
                for (i in 2:nmonths){
                    epsilon[i, j, k, l] ~ dnorm(phi * epsilon[i-1, j, k, l], tau.corr)
                }
            }
        }
    }
       
    tau.corr ~ dgamma(0.001, 0.001)  # THINK ABOUT ME!
    phi ~ dunif(-1, 1)
    
    # priors for Se/Sp
    for (l in 1:nspecies) {
        Se[l] ~ dbeta(20.833, 4.148)
        Sp[l] ~ dbeta(8.403, 1.001)
    }
}", fill = TRUE)
sink()
