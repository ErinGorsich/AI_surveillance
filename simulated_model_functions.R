##########################################################
##########################################################
#functions for simulating data
##########################################################
##########################################################

base.model <- function (n, lambda) {
  app.prev <- array(dim=c(nmonths, nyears, nhucs), NA)
  for (i in 1:nmonths) {
    for (j in 1:nyears) {
      for (k in 1:nhucs) {
        Se <- rbeta(1, 20.833, 4.148)
        Sp <- rbeta(1, 8.403, 1.001)
        # y[i,j,k] ~ rbinom(n[i,j,k], lambda[i,j,k])
        app.prev[i,j,k] = Se*lambda[i,j,k] + Sp*(1-lambda[i,j,k])
      }
    }
  }
  
  return(app.prev)
}

base.model <- function (n, lambda) {
  app.prev <- array(dim=c(nmonths, nyears, nhucs), NA)
  Se <- array(dim=c(nmonths, nyears, nhucs), NA)
  Sp <- array(dim=c(nmonths, nyears, nhucs), NA)
  for (i in 1:nmonths){
    for (j in 1:nyears){
      for (k in 1:nhucs){
        Se[i,j,k] <- rbeta(1, 20.833, 4.148)
        Sp[i,j,k] <- rbeta(1, 8.403, 1.001)
        omega <- rbinom(1, n[i,j,k], lambda[i,j,k])
        sigma <- omega/n[i,j,k]
        app.prev[i,j,k] <- Se[i,j,k]*sigma + Sp[i,j,k]*(1-sigma)
      }
    }
  }
  list <- list(app.prev = app.prev, Se=Se, Sp=Sp)
  return(list)
}
