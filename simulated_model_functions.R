##########################################################
##########################################################
#functions for simulating data
##########################################################
##########################################################

base.model <- function (n, lambda, Sp, Se) {
  app.prev <- array(dim=c(nmonths, nyears, nhucs), NA)
  for (i in 1:nmonths) {
    for (j in 1:nyears) {
      for (k in 1:nhucs) {
        y[i,j,k] ~ rbinom(n[i,j,k], lambda[i,j,k])
        app.prev[i,j,k] = Se*lambda[i,j,k] + Sp*(1-lambda[i,j,k])
      }
    }
  }
  
  return(app.prev)
}
