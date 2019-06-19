#########################################################################################
#########################################################################################
#plot base model fit data from simulated data
#########################################################################################
#########################################################################################

base.fit <- readRDS("/home/webblab/Documents/HP/base_model_fit_sim_6.rds")
lambda <- readRDS('/home/webblab/Documents/HP/sim_data_6_lamda.rds')

nmonths <- 12
nyears <- 5
nhuc <- 191
nparam <- 2+(12*5*191)
nchains <- 3

#calculate a mean for each month, year, huc combination
y <- matrix(ncol=nparam, nrow=nchains, NA)
z <- matrix(ncol=nparam, nrow=1, NA)
for (l in 1:nparam) {
  for (m in 1:nchains){
    y[m,l] <- quantile(base.fit[[m]][,l], probs=0.5)
  }
}
for (l in 1:nparam){
  for (m in 1:nchains){
    z[1,l] <- quantile(y[,l], probs=0.5)
  }
}
mean <- array(dim=c(nmonths, nyears, nhuc), NA)
for (i in 1:nmonths) {
  for (j in 1:nyears) {
    for (k in 1:nhucs) {
      l <- (k-1)*60 + (j-1)*12 + (i+2)
      mean[i,j,k] <- quantile(z[,l], probs=0.5)
    }
  }
}
mean.huc <- matrix(nrow=nhuc, ncol=1, NA)
for (k in 1:nhuc){
  mean.huc[k,1] <- quantile(mean[,,k], probs=0.5)
}
mean.all <- quantile(mean.huc[,1], probs=0.5)
stand.dev <- sd(mean.huc[,1])
outlier <- matrix(nrow=nhuc, ncol=1, NA)
for (k in 1:nhuc){
  if (mean.huc[k,1]>mean.all+stand.dev) {
    outlier[k,1] <- TRUE
  } else if (mean.huc[k,1]<mean.all-stand.dev) {
    outlier[k,1] <- TRUE
  } else {
    outlier[k,1] <- FALSE
  }
}
outlier <- as.data.frame(outlier)

difference <- array(dim=c(nmonths, nyears, nhuc), NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      difference[i,j,k] <- lambda[i,j,k]-mean[i,j,k]
    }
  }
}

plot <- data.frame(month = rep(seq(1,12), times = nyears*nhuc),
         year = rep(seq(1,5), each = nmonths, times=nhuc),
         huc = rep(seq(1,191), each = nmonths*nyears),
         lambda.true = NA, lambbda.est = NA, difference = NA, mean = NA, outlier = NA)
for (i in 1:nmonths) {
  for (j in 1:nyears) {
    for (k in 1:nhuc) {
      plot$lambda.true[plot$month==i & plot$year==j & plot$huc==k] <- lambda[i,j,k]
      plot$lambbda.est[plot$month==i & plot$year==j & plot$huc==k] <- mean[i,j,k]
      plot$difference[plot$month==i & plot$year==j & plot$huc==k] <- difference[i,j,k]
      plot$mean[plot$huc==k] <- mean.huc[k,1]
      plot$outlier[plot$huc==k] <- outlier$V1[k]
    }
  }
}
##################################################################################################
#make plots
##################################################################################################

plot.1 <- plot[plot$outlier==TRUE, ]

outlier.hucs <- ggplot(plot, aes(x=huc, y=da))
