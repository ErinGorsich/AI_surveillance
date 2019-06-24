#########################################################################################
#########################################################################################
#plot base model fit data from simulated data
#########################################################################################
#########################################################################################
library(ggplot2)
library(dplyr)

base.fit <- readRDS("/home/webblab/Documents/HP/simulated_results/base_model_fit_sim_6_short20.rds")
lambda <- readRDS('/home/webblab/Documents/HP/sim_data_6_lamda_short20.rds')
# base.fit <- readRDS('~/Github/base_model_fit_sim_6.rds')
# lambda <- readRDS('~/Github/sim_data_6_lamda.rds')

nmonths <- dim(lambda)[1]
nyears <- dim(lambda)[2]
nhucs <- dim(lambda)[3]
nparam <- 2+(nmonths*nhucs*nyears)
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
mean <- array(dim=c(nmonths, nyears, nhucs), NA)
for (i in 1:nmonths) {
  for (j in 1:nyears) {
    for (k in 1:nhucs) {
      l <- (k-1)*60 + (j-1)*12 + (i+2)
      mean[i,j,k] <- quantile(z[,l], probs=0.5)
    }
  }
}
#calculate a huc level mean of the lambda estimate
mean.huc <- matrix(nrow=nhucs, ncol=1, NA)
stand.dev.huc <- matrix(nrow=nhucs, ncol=1, NA)
for (k in 1:nhucs){
  mean.huc[k,1] <- quantile(mean[,,k], probs=0.5)
  stand.dev.huc[k,1] <- sd(mean[,,k])
}
mean.all <- quantile(mean.huc[,1], probs=0.5)
#determine if a huc is an outlier
#outlier = outside +/- overall standard deviation
stand.dev.all <- sd(mean.huc[,1])
outlier <- matrix(nrow=nhucs, ncol=1, NA)
for (k in 1:nhucs){
  if (mean.huc[k,1]>mean.all+stand.dev.all) {
    outlier[k,1] <- TRUE
  } else if (mean.huc[k,1]<mean.all-stand.dev.all) {
    outlier[k,1] <- TRUE
  } else {
    outlier[k,1] <- FALSE
  }
}
outlier <- as.data.frame(outlier)
#calculate a mean "true" lambda for each huc
mean.huc.true <- matrix(nrow=nhucs, ncol=1, NA)
sd.huc.true <- matrix(nrow=nhucs, ncol=1, NA)
max.huc.true <- matrix(nrow=nhucs, ncol=1, NA)
for (k in 1:nhucs){
  mean.huc.true[k,1] <- quantile(lambda[,,k], probs=0.5)
  sd.huc.true[k,1] <- sd(lambda[,,k])
  max.huc.true[k,1] <- max(lambda[,,k])
}
#find the difference between the "true" lambda mean and 
# estimated lambda mean for each huc
# difference <- array(dim=c(nmonths, nyears, nhucs), NA)
# for (i in 1:nmonths){
#   for (j in 1:nyears){
#     for (k in 1:nhucs){
#       difference[i,j,k] <- lambda[i,j,k]-mean[i,j,k]
#     }
#   }
# }
mean.difference <- matrix(ncol=1, nrow=nhucs, NA)
for (k in 1:nhucs){
  mean.difference[k,1] <- mean.huc[k,1] - mean.huc.true[k,1]
}

plot <- data.frame(huc = as.factor(seq(1, nhucs)), lambda.true = mean.huc.true, 
                   lambda.true.max = max.huc.true, lambda.est = mean.huc, 
                   sd.true = sd.huc.true, sd.est = stand.dev.huc, 
                   difference = mean.difference, outlier = outlier$V1)

##################################################################################################
#make plots
##################################################################################################

plot.1 <- plot[plot$outlier==TRUE, ]
for (i in 1:length(plot.1$lambda.est)){
  plot.1$ci.95.min[i] <- 2*plot.1$sd.est[i]+plot.1$lambda.est[i]
  plot.1$ci.95.max[i] <- 2*plot.1$sd.est[i]-plot.1$lambda.est[i]
}

outlier.hucs <- ggplot(plot.1, aes(x=huc)) +
  geom_point(aes(y = lambda.est), color = "darkblue") +
  geom_errorbar(aes(ymin=ci.95.min, ymax=ci.95.max), color="darkblue") +
  geom_point(aes(y= lambda.true), color = "darkred") +
  geom_errorbar(aes(ymin=0, ymax=lambda.true.max), color="darkred")+
  geom_hline(aes(yintercept=mean.all)) +
  theme_minimal() +
  ggtitle("Comparison of Simulated and Estimated True Prevalence for Outlier Watersheds")+
  xlab("Watershed") +
  ylab("Prevalence")

plot <- plot %>% arrange(difference)

plot.2 <- plot[1:20, ]  

smallest.diff <- ggplot(plot.2, aes(x=huc)) +
  geom_point(aes(y = difference), color="darkgreen") +
  geom_errorbar(aes(ymin=lambda.true, ymax=lambda.est), color="darkgreen")+
  theme_minimal()+
  ggtitle("20 Watersheds with the Smallest Difference between True and Estimated Prevalence")+
  xlab("Watershed") +
  ylab("Difference")

plot.3 <- plot[(nhucs-20):nhucs, ]

largest.diff <- ggplot(plot.3, aes(x=huc)) +
  geom_point(aes(y=difference), color="darkorange")+
  geom_errorbar(aes(ymin=lambda.true, ymax=lambda.est), color="darkorange")+
  theme_minimal()+
  ggtitle("20 Watersheds with the Largest Difference between True and Estimated Prevalence")+
  xlab("Watershed")+
  ylab("Difference")

plot.4 <- data.frame(month=as.factor(rep(seq(1,12), times=5)), year=as.factor(rep(seq(1,5), each=12)),
                     fitted.prev = NA, fitted.prev.min=NA, fitted.prev.max=NA,
                     lambda = NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    hold <- mean[i,j, ]
    plot.4$fitted.prev[plot.4$month==i & plot.4$year==j] <- quantile(hold, probs=0.5)
    plot.4$fitted.prev.min[plot.4$month==i & plot.4$year==j] <- min(hold)
    plot.4$fitted.prev.max[plot.4$month==i & plot.4$year==j] <- max(hold)
    hold.1 <- lambda[i,j,]
    plot.4$lambda[plot.4$month==i & plot.4$year==j] <- quantile(hold.1, probs=0.5)
  }
}

fitted.month.mean <- ggplot(plot.4, aes(x=month, y=fitted.prev))+
  geom_point(color="darkorchid4")+
  geom_errorbar(aes(ymin=fitted.prev.min, ymax=fitted.prev.max), color="darkorchid4")+
  facet_grid(.~ year)+
  theme_minimal()+
  ggtitle('Fitted Data Mean True Prevalence by Month')+
  ylab('Mean Fitted True Prevalence')

compare.means <- ggplot(plot.4, aes(x=lambda, y=fitted.prev))+
  geom_point(color="dodgerblue4")+
  theme_minimal()+
  ylim(0, 1)+
  ggtitle('Given vs. Fitted True Prevalence (Monthly Mean)')

plot.5 <- data.frame(month=as.factor(rep(seq(1,12), times=nyears*nhucs)), year=as.factor(rep(seq(1,5), each=nmonths, times=nhucs)),
                     huc=rep(seq(1,nhucs), each=nyears*nmonths), fitted.prev=NA, lambda=NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      plot.5$fitted.prev[plot.5$month==i & plot.5$year==j & plot.5$huc==k] <- mean[i,j,k]
      plot.5$lambda[plot.5$month==i & plot.5$year==j & plot.5$huc==k] <- lambda[i,j,k]
    }
  }
}

compare <- ggplot(plot.5, aes(x=lambda, y=fitted.prev))+
  geom_point(color="dodgerblue4")+
  theme_minimal()+
  ggtitle('Given vs. Fitted True Prevalence')+
  ylim(0,1)+
  ylab('Estimated True Prevalence')+
  xlab('Given True Prevalence')

