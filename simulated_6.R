############################################################################################
############################################################################################
#Simulated Data Attempt 6
#create a data set of apparent prevalence based on true prevalence given to the model
############################################################################################
############################################################################################

data <- readRDS("/home/webblab/Documents/HP/locationsamplingevent_n_y_speciesgroup.rds")

data <- data[data$species.group %in% c(1,2,4), ]
for (i in 1:length(data$species.group)) {
  if (data$species.group[i] == 4) {
    data$species.group[i] <- 3
  } else {data$species.group[i] <- data$species.group[i]}
}
data$YEAR <- NA
for (i in 1:length(data$year)) {
  if (data$year[i] == "2007") {
    data$YEAR[i] <- 1
  } else if (data$year[i] == "2008") {
    data$YEAR[i] <- 2
  } else if (data$year[i] == "2009") {
    data$YEAR[i] <- 3
  } else if (data$year[i] == "2010") {
    data$YEAR[i] <- 4
  } else {
    data$YEAR[i] <- 5
  }
}
data$month <- as.numeric(data$month)

huc <- data.frame(watershed = unique(data$watershed), huc = seq(1,length(unique(data$watershed))))

for (i in 1:length(data$watershed)) {
  data$HUC[i] <- huc$huc[huc$watershed == data$watershed[i]]
}

nmonths <- length(unique(data$month))
nyears <- length(unique(data$YEAR))
nhucs <- length(unique(data$HUC))

tally <- data.frame(species.group = rep(seq(1,3), each = nyears*nmonths*nhucs),
                    month = rep(seq(1,12), times = nyears*nhucs*3),
                    year = rep(seq(1,5), each = 12, times = nhucs*3),
                    huc = rep(seq(1, nhucs), each = nyears*nmonths, rep = 3), 
                    n = NA, y = NA, app.prev = NA)
for (h in 1:3) {
  for (i in 1:nmonths) {
    for (j in 1:nyears) {
      for (k in 1:nhucs) {
        hold <- data$n[data$species.group == h & data$month == i & data$YEAR == j & data$HUC == k]
        tally$n[tally$species.group == h & tally$month == i & tally$year == j & tally$huc == k] <-
          sum(hold)
        hold <- data$y[data$species.group == h & data$month == i & data$YEAR == j & data$HUC ==k]
        tally$y[tally$species.group == h & tally$month == i & tally$year == j & tally$huc ==k] <- 
          sum(hold)
      }
    }
  }
}
for (i in 1:length(tally$species.group)) {
  tally$app.prev[i] <- tally$y[i]/tally$n[i]
  if (is.na(tally$app.prev[i] == TRUE)) {
    tally$app.prev[i] <- 0
  }
}

#remove hucs with no dabbling duck data and renumber hucs
n.zero <- matrix(ncol=1, nrow=nhucs, NA)
for (k in 1:nhucs) {
  hold <- tally[tally$huc == k & tally$species.group==1, ]
  n.zero[k,1] <- all(hold$n == 0)
}
n.zero <- as.data.frame(n.zero)
n.zero$huc <- seq(1:194)
zero.hucs <- n.zero$huc[n.zero$V1 == TRUE]
tally <- tally[!(tally$huc %in% zero.hucs),]
renumber <- data.frame(huc = unique(tally$huc), renumber = seq(1:191))
tally$huc.new <- NA
for (i in 1:length(tally$huc)) {
  tally$huc.new[i] <- renumber$renumber[renumber$huc==tally$huc[i]]
}
tally$huc <- NULL
tally$huc <- tally$huc.new

tally$true.prev <- NA
spec = 0.999
sens = 0.83
for (i in 1:length(tally$species.group)){
  tally$true.prev[i] <- (tally$app.prev[i]+(spec-1))/(spec+(sens-1))
  if (tally$true.prev[i] > 1) {
    tally$true.prev[i] <- 1
  } else if (tally$true.prev[i] < 0) {
    tally$true.prev[i] <- abs(tally$true.prev[i])
  } else {
    tally$true.prev[i] <- tally$true.prev[i]
  }
}
for (i in 1:length(tally$true.prev)) {
  if (tally$app.prev[i]<1 & tally$true.prev[i]==1) {
    tally$true.prev[i] <- tally$true.prev[i] - (1-tally$app.prev[i])
  } else {
    tally$true.prev[i] <- tally$true.prev[i]
  }
}

############################################################################################################
#simulate data using model function
#species.group==1
############################################################################################################
nmonths <- length(unique(tally$month))
nyears <- length(unique(tally$year))
nhucs <- length(unique(tally$huc))

n <- array(dim=c(nmonths, nyears, nhucs), NA)
lambda <- array(dim=c(nmonths, nyears, nhucs), NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      n[i,j,k] <- tally$n[tally$species.group==1 & tally$month==i & tally$year==j & tally$huc==k]
      lambda[i,j,k] <- tally$true.prev[tally$species.group==1 & tally$month==i & tally$year==j & tally$huc==k]
    }
  }
}
# Se <- rbeta(nmonths*nyears*nhucs, 20.833, 4.148)
# Sp <- rbeta(nmonths*nyears*nhucs, 8.403, 1.001)
# Se <- mean(Se)
# Sp <- mean(Sp)
# #Sp=0.8941578, Se=0.8340683

simulated <- base.model(n=n, lambda=lambda)
app.prev <- simulated$app.prev
Se <- simulated$Se
Sp <- simulated$Sp

# saveRDS(app.prev, "/home/webblab/Documents/HP/sim_data_attempt6.rds")
# saveRDS(n, "/home/webblab/Documents/HP/sim_data_6_n.rds")
# saveRDS(lambda, "/home/webblab/Documents/HP/sim_data_6_lamda.rds")

#######################################################################################################
#create a y matrix for model fittiing based on the generated apparent prevalences
#######################################################################################################
y <- array(dim=c(nmonths, nyears, nhucs), NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      positive <- app.prev[i,j,k]*n[i,j,k]
      y[i,j,k] <- round(positive, digits = 0)
      if (is.na(y[i,j,k]) == TRUE) {
        y[i,j,k] <- 0
      } else {
        y[i,j,k] <- y[i,j,k]
      }
    }
  }
}

# saveRDS(y, "/home/webblab/Documents/HP/sim_data_6_y.rds")

#####################################################################################################
#find mean Sp and Se
#####################################################################################################
means.sp <- matrix(ncol = 1, nrow = nhucs, NA)
means.se <- matrix(ncol = 1, nrow = nhucs, NA)
for (k in 1:nhucs) {
    means.sp[k,1] <- quantile(Sp[ , ,k], probs=0.5)
    means.se[k,1] <- quantile(Se[ , ,k], probs=0.5)
}
sp.mean <- quantile(means.sp, probs=0.5)
se.mean <- quantile(means.se, probs=0.5)
    
###################################################################################################
#trying to figure out what was wrong
###################################################################################################
# for (i in 1:nmonths) {
#   for (j in 1:nyears) {
#     for (k in 1:nhucs) {
#       hold <- tally[tally$huc == k,]
#       all(hold$app.prev == 0 & hold$n == 0)
#       all(app.prev[,,k] == 1)
#     }
#   }
# }
# 
# answer <- matrix(nrow=194, ncol=1, NA)
# for (k in 1:nhucs){
#   hold <- app.prev[,,k]
#   for (i in 1:nmonths){
#     for (j in 1:nyears){
#       if (is.na(hold[i,j])==TRUE) {
#         hold[i,j] <- 2
#       } else {
#         hold[i,j] <- hold[i,j]
#       }
#     }
#   }
#   answer[k,] <- all(hold == 2)
# }
# 
# answer <- as.data.frame(answer)
# answer$huc <- seq(1,194)
# see <- answer[answer$V1==TRUE, ]

#############################################################################################
#plot simulated data
#############################################################################################

plot <- data.frame(month = as.factor(rep(seq(1,12), times = nyears*nhucs)), 
                   year = as.factor(rep(seq(1,5), each = nmonths, times=nhucs)),
                   huc = rep(seq(1,nhucs), each=nmonths*nyears), n=NA, y=NA, lambda=NA, p=NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      plot$n[plot$month==i & plot$year==j & plot$huc==k] <- n[i,j,k]
      plot$y[plot$month==i & plot$year==j & plot$huc==k] <- y[i,j,k]
      plot$lambda[plot$month==i & plot$year==j & plot$huc==k] <- lambda[i,j,k]
      plot$p[plot$month==i & plot$year==j & plot$huc==k] <- app.prev[i,j,k]
    }   
  }
}

plot.2 <- data.frame(month = plot$month, year = plot$year, n=NA)

for (i in 1:nmonths){
  for (j in 1:nyears){
    hold <- plot$n[plot$month==i & plot$year==j]
    plot.2$n[plot.2$month==i & plot.2$year==j] <- sum(hold)
  }
}
palette <- brewer.pal(12, "Set3")

total.samples.plot <- ggplot(plot.2, aes(x=year, y=n, fill=month))+
  geom_bar(stat="identity", position="dodge") +
  theme_minimal()+
  scale_fill_manual(values=palette) +
  # ylim(0,10000)+
  ggtitle("Total Number of Dabbling Duck Samples by Month")

plot.3 <- data.frame(month=plot$month, year=plot$year, n=plot.2$n, 
                     lambda.mean=NA, lambda.min=NA, lambda.max=NA,
                     p.mean=NA, p.min=NA, p.max=NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    hold <- plot$lambda[plot$month==i & plot$year==j]
    for (l in 1:length(hold)) {
      if (is.na(hold[l]==TRUE)) {
        hold[l] <- 0
      } else {
        hold[l] <- hold[l]
      }
    }
    plot.3$lambda.mean[plot.3$month==i & plot.3$year==j] <- quantile(hold, prob=0.5)
    plot.3$lambda.max[plot.3$month==i & plot.3$year==j] <- max(hold)
    plot.3$lambda.min[plot.3$month==i & plot.3$year==j] <- min(hold)
    hold.1 <- plot$p[plot$month==i & plot$year==j]
    for (l in 1:length(hold.1)) {
      if (is.na(hold.1[l]==TRUE)) {
        hold.1[l] <- 0
      } else {
        hold.1[l] <- hold[l]
      }
    }
    plot.3$p.mean[plot.3$month==i & plot.3$year==j] <- quantile(hold.1, prob=0.5)
    plot.3$p.min[plot.3$month==i & plot.3$year==j] <- min(hold.1)
    plot.3$p.max[plot.3$month==i & plot.3$year==j] <- max(hold.1)
  }
}

lambda.plot <- ggplot(plot.3, aes(x=month, y=lambda.mean)) +
  geom_point(color="purple4")+
  geom_errorbar(aes(ymin=lambda.min, ymax=lambda.max), color="purple4")+
  facet_grid(.~ year) +
  theme_minimal() +
  ylim(0,1)+
  ggtitle("Dabbling Duck Given True Prevalences by Month")+
  ylab("Mean Given True Prevalence")

p.plot <- ggplot(plot.3, aes(x=month, y=p.mean)) +
  geom_point(color="violetred4") +
  geom_errorbar(aes(ymin=p.min, ymax=p.max), color="violetred4")+
  facet_grid(.~ year)+
  theme_minimal()+
  ggtitle("Dabbling Duck Mean Simulated Apparent Prevalences by Month") +
  ylab("Mean Simulated Apparent Prevalence")

#########################################################################################################
#remove hucs that do not have sufficient data
#########################################################################################################
#what is sufficient data?
  #>5 combinations with n>0
# how to find these hucs?
  # start with n array
  # 0 -> NA
# dummy array is.na TRUE/FALSE
  # TRUE = no data
  # FALSE = data
#remove hucs
  #see above for previous fix
  #remove pieces of dummy array based with >55 TRUES
    #make n.na a data.frame
    #tally number of TRUES for each huc
    #remove hucs with >55 TRUES
  #remove pieces of n based on dummy array
  #renumber hucs
#remake n, lambda, app.prev, and y matrices

n.new <- array(dim=c(nmonths, nyears, nhucs), NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      if (n[i,j,k]==0){
        n.new[i,j,k] <- NA
      } else {
        n.new[i,j,k] <- n[i,j,k]
      }
    }
  }
}

n.na <- array(dim=c(nmonths, nyears, nhucs), NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      n.na[i,j,k] <- is.na(n.new[i,j,k])
    }
  }
}

na.data <- data.frame(month=rep(seq(1,12), times=nyears*nhucs), year=rep(seq(1,5), each=nmonths, times=nhucs),
                      huc=rep(seq(1,nhucs), each=nmonths*nyears), n.na = NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      na.data$n.na[na.data$month==i & na.data$year==j & na.data$huc==k] <- n.na[i,j,k]
    }
  }
}
na.true <- na.data[na.data$n.na==TRUE, ]
na.count <- count(na.data[na.data$n.na==TRUE,], huc, sort = TRUE)
na.sufficient <- na.count[!(na.count$n>=50), ]
na.sufficient$huc.2 <- NA

renumber <- data.frame(huc.new = seq(1,length(na.sufficient$huc)), huc = na.sufficient$huc)

for (l in 1:length(na.sufficient)){
  m <- na.sufficient$huc[l]
  na.sufficient$huc.2[l] <- renumber$huc.new[renumber$huc==m] 
}

tally.mall <- tally[tally$species.group==1, ]
tally.mall$huc.2 <- NA
for (l in 1:length(tally.mall$huc.2)){
  m <- tally.mall$huc.new[l]
  o <- renumber$huc.new[renumber$huc==m]
  if (length(o)>0){
    o <- o
  } else {
    o <- NA
  }
  tally.mall$huc.2[l] <- o
}
tally.mall.short <- tally.mall[!(is.na(tally.mall$huc.2)), ]

nmonths <- length(unique(tally.mall.short$month))
nyears <- length(unique(tally.mall.short$year))
nhucs <- length(unique(tally.mall.short$huc.2))

n.short <- array(dim=c(nmonths, nyears, nhucs), NA)
lambda.short <- array(dim=c(nmonths, nyears, nhucs), NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      n.short[i,j,k] <- tally.mall.short$n[tally.mall.short$species.group==1 & tally.mall.short$month==i &
                                             tally.mall.short$year==j & tally.mall.short$huc.2==k]
      lambda.short[i,j,k] <- tally.mall.short$true.prev[tally.mall.short$species.group==1 & tally.mall.short$month==i &
                                                          tally.mall.short$year==j & tally.mall.short$huc.2==k]
    }
  }
}

simulated <- base.model(n=n.short, lambda=lambda.short)
app.prev <- simulated$app.prev
Se <- simulated$Se
Sp <- simulated$Sp

saveRDS(app.prev, "/home/webblab/Documents/HP/sim_data_attempt6_short10.rds")
saveRDS(n.short, "/home/webblab/Documents/HP/sim_data_6_n_short10.rds")
saveRDS(lambda.short, "/home/webblab/Documents/HP/sim_data_6_lamda_short10.rds")

y.short <- array(dim=c(nmonths, nyears, nhucs), NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      positive <- app.prev[i,j,k]*n.short[i,j,k]
      y.short[i,j,k] <- round(positive, digits = 0)
      if (is.na(y.short[i,j,k]) == TRUE) {
        y.short[i,j,k] <- 0
      } else {
        y.short[i,j,k] <- y.short[i,j,k]
      }
    }
  }
}

saveRDS(y.short, "/home/webblab/Documents/HP/sim_data_6_y_short10.rds")

#############################################################################################
#plot short simulated data
#############################################################################################

plot.short <- data.frame(month = as.factor(rep(seq(1,12), times = nyears*nhucs)), 
                   year = as.factor(rep(seq(1,5), each = nmonths, times=nhucs)),
                   huc = rep(seq(1,nhucs), each=nmonths*nyears), n=NA, y=NA, p=NA)
for (i in 1:nmonths){
  for (j in 1:nyears){
    for (k in 1:nhucs){
      plot.short$n[plot.short$month==i & plot.short$year==j & plot.short$huc==k] <- n.short[i,j,k]
      plot.short$y[plot.short$month==i & plot.short$year==j & plot.short$huc==k] <- y.short[i,j,k]
      plot.short$p[plot.short$month==i & plot.short$year==j & plot.short$huc==k] <- app.prev[i,j,k]
    }   
  }
}

plot.2.short <- data.frame(month = plot.short$month, year = plot.short$year, n=NA)

for (i in 1:nmonths){
  for (j in 1:nyears){
    hold <- plot.short$n[plot.short$month==i & plot.short$year==j]
    plot.2.short$n[plot.2.short$month==i & plot.2.short$year==j] <- sum(hold)
  }
}
palette <- brewer.pal(12, "Set3")

total.samples.plot <- ggplot(plot.2.short, aes(x=year, y=n, fill=month))+
  geom_bar(stat="identity", position="dodge") +
  theme_minimal()+
  scale_fill_manual(values=palette) +
  ylim(0,10000)+
  ggtitle("Total Number of Dabbling Duck Samples by Month")
