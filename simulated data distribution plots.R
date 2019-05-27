###################################################################################
###################################################################################
#Check distribution of simulated data across watersheds
###################################################################################
###################################################################################
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

simulated <- readRDS("/home/webblab/Documents/HP/simulated_data.rds")

#break down data by species group
sim.mall <- simulated[simulated$species.group == 1, ]
sim.diving <- simulated[simulated$species.group == 2, ]
sim.geese <- simulated[simulated$species.group == 3, ]

###################################################################################
#plot distribution of Dabbling ducks data by month and year
###################################################################################
total <- matrix(nrow = 5, ncol= 12)
for (i in 1:5) {
  for (j in 1:12) {
    y <- as.numeric(sim.mall %>% tally(sim.mall$total[sim.mall$month == j & sim.mall$year == i]))
    total[i, j] <- y
  }
}

total.mall <- data.frame(year = as.factor(rep(seq(1,5), each = 12)),
                         month = as.factor(rep(seq(1,12), times = 5)))
total.mall$total <- NA
for (i in 1:5) {
  for (j in 1:12) {
    y <- total[i,j]
    total.mall$total[total.mall$year == i & total.mall$month == j] <- y
  }
}

palette <- brewer.pal(12, "Set3")

mall.plot <- ggplot(total.mall, aes(x=year, y = total, fill = month)) +
  geom_bar(stat = "identity", position = "dodge")+
  theme_minimal() +
  scale_fill_manual(values = palette) +
  ylim(0, 10000) +
  ggtitle("Total Number of Dabbling Duck Samples by Month")

jpeg('/home/webblab/Documents/HP/distribution_plots/simulated_dabbling_samplenumber.jpeg')
  mall.plot
dev.off()
  
###################################################################################
#plot distribution of Diving Ducks data by month and year
###################################################################################
total <- matrix(nrow = 5, ncol= 12)
for (i in 1:5) {
  for (j in 1:12) {
    y <- as.numeric(sim.diving %>% tally(sim.diving$total[sim.diving$month == j & sim.diving$year == i]))
    total[i, j] <- y
  }
}

total.diving <- data.frame(year = as.factor(rep(seq(1,5), each = 12)),
                         month = as.factor(rep(seq(1,12), times = 5)))
total.diving$total <- NA
for (i in 1:5) {
  for (j in 1:12) {
    y <- total[i,j]
    total.diving$total[total.diving$year == i & total.diving$month == j] <- y
  }
}

diving.plot <- ggplot(total.diving, aes(x=year, y = total, fill = month)) +
  geom_bar(stat = "identity", position = "dodge")+
  theme_minimal() +
  scale_fill_manual(values = palette) +
  ylim(0, 1000) +
  ggtitle("Total Number of Diving Duck Samples by Month")

jpeg('/home/webblab/Documents/HP/distribution_plots/simulated_diving_samplenumber.jpeg')
  diving.plot
dev.off()

###################################################################################
#plot distribution of Anserinae data by month and year
###################################################################################
total <- matrix(nrow = 5, ncol= 12)
for (i in 1:5) {
  for (j in 1:12) {
    y <- as.numeric(sim.geese %>% tally(sim.geese$total[sim.geese$month == j & sim.geese$year == i]))
    total[i, j] <- y
  }
}

total.geese <- data.frame(year = as.factor(rep(seq(1,5), each = 12)),
                         month = as.factor(rep(seq(1,12), times = 5)))
total.geese$total <- NA
for (i in 1:5) {
  for (j in 1:12) {
    y <- total[i,j]
    total.geese$total[total.geese$year == i & total.geese$month == j] <- y
  }
}

geese.plot <- ggplot(total.geese, aes(x=year, y = total, fill = month)) +
  geom_bar(stat = "identity", position = "dodge")+
  theme_minimal() +
  scale_fill_manual(values = palette) +
  ylim(0, 3500) +
  ggtitle("Total Number of Anserinae Samples by Month")

jpeg('/home/webblab/Documents/HP/distribution_plots/simulated_anserinae_samplenumber.jpeg')
  geese.plot
dev.off()

##################################################################################
#plot apparent prevalence distribution by month and year
#dabbling ducks
##################################################################################
apparent.mean <- matrix(nrow = 5, ncol = 12)
apparent.sd <- matrix(nrow = 5, ncol = 12)
for (i in 1:5) {
  for (j in 1:12) {
    hold <- sim.mall$apparent.prevalance[sim.mall$year == i & sim.mall$month == j]
    mean <- quantile(hold, probs = 0.5)
    stand <- sd(hold)
    apparent.mean[i, j] <- mean
    apparent.sd[i, j] <- stand
  }
}

apparent.mall <- data.frame(year = as.factor(rep(seq(1,5), each = 12)), 
                            month = as.factor(rep(seq(1, 12), times = 5)),
                            mean = NA, stand.dev = NA, ymin = NA, ymax = NA)
for (i in 1:5) {
  for (j in 1:12) {
    apparent.mall$mean[apparent.mall$year == i & apparent.mall$month == j] <-
      apparent.mean[i, j]
    apparent.mall$stand.dev[apparent.mall$year == i & apparent.mall$month == j] <-
      apparent.sd[i, j]
    apparent.mall$ymin[apparent.mall$year == i & apparent.mall$month == j] <-
      (apparent.mean[i, j] - apparent.sd[i,j])
    apparent.mall$ymax[apparent.mall$year == i & apparent.mall$month == j] <-
      (apparent.mean[i, j] + apparent.sd[i, j])
  }
}
for (i in 1:5) {
  for (j in 1:12) {
    if (apparent.mall$ymin[apparent.mall$year == i & apparent.mall$month == j] < 0) {
      apparent.mall$ymin[apparent.mall$year == i & apparent.mall$month == j] <- 0
    } else {
      apparent.mall$ymin[apparent.mall$year == i & apparent.mall$month == j] <-
        apparent.mall$ymin[apparent.mall$year == i & apparent.mall$month == j]
    }
  }
}

mall.apparent.plot <- ggplot(apparent.mall, aes(x=month, y = mean)) +
  geom_point(color = "darkblue") +
  geom_errorbar(aes( ymin = ymin, ymax = ymax), color = "darkblue") +
  facet_grid(. ~ year) +
  theme_minimal() +
  ylim(0, 0.3) +
  ggtitle("Dabbling Duck Mean Apparent Prevalence")
  
jpeg('/home/webblab/Documents/HP/distribution_plots/simulated_dabbling_appprev.jpeg')
  mall.apparent.plot
dev.off()

##################################################################################
#plot apparent prevalence distribution by month and year
#diving ducks
##################################################################################
apparent.mean <- matrix(nrow = 5, ncol = 12)
apparent.sd <- matrix(nrow = 5, ncol = 12)
for (i in 1:5) {
  for (j in 1:12) {
    hold <- sim.diving$apparent.prevalance[sim.diving$year == i & sim.diving$month == j]
    mean <- quantile(hold, probs = 0.5)
    stand <- sd(hold)
    apparent.mean[i, j] <- mean
    apparent.sd[i, j] <- stand
  }
}

apparent.diving <- data.frame(year = as.factor(rep(seq(1,5), each = 12)), 
                            month = as.factor(rep(seq(1, 12), times = 5)),
                            mean = NA, stand.dev = NA, ymin = NA, ymax = NA)
for (i in 1:5) {
  for (j in 1:12) {
    apparent.diving$mean[apparent.diving$year == i & apparent.diving$month == j] <-
      apparent.mean[i, j]
    apparent.diving$stand.dev[apparent.diving$year == i & apparent.diving$month == j] <-
      apparent.sd[i, j]
    apparent.diving$ymin[apparent.diving$year == i & apparent.diving$month == j] <-
      (apparent.mean[i, j] - apparent.sd[i,j])
    apparent.diving$ymax[apparent.diving$year == i & apparent.diving$month == j] <-
      (apparent.mean[i, j] + apparent.sd[i, j])
  }
}
for (i in 1:5) {
  for (j in 1:12) {
    if (apparent.diving$ymin[apparent.diving$year == i & apparent.diving$month == j] < 0) {
      apparent.diving$ymin[apparent.diving$year == i & apparent.diving$month == j] <- 0
    } else {
      apparent.diving$ymin[apparent.diving$year == i & apparent.diving$month == j] <-
        apparent.diving$ymin[apparent.diving$year == i & apparent.diving$month == j]
    }
  }
}

diving.apparent.plot <- ggplot(apparent.diving, aes(x=month, y = mean)) +
  geom_point(color = "darkblue") +
  geom_errorbar(aes( ymin = ymin, ymax = ymax), color = "darkblue") +
  facet_grid(. ~ year) +
  theme_minimal() +
  ylim(0, 0.3) +
  ggtitle("Diving Duck Mean Apparent Prevalence")

jpeg('/home/webblab/Documents/HP/distribution_plots/simulated_diving_appprev.jpeg')
  diving.apparent.plot
dev.off()

##################################################################################
#plot apparent prevalence distribution by month and year
#Anserinae
##################################################################################
apparent.mean <- matrix(nrow = 5, ncol = 12)
apparent.sd <- matrix(nrow = 5, ncol = 12)
for (i in 1:5) {
  for (j in 1:12) {
    hold <- sim.geese$apparent.prevalance[sim.geese$year == i & sim.geese$month == j]
    mean <- quantile(hold, probs = 0.5)
    stand <- sd(hold)
    apparent.mean[i, j] <- mean
    apparent.sd[i, j] <- stand
  }
}

apparent.geese <- data.frame(year = as.factor(rep(seq(1,5), each = 12)), 
                            month = as.factor(rep(seq(1, 12), times = 5)),
                            mean = NA, stand.dev = NA, ymin = NA, ymax = NA)
for (i in 1:5) {
  for (j in 1:12) {
    apparent.geese$mean[apparent.geese$year == i & apparent.geese$month == j] <-
      apparent.mean[i, j]
    apparent.geese$stand.dev[apparent.geese$year == i & apparent.geese$month == j] <-
      apparent.sd[i, j]
    apparent.geese$ymin[apparent.geese$year == i & apparent.geese$month == j] <-
      (apparent.mean[i, j] - apparent.sd[i,j])
    apparent.geese$ymax[apparent.geese$year == i & apparent.geese$month == j] <-
      (apparent.mean[i, j] + apparent.sd[i, j])
  }
}
for (i in 1:5) {
  for (j in 1:12) {
    if (apparent.geese$ymin[apparent.geese$year == i & apparent.geese$month == j] < 0) {
      apparent.geese$ymin[apparent.geese$year == i & apparent.geese$month == j] <- 0
    } else {
      apparent.geese$ymin[apparent.geese$year == i & apparent.geese$month == j] <-
        apparent.geese$ymin[apparent.geese$year == i & apparent.geese$month == j]
    }
  }
}

geese.apparent.plot <- ggplot(apparent.geese, aes(x=month, y = mean)) +
  geom_point(color = "darkblue") +
  geom_errorbar(aes( ymin = ymin, ymax = ymax), color = "darkblue") +
  facet_grid(. ~ year) +
  theme_minimal() +
  ylim(0, 0.3) +
  ggtitle("Anserinae Mean Apparent Prevalence")

jpeg('/home/webblab/Documents/HP/distribution_plots/simulated_anserinae_appprev.jpeg')
  geese.apparent.plot
dev.off()

#####################################################################################
#number of samples by huc
#####################################################################################

nhuc <- 222

spatial <- matrix(nrow = 3, ncol = nhuc)
for (i in 1:3) {
  for (j in 1:nhuc) {
    spatial [i, j] <- sum(simulated$total[simulated$species.group == i & 
                                            simulated$huc == j])
  }
}

samples.huc <- data.frame(species.group = as.factor(rep(seq(1,3), each = nhuc)),
                          huc = seq(1,nhuc), samples = NA)

for (i in 1:3) {
  for (j in 1:nhuc) {
    samples.huc$samples[samples.huc$species.group == i & samples.huc$huc == j] <-
      spatial[i, j]
  }
} 

palette <- brewer.pal(3, "YlGnBu")

samples.plot <- ggplot(samples.huc, aes(x=huc, y=samples, fill = species.group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  scale_fill_manual(values = palette) +
  facet_grid(. ~ species.group) +
  ggtitle("Total Number of Samples by Huc")

jpeg('home/webblab/Documents/HP/distribution_plots/simulated_totalsamples.jpeg')
  samples.plot
dev.off()