###################################################################################
###################################################################################
#Check distribution of actual data across watersheds
###################################################################################
###################################################################################
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

ai <- readRDS("/home/webblab/Documents/HP/locationsamplingevent_n_y_speciesgroup.rds")
ai <- ai[ai$species.group %in% c(1,2,4), ]
ai$sample.event <- NULL
  
ai$year.number <- NA
for (i in 1:length(ai$year)){ #length(ai$year)
  if (ai$year[i] == '2007') {
    ai$year.number[i] <- 1
  } else if (ai$year[i] == '2008') {
    ai$year.number[i] <- 2
  } else if (ai$year[i] == '2009') {
    ai$year.number[i] <- 3
  } else if (ai$year[i] == '2010') {
    ai$year.number[i] <- 4
  } else {
    ai$year.number[i] <- 5
  }
}
ai$year <- NULL

for (i in 1:length(ai$species.group)) { #length(ai$species.group)
  if (ai$species.group[i] == 4) {
    ai$species.group[i] <- 3
  }
}

ai$huc <- as.numeric(ai$watershed)
ai$watershed <- NULL

tally <- data.frame(species.group = as.factor(rep(seq(1,3), each = 60*222)),
                    huc = as.factor(rep(seq(1,222), each = 60, times = 3)),
                    year = as.factor(rep(seq(1,5), each = 12, times = 666)),
                    month = as.factor(rep(seq(1,12), times = 5*222)),
                    total = NA, positive = NA, apparent.prev = NA)

total <- ai %>% group_by (species.group, huc, month, year.number)
total$month <- as.numeric(total$month)
for (i in 1:3) {
  for (j in 1:222) {
    for (k in 1:5) {
      for (l in 1:12) {
        tally$total[tally$species.group == i & tally$huc == j & tally$year == k & 
                    tally$month == l] <-
          sum(total$n[total$species.group == i & total$huc == j & total$year.number == k &
                      total$month == l])
        tally$positive[tally$species.group == i & tally$huc == j & tally$year == k &
                      tally$month == l] <-
          sum(total$y[total$species.group == i & total$huc == j & total$year.number == k &
                        total$month == l])
      }  
    }
  }
}

tally$apparent.prev <- tally$positive/tally$total
tally$apparent.prev[tally$apparent.prev == "NaN"] <- 0

real.mall <- ai[ai$species.group == 1, ]
real.diving <- ai[ai$species.group ==2, ]
real.geese <- ai[ai$species.group == 3, ]
palette <- brewer.pal(12, "Paired")

#################################################################################
#plot distribution of dabbling duck data by month and year
#################################################################################
real.mall$species.group <- NULL
real.mall$month <- as.numeric(real.mall$month)
total <- matrix(nrow = 5, ncol = 12, NA)
for (i in 1:5) {
  for (j in 1:12) {
    total[i,j] <- sum(real.mall$n[real.mall$year.number == i & real.mall$month == j])
  }
}

total.mall <- data.frame(year = as.factor(rep(seq(1,5), each =12)),
                         month = as.factor(rep(seq(1,12), times = 5)),
                         total = NA)

for (i in 1:5) {
  for (j in 1:12) {
    total.mall$total[total.mall$year == i & total.mall$month == j] <- total[i,j]
  }
}

real.mall.plot <- ggplot(total.mall, aes(x=year, y = total, fill = month)) +
  geom_bar(stat = "identity", position = "dodge")+
  theme_minimal() +
  scale_fill_manual(values = palette) +
  ylim(0, 10000) +
  ggtitle("Total Number of Dabbling Duck Samples by Month (Actual Data)")

jpeg('/home/webblab/Documents/HP/distribution_plots/real_dabbling_samplenumber.jpeg')
  real.mall.plot
dev.off()

#################################################################################
#plot distribution of diving duck data by month and year
#################################################################################

real.diving$species.group <- NULL
real.diving$month <- as.numeric(real.diving$month)

total <- matrix(nrow = 5, ncol = 12, NA)
for (i in 1:5) {
  for (j in 1:12) {
    total[i,j] <- sum(real.diving$n[real.diving$year.number == i & real.diving$month == j])
  }
}

total.diving <- data.frame(year = as.factor(rep(seq(1,5), each =12)),
                         month = as.factor(rep(seq(1,12), times = 5)),
                         total = NA)

for (i in 1:5) {
  for (j in 1:12) {
    total.diving$total[total.diving$year == i & total.diving$month == j] <- total[i,j]
  }
}

real.diving.plot <- ggplot(total.diving, aes(x=year, y = total, fill = month)) +
  geom_bar(stat = "identity", position = "dodge")+
  theme_minimal() +
  scale_fill_manual(values = palette) +
  ylim(0, 1000) +
  ggtitle("Total Number of Diving Duck Samples by Month (Actual Data)")

jpeg('/home/webblab/Documents/HP/distribution_plots/real_diving_samplenumber.jpeg')
  real.diving.plot
dev.off()

#################################################################################
#plot distribution of Anserinae data by month and year
#################################################################################

real.geese$species.group <- NULL
real.geese$month <- as.numeric(real.geese$month)

total <- matrix(nrow = 5, ncol = 12, NA)
for (i in 1:5) {
  for (j in 1:12) {
    total[i,j] <- sum(real.geese$n[real.geese$year.number == i & real.geese$month == j])
  }
}

total.geese <- data.frame(year = as.factor(rep(seq(1,5), each =12)),
                           month = as.factor(rep(seq(1,12), times = 5)),
                           total = NA)

for (i in 1:5) {
  for (j in 1:12) {
    total.geese$total[total.geese$year == i & total.geese$month == j] <- total[i,j]
  }
}

real.geese.plot <- ggplot(total.geese, aes(x=year, y = total, fill = month)) +
  geom_bar(stat = "identity", position = "dodge")+
  theme_minimal() +
  scale_fill_manual(values = palette) +
  ylim(0, 3500) +
  ggtitle("Total Number of Anserinae Samples by Month (Actual Data)")

jpeg('/home/webblab/Documents/HP/distribution_plots/real_anserinae_samplenumber.jpeg')
  real.geese.plot
dev.off()

######################################################################################
#plot apparent prevalence of dabbling duck data by month and year
##################################################################################

tally.mall <- tally[tally$species.group == 1, ]

mean.apparent <- matrix(nrow = 5, ncol = 12)
sd.apparent <- matrix(nrow = 5, ncol = 12)
for (i in 1:5) {
  for (j in 1:12) {
    hold <- tally.mall$apparent.prev[tally.mall$year == i & tally.mall$month == j]
    mean.apparent[i, j] <- quantile(hold, probs = 0.5)
    sd.apparent[i, j] <- sd(hold)
  }
}

apparent.mall <- data.frame(year = as.factor(rep(seq(1,5), each = 12)), 
                            month = as.factor(rep(seq(1, 12), times = 5)),
                            mean = NA, stand.dev = NA, ymin = NA, ymax = NA)
for (i in 1:5) {
  for (j in 1:12) {
    apparent.mall$mean[apparent.mall$year == i & apparent.mall$month == j] <-
      mean.apparent[i, j]
    apparent.mall$stand.dev[apparent.mall$year == i & apparent.mall$month == j] <-
      sd.apparent[i, j]
    apparent.mall$ymin[apparent.mall$year == i & apparent.mall$month == j] <-
      (mean.apparent[i, j] - sd.apparent[i,j])
    apparent.mall$ymax[apparent.mall$year == i & apparent.mall$month == j] <-
      (mean.apparent[i, j] + sd.apparent[i, j])
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

real.mall.apparent.plot <- ggplot(apparent.mall, aes(x=month, y = mean)) +
  geom_point(color = "darkred") +
  geom_errorbar(aes( ymin = ymin, ymax = ymax), color = "darkred") +
  facet_grid(. ~ year) +
  theme_minimal() +
  ylim(0, 0.3) +
  ggtitle("Dabbling Duck Mean Apparent Prevalence (Actual Data)")

jpeg('/home/webblab/Documents/HP/distribution_plots/real_dabbling_appprev.jpeg')
  real.mall.apparent.plot
dev.off()

######################################################################################
#plot apparent prevalence of diving duck data by month and year
##################################################################################

tally.diving <- tally[tally$species.group == 2, ]

mean.apparent <- matrix(nrow = 5, ncol = 12)
sd.apparent <- matrix(nrow = 5, ncol = 12)
for (i in 1:5) {
  for (j in 1:12) {
    hold <- tally.diving$apparent.prev[tally.diving$year == i & tally.diving$month == j]
    mean.apparent[i, j] <- quantile(hold, probs = 0.5)
    sd.apparent[i, j] <- sd(hold)
  }
}

apparent.diving <- data.frame(year = as.factor(rep(seq(1,5), each = 12)), 
                            month = as.factor(rep(seq(1, 12), times = 5)),
                            mean = NA, stand.dev = NA, ymin = NA, ymax = NA)
for (i in 1:5) {
  for (j in 1:12) {
    apparent.diving$mean[apparent.diving$year == i & apparent.diving$month == j] <-
      mean.apparent[i, j]
    apparent.diving$stand.dev[apparent.diving$year == i & apparent.diving$month == j] <-
      sd.apparent[i, j]
    apparent.diving$ymin[apparent.diving$year == i & apparent.diving$month == j] <-
      (mean.apparent[i, j] - sd.apparent[i,j])
    apparent.diving$ymax[apparent.diving$year == i & apparent.diving$month == j] <-
      (mean.apparent[i, j] + sd.apparent[i, j])
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

real.diving.apparent.plot <- ggplot(apparent.diving, aes(x=month, y = mean)) +
  geom_point(color = "darkred") +
  geom_errorbar(aes( ymin = ymin, ymax = ymax), color = "darkred") +
  facet_grid(. ~ year) +
  theme_minimal() +
  ylim(0, 0.3) +
  ggtitle("Diving Duck Mean Apparent Prevalence (Actual Data)")

jpeg('/home/webblab/Documents/HP/distribution_plots/real_diving_appprev.jpeg')
  real.diving.apparent.plot
dev.off()

######################################################################################
#plot apparent prevalence of Anserinae data by month and year
##################################################################################

tally.geese <- tally[tally$species.group == 3, ]

mean.apparent <- matrix(nrow = 5, ncol = 12)
sd.apparent <- matrix(nrow = 5, ncol = 12)
for (i in 1:5) {
  for (j in 1:12) {
    hold <- tally.geese$apparent.prev[tally.geese$year == i & tally.geese$month == j]
    mean.apparent[i, j] <- quantile(hold, probs = 0.5)
    sd.apparent[i, j] <- sd(hold)
  }
}

apparent.geese <- data.frame(year = as.factor(rep(seq(1,5), each = 12)), 
                              month = as.factor(rep(seq(1, 12), times = 5)),
                              mean = NA, stand.dev = NA, ymin = NA, ymax = NA)
for (i in 1:5) {
  for (j in 1:12) {
    apparent.geese$mean[apparent.geese$year == i & apparent.geese$month == j] <-
      mean.apparent[i, j]
    apparent.geese$stand.dev[apparent.geese$year == i & apparent.geese$month == j] <-
      sd.apparent[i, j]
    apparent.geese$ymin[apparent.geese$year == i & apparent.geese$month == j] <-
      (mean.apparent[i, j] - sd.apparent[i,j])
    apparent.geese$ymax[apparent.geese$year == i & apparent.geese$month == j] <-
      (mean.apparent[i, j] + sd.apparent[i, j])
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

real.geese.apparent.plot <- ggplot(apparent.diving, aes(x=month, y = mean)) +
  geom_point(color = "darkred") +
  geom_errorbar(aes( ymin = ymin, ymax = ymax), color = "darkred") +
  facet_grid(. ~ year) +
  theme_minimal() +
  ylim(0, 0.3) +
  ggtitle("Anserinae Mean Apparent Prevalence (Actual Data)")

jpeg('/home/webblab/Documents/HP/distribution_plots/real_anserinae_appprev.jpeg')
  real.geese.apparent.plot
dev.off()

#####################################################################################
#number of samples by huc and species group
#####################################################################################

nhuc <- length(unique(ai$huc))

real.spatial <- matrix(nrow = 3, ncol = nhuc)
for (i in 1:3) {
  for (j in 1:nhuc) {
    spatial[i,j] <- sum(tally$total[tally$species.group == i & tally$huc == j])
  }
}

real.samples.huc <- data.frame(species.group = as.factor(rep(seq(1,3), each = nhuc)),
                               huc = seq(1, nhuc), samples = NA)
for (i in 1:3) {
  for (j in 1:nhuc) {
    real.samples.huc$samples[real.samples.huc$species.group == i & 
                               real.samples.huc$huc == j] <- spatial[i,j]
  }
}

real.samples.plot <- ggplot(real.samples.huc, aes(x = huc, y = samples, 
                                                  fill = species.group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  scale_fill_manual(values = c("deepskyblue4", "deeppink4", "mediumpurple4")) +
  facet_grid(. ~ species.group) +
  ggtitle("Total Number of Samples by Huc (Actual Data)")

jpeg('home/webblab/Documents/HP/distribution_plots/real_totalsample.jpeg')
  real.samples.plot
dev.off()