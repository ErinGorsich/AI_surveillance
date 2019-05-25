###################################################################################
###################################################################################
#Check distribution of actual data across watersheds
###################################################################################
###################################################################################
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

ai <- readRDS("~/Github/locationsamplingevent_n_y_speciesgroup.rds")
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
  ggtitle("Total Number of Dabbling Duck Samples by Month (Actual Data)")

# jpeg()
#   real.mall.plot
# dev.off()

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
  ggtitle("Total Number of Diving Duck Samples by Month (Actual Data)")

# jpeg()
#   real.diving.plot
# dev.off()

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
  ggtitle("Total Number of Anserinae Samples by Month (Actual Data)")

# jpeg()
#   real.geese.plot
# dev.off()