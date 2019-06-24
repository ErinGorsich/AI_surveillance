###############################################################
###############################################################
#TEST SIMSTUDY PROGRAM#########################################
###############################################################
###############################################################

library(simstudy)
library(dplyr)

see <- readRDS('/home/webblab/Documents/HP/AVHS_samplingevent_speciesgroup.rds')

#create data set definitions
table <- defData(varname = "Species Group", 
                 formula = "0.80;0.05;0.15", dist = "categorical")
table <- defData(table, varname = "huc", formula = "1;194", dist = "uniform")

counts <- table(see$collection.month)
total <- sum(counts)
probs <- matrix(NA, nrow = 1, ncol = 12)
for (i in 1:length(counts)) {
  probs[i] <- counts[i]/total
}
probs
table <- defData(table, varname = "month",
                 formula = "0.125;0.021;0.011;0.005;0.013;0.054;0.041;0.103;0.137;0.187;0.148;0.154",
                 dist = "categorical")

hold <- summary(see$collection.year)
total <- sum(hold)
probs <- matrix(NA, nrow = 1, ncol = 5)
for (i in 1:length(hold)) {
  probs[i] <- hold[i]/total
}
probs
table <- defData(table, varname = "year", formula = "0.219;0.264;0.213;0.135;0.169", 
                 dist = "categorical")
#2007=1, 2008=2, 2009=3, 2010=4, 2015=15

results <- summary(see$AIpcr_susneg)
results <- results[1:2]
total <- sum(results)
probs <- matrix(NA, nrow = 1, ncol = 2)
for (i in 1:length(results)) {
  probs[i] <- results[i]/total
}
probs
table <- defData(table, varname = "result", formula = "0.882;0.118", dist = "categorical")
#positive = 2, negative = 1

#generate a dataset
data <- genData(220000, table)

#trucate huc numbers to whole numbers
data$HUC <- NA
for (i in 1:length(data$huc)) {
  data$HUC[i] <- trunc(data$huc[i])
}

#convert year numbers to actual years
data$YEAR <- NA
for (i in 1:length(data$year)) {
  if (data$year[i] == 1) {
    data$YEAR[i] <- 2007
  } else if (data$year[i] == 2) {
    data$YEAR[i] <- 2008
  } else if (data$year[i] == 3) {
    data$YEAR[i] <- 2009
  } else if (data$year[i] == 4) {
    data$YEAR[i] <- 2010
  } else {
    data$YEAR[i] <- 2015
  }
}

#convert results to positive, negative, unknown
#data$RESULT <- NA
#for (i in 1:length(data$result)) {
#  if (data$result[i] == 1) {
 #   data$RESULT[i] <- "negative"
  #} else if (data$result[i] == 2) {
   # data$RESULT[i] <- "positive" 
#}

#create a table for use in the model code
chart <- data.frame(species.group = rep(seq(1,3), each = 222*5*12),
                    huc = rep(seq(1,222), each = 60, times = 3),
                    month = rep(seq(1,12), times = 5*222*3),
                    year = rep(seq(1,5), each = 12, times = 222*3))
count <- data %>% group_by (Species.Group, HUC, month, year, result) %>% tally()

chart$total <- NA
chart$positive <- NA
for (i in 1:3) {
  for (j in 1:222) {
    for (k in 1:12) {
      for (l in 1:5) {
        chart$positive[chart$species.group==i & chart$huc==j & chart$month==k & chart$year==l] <-
          count$n[count$Species.Group==i & count$HUC==j & count$month==k & count$year==l & count$result == 1]
      }
    }
  }
}

for (i in 1:3) {
  for (j in 1:222) {
    for (k in 1:12) {
      for (l in 1:5) {
        m <- count$n[count$Species.Group==i & count$HUC==j & count$month==k & count$year==l & count$result == 1]
        n <- count$n[count$Species.Group==i & count$HUC==j & count$month==k & count$year==l & count$result == 2]
        if (is.na(m[1]==TRUE)) {
          m[1] <- 0
        } 
        if (is.na(n[1]==TRUE)) {
          n[1] <- 0
        }
        o <- m[1] + n[1]
        chart$positive[chart$species.group==i & chart$huc==j & chart$month==k & chart$year==l] <- n
        chart$total[chart$species.group==i & chart$huc==j & chart$month==k & chart$year==l] <- o
      }
    }
  }
}

chart$apparent.prevalance <- NA
for (i in 1:length(chart$total)) {
  chart$apparent.prevalance[i] <- chart$positive[i]/chart$total[i]
}
chart$apparent.prevalance[is.na(chart$apparent.prevalance)] <- 0

saveRDS(chart, '/home/webblab/Documents/HP/simulated_data.rds')

#create an n (total samples) and y (positive samples) matrix for model fitting
#NO SPECIES GROUP!!!!! set at species.group = 1 (mallards)

n <- array(NA, dim = c(12, 5, 222))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:222) {
      n[i, j, k] <- chart$total[chart$species.group == 1 & chart$month == i & chart$year == j & chart$huc == k]
    }
  }
}

y <- array(NA, dim = c(12, 5, 222))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:222) {
      y[i, j, k] <- chart$positive[chart$species.group == 1 & chart$month == i & chart$year == j & chart$huc == k]
    }
  }
}

saveRDS(n, '/home/webblab/Documents/HP/simulated_n_matrix.rds')
saveRDS(y, '/home/webblab/Documents/HP/simulated_y_matrix.rds')
