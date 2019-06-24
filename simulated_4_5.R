######################################################################################
######################################################################################
#Attempt 4 and 5 data simulation
#decreased number of total samples (4) and decreased apparent prevalence (5)
######################################################################################
######################################################################################

long <- readRDS("/home/webblab/Documents/HP/sim_data_2.rds")

long$halved.samples <- NA
long$halved.prev <- NA
for (i in 1:length(long$app.prev)) {
  long$halved.samples[i] <- 0.5*(long$total.samples[i])
  long$halved.samples[i] <- round(long$halved.samples[i])
  long$halved.prev[i] <- 0.5*(long$app.prev[i])
}

long$sim.pos.4 <- NA
long$sim.app.prev.4 <- NA
for (i in 1:length(long$app.prev)) {
  hold <- rbinom(long$halved.samples[i], 1, long$app.prev[i])
  long$sim.pos.4[i] <- length(hold[hold == 1])
  long$sim.app.prev.4[i] <- long$sim.pos.4[i]/long$halved.samples[i]
}

long$sim.pos.5 <- NA
long$sim.app.prev.5 <- NA
for (i in 1:length(long$app.prev)) {
  hold <- rbinom(long$total.samples[i], 1, long$halved.prev[i])
  long$sim.pos.5[i] <- length(hold[hold == 1])
  long$sim.app.prev.5[i] <- long$sim.pos.5[i]/long$total.samples[i]
}

#remove watersheds with no samples, but maintain species groups
remove <- long[!(long$total.samples == 0), ]
remove.mall <- remove[remove$species.group == 1, ]
shortened.mall <- long[long$species.group ==1 & long$huc %in% unique(remove.mall$huc), ]
remove.diving <- remove[remove$species.group ==2, ]
shortened.diving <- long[long$species.group ==2 & long$huc %in% unique(remove.diving$huc), ]
remove.geese <- long[long$species.group == 3, ]
shortened.geese <- long[long$species.group == 3 & long$huc %in% unique(remove.geese$huc), ]

#bind seperate species group shortened sets into a single shortened dataset
shortened <- rbind(shortened.mall, shortened.diving, shortened.geese)

#save attempt 3
saveRDS(shortened, "/home/webblab/Documents/HP/sim_data_4_5.rds")

####################################################################################
#create n and y matrices for Attempt 4
#species.group ==1 
####################################################################################

nhucs = length(unique(shortened$huc[shortened$species.group ==1]))

#renumber hucs to go in chronological order for model fitting
number <- data.frame(huc = unique(shortened$huc[shortened$species.group == 1]),
                     new = seq(1,nhucs))
# shortened.mall <- shortened[shortened$species.group == 1, ]
for (i in 1:length(shortened.mall$huc)) {
  shortened.mall$huc.new[i] <- number$new[number$huc == shortened.mall$huc[i]]
}

n <- array(NA, dim = c(12, 5, nhucs))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:nhucs) {
      n[i, j, k] <- shortened.mall$halved.samples[shortened.mall$species.group == 1 & shortened.mall$month == i & 
                                                   shortened.mall$year == j & shortened.mall$huc.new == k]
    }
  }
}

y <- array(NA, dim = c(12, 5, nhucs))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:nhucs) {
      y[i, j, k] <- shortened.mall$sim.pos.4[shortened.mall$species.group == 1 & shortened.mall$month == i & 
                                             shortened.mall$year == j & shortened.mall$huc.new == k]
    }
  }
}

saveRDS(n, "/home/webblab/Documents/HP/sim_data_n_4.rds")
saveRDS(y, "/home/webblab/Documents/HP/sim_data_y_4.rds")

####################################################################################
#create n and y matrices for Attempt 5
#species.group ==1 
####################################################################################
# nhucs = length(unique(shortened$huc[shortened$species.group ==1]))
# 
# #renumber hucs to go in chronological order for model fitting
# number <- data.frame(huc = unique(shortened$huc[shortened$species.group == 1]),
#                      new = seq(1,nhucs))
# # shortened.mall <- shortened[shortened$species.group == 1, ]
# for (i in 1:length(shortened.mall$huc)) {
#   shortened.mall$huc.new[i] <- number$new[number$huc == shortened.mall$huc[i]]
# }

n <- array(NA, dim = c(12, 5, nhucs))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:nhucs) {
      n[i, j, k] <- shortened.mall$total.samples[shortened.mall$species.group == 1 & shortened.mall$month == i & 
                                                    shortened.mall$year == j & shortened.mall$huc.new == k]
    }
  }
}

y <- array(NA, dim = c(12, 5, nhucs))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:nhucs) {
      y[i, j, k] <- shortened.mall$sim.pos.5[shortened.mall$species.group == 1 & shortened.mall$month == i & 
                                               shortened.mall$year == j & shortened.mall$huc.new == k]
    }
  }
}

saveRDS(n, "/home/webblab/Documents/HP/sim_data_n_5.rds")
saveRDS(y, "/home/webblab/Documents/HP/sim_data_y_5.rds")