#################################################################
#################################################################
#simulated data Attempt 3
#remove hucs from long with no data
#################################################################
#################################################################

long <- readRDS("/home/webblab/Documents/HP/sim_data_2.rds")

#create a dataset that has only rows where long$total.samples>0
remove <- long[!(long$total.samples == 0), ]

#need all month/year combinations for hucs with samples in >0 combinations
#need to keep species group seperate so make 3 seperate datasets
#match hucs in remove$hucs to long$hucs

remove.mall <- remove[remove$species.group == 1, ]
shortened.mall <- long[long$species.group ==1 & long$huc %in% unique(remove.mall$huc), ]
remove.diving <- remove[remove$species.group ==2, ]
shortened.diving <- long[long$species.group ==2 & long$huc %in% unique(remove.diving$huc), ]
remove.geese <- long[long$species.group == 3, ]
shortened.geese <- long[long$species.group == 3 & long$huc %in% unique(remove.geese$huc), ]

#bind seperate species group shortened sets into a single shortened dataset
shortened <- rbind(shortened.mall, shortened.diving, shortened.geese)

#save attempt 3
saveRDS(shortened, "/home/webblab/Documents/HP/sim_data_3.rds")

###########################################################################
#create n and y matrices for model fitting
#species.group == 1
###########################################################################
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
      n[i, j, k] <- shortened.mall$total.samples[shortened.mall$species.group == 1 & shortened.mall$month == i & 
                                                   shortened.mall$year == j & shortened.mall$huc.new == k]
    }
  }
}

y <- array(NA, dim = c(12, 5, nhucs))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:nhucs) {
      y[i, j, k] <- shortened.mall$sim.pos[shortened.mall$species.group == 1 & shortened.mall$month == i & 
                                             shortened.mall$year == j & shortened.mall$huc.new == k]
    }
  }
}

saveRDS(n, "/home/webblab/Documents/HP/sim_data_3_n.rds")
saveRDS(y, "/home/webblab/Documents/HP/sim_data_3_y.rds")

###########################################################################
#create n and y matrices for model fitting
#species.group == 3
###########################################################################
nhucs = length(unique(shortened$huc[shortened$species.group ==3]))

#renumber hucs to go in chronological order for model fitting
number <- data.frame(huc = unique(shortened$huc[shortened$species.group == 3]),
                     new = seq(1,nhucs))
shortened.geese <- shortened[shortened$species.group == 3, ]
for (i in 1:length(shortened.geese$huc)) {
  shortened.geese$huc.new[i] <- number$new[number$huc == shortened.geese$huc[i]]
}

n <- array(NA, dim = c(12, 5, nhucs))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:nhucs) {
      n[i, j, k] <- shortened.geese$total.samples[shortened.geese$species.group == 3 & shortened.geese$month == i & 
                                                   shortened.geese$year == j & shortened.geese$huc.new == k]
    }
  }
}

y <- array(NA, dim = c(12, 5, nhucs))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:nhucs) {
      y[i, j, k] <- shortened.geese$sim.pos[shortened.geese$species.group == 3 & shortened.geese$month == i & 
                                             shortened.geese$year == j & shortened.geese$huc.new == k]
    }
  }
}

saveRDS(n, "/home/webblab/Documents/HP/sim_data_3_n_geese.rds")
saveRDS(y, "/home/webblab/Documents/HP/sim_data_3_y_geese.rds")

###########################################################################
#create n and y matrices for model fitting
#species.group == 2
###########################################################################
nhucs = length(unique(shortened$huc[shortened$species.group == 2]))

#renumber hucs to go in chronological order for model fitting
number <- data.frame(huc = unique(shortened$huc[shortened$species.group == 2]),
                     new = seq(1,nhucs))
shortened.diving <- shortened[shortened$species.group == 2, ]
for (i in 1:length(shortened.diving$huc)) {
  shortened.diving$huc.new[i] <- number$new[number$huc == shortened.diving$huc[i]]
}

n <- array(NA, dim = c(12, 5, nhucs))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:nhucs) {
      n[i, j, k] <- shortened.diving$total.samples[shortened.diving$species.group == 2 & shortened.diving$month == i & 
                                                    shortened.diving$year == j & shortened.diving$huc.new == k]
    }
  }
}

y <- array(NA, dim = c(12, 5, nhucs))
for (i in 1:12) {
  for (j in 1:5) {
    for (k  in 1:nhucs) {
      y[i, j, k] <- shortened.diving$sim.pos[shortened.diving$species.group == 2 & shortened.diving$month == i & 
                                              shortened.diving$year == j & shortened.diving$huc.new == k]
    }
  }
}

saveRDS(n, "/home/webblab/Documents/HP/sim_data_3_n_diving.rds")
saveRDS(y, "/home/webblab/Documents/HP/sim_data_3_y_diving.rds")
