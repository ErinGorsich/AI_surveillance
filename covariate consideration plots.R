##############################################################################
##############################################################################
#plots to determine if a species covariate is needed
##############################################################################
##############################################################################


library(ggplot2)
library(RColorBrewer)
library(gridExtra)

setwd("~/Github")
data <- readRDS("locationsamplingevent_n_y_speciesgroup.rds")

#limit data to dabbling ducks, diving ducks, and geese and swans
data <- data[data$species.group %in% c("1", "2", "4"), ]
data$species.group[data$species.group == "4"] <- 3 
#convert years into numbers 1-5
data$year.number <- NA
for (i in 1:length(data$year)) {
  if (data$year[i] == "2007") {
    data$year.number[i] <- 1
  } else if (data$year[i] == "2008") {
    data$year.number[i] <- 2
  } else if (data$year[i] == "2009") {
    data$year.number[i] <- 3
  } else if (data$year[i] == "2010") {
    data$year.number[i] <- 4
  } else {
    data$year.number[i] <- 5
  }
}  

data$month.number <- as.integer(data$month)
data$huc.number <- as.(data$watershed)

#seperate species groups
dabbling <- data[data$species.group == "1", ]
diving <- data[data$species.group == "2", ]
geese <- data[data$species.group == "3", ]

nhucs <- length(unique(data$watershed))
nyears <- length(unique(data$year))
nmonths <- length(unique(data$month))

#create a matrix with the total number of samples for each watershed in each month for each species
total.dabbling <- matrix(data = NA, nrow = nhucs, ncol = nyears*nmonths)
positive.dabbling <- matrix(data = NA, nrow = nhucs, ncol = nyears*nmonths)
for (i in 1:length(unique(data$month.number))) {
  for (j in 1:length(unique(data$year.number))) {
    for (k in 1:length(unique(data$watershed))) {
      m <- unique(data$watershed)[k]
      l <- (nmonths*(j-1))+i
      total.dabbling[k, l] <- sum(dabbling$n[dabbling$month.number == i & dabbling$year.number == j & dabbling$watershed == m])
      positive.dabbling[k, l] <- sum(dabbling$y[dabbling$month.number == i & dabbling$year.number == j & dabbling$watershed == m])
    }
  }
}  

total.diving <- matrix(data = NA, nrow = nhucs, ncol = nyears*nmonths)
positive.diving <- matrix(data = NA, nrow = nhucs, ncol = nyears*nmonths)
for (i in 1:length(unique(data$month))) {
  for (j in 1:length(unique(data$year.number))) {
    for (k in 1:length(unique(data$watershed))) {
      m <- unique(data$watershed)[k]
      l <- (nmonths*(j-1))+i
      total.diving[k, l] <- sum(diving$n[diving$month == i & diving$year.number == j & diving$watershed == m])
      positive.diving[k, l] <- sum(diving$y[diving$month == i & diving$year.number == j & diving$watershed == m])
    }
  }
}  

total.geese <- matrix(data = NA,  nrow = nhucs, ncol = nyears*nmonths)
positive.geese <- matrix(data = NA,  nrow = nhucs, ncol = nyears*nmonths)
for (i in 1:length(unique(data$month))) {
  for (j in 1:length(unique(data$year.number))) {
    for (k in 1:length(unique(data$watershed))) {
      m <- unique(data$watershed)[k]
      l <- (nmonths*(j-1))+i
      total.geese[k, l] <- sum(geese$n[geese$month == i & geese$year.number == j & geese$watershed == m])
      positive.geese[k, l] <- sum(geese$y[geese$month == i & geese$year.number == j & geese$watershed == m])
    }
  }
}  

#create a data frame for plotting with ggplot2

watershed <- data.frame(watershed = rep(unique(data$watershed), each=nmonths, times = nyears),
                        month = rep(seq(1,12), times = nhucs), year = rep(seq(1,5), each = nmonths*nhucs),
                        total.dabbling = NA, positive.dabbling = NA, apparent.dabbling = NA,
                        total.diving = NA, positive.diving = NA, apparent.diving = NA,
                        total.geese = NA, positive.geese = NA, apparent.geese = NA)

for (i in 1:nmonths) {
  for (j in 1:nyears) {
    for (k in 1:nhucs) {
      m <- unique(data$watershed)[k]
      watershed$total.dabbling[watershed$watershed == m & watershed$month == i & watershed$year == j] <- total.dabbling[k, (nmonths*(j-1))+i]
      watershed$total.diving[watershed$watershed == m & watershed$month == i & watershed$year == j] <- total.diving[k, (nmonths*(j-1))+i]
      watershed$total.geese[watershed$watershed == m & watershed$month == i & watershed$year == j] <- total.geese[k, (nmonths*(j-1))+i]
      watershed$positive.dabbling[watershed$watershed == m & watershed$month == i & watershed$year == j] <- positive.dabbling[k, (nmonths*(j-1))+i]
      watershed$positive.diving[watershed$watershed == m & watershed$month == i & watershed$year == j] <- positive.diving[k, (nmonths*(j-1))+i]
      watershed$positive.geese[watershed$watershed == m & watershed$month == i & watershed$year == j] <- positive.geese[k, (nmonths*(j-1))+i]
    }
  }
}

watershed <- watershed[!(watershed$total.dabbling == 0 & watershed$positive.dabbling == 0 
                          & watershed$total.diving == 0 & watershed$positive.diving == 0
                          & watershed$total.geese == 0 & watershed$positive.geese == 0), ]
watershed$apparent.dabbling <- watershed$positive.dabbling/watershed$total.dabbling
watershed$apparent.diving <- watershed$positive.diving/watershed$total.diving
watershed$apparent.geese <- watershed$positive.geese/watershed$total.geese

#make a plot

n <- length(unique(watershed$watershed))
ball <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pal <- sample(ball, n)

dabbling.diving <- ggplot(data = watershed) +
                    geom_point(aes(x=apparent.dabbling, y=apparent.diving, color = watershed)) +
                    scale_color_manual(values = pal) +
                    xlab("Dabbling Duck Apparent Prevlence") +
                    ylab("Diving Duck Apparent Prevalence") +
                    theme(panel.background=element_blank(),
                    axis.line = element_line(),
                    legend.position = "none")

dabbling.geese <- ggplot(data = watershed) +
                    geom_point(aes(x=apparent.dabbling, y=apparent.geese, color = watershed)) +
                    scale_color_manual(values = pal) +
                    xlab("Dabbling Duck Apparent Prevalence") +
                    ylab("Geese and Swan Apparent Prevalence") +
                    theme(panel.background=element_blank(),
                    axis.line = element_line(),
                    legend.position = "none")
jpeg("~/Honors Thesis/Thesis Work/Figures/Covariate Determination.jpeg")
grid.arrange(dabbling.diving, dabbling.geese, ncol = 2)
dev.off()

cor(watershed$apparent.dabbling, watershed$apparent.diving, use = "complete.obs")
cor(watershed$apparent.dabbling, watershed$apparent.geese, use = "complete.obs")
