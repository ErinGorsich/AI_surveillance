##################################################################################
##################################################################################
#Figures to check assumptions of the models
##################################################################################
##################################################################################

library(ggplot2)
library(grDevices)
library(dplyr)

# setwd("~/Github")
setwd("~/HP/Data")
data <- readRDS("AVHS_samplingevent.rds")

##################################################################################
#assumption 1
#all sampling events within a single huc have true prevalences from the a huc-wide
#common Beta distribution
##################################################################################

#create a data set with sampling event, number of samples, number of positive
#samples, and the apparent prevalence
sampling.events <- data.frame(sampling.event = unique(data$event.number.week),
                              watershed = NA, n = NA, y = NA,
                              apparent.prevalence = NA)

for (i in 1:length(sampling.events$sampling.event)) {
  hold <- data[data$event.number.week == i, ]
  sampling.events[sampling.events$sampling.event == i, ]$watershed <- hold$huc4[1]
  hold$tally <- 1
  sampling.events[sampling.events$sampling.event == i, ]$n <- sum(hold$tally)
  hold.y <- hold[hold$AIpcr_susneg == "positive", ]
  sampling.events[sampling.events$sampling.event == i, ]$y <- sum(hold.y$tally)
}

sampling.events$apparent.prevalence <- sampling.events$y/sampling.events$n

#create a data set with watershed, number of samples, number of positive
#samples, and the apparent prevalence

data$huc4.b <- as.integer(data$huc4)
watershed <- data.frame(huc = unique(data$huc4), watershed = unique(data$huc4.b), 
                        n = NA, y = NA, apparent.prevalence = NA)

for ( i in 1:length(watershed$watershed)) {
  j <- unique(watershed$watershed[i])
  hold <- data[data$huc4.b == j, ]
  hold$tally <- 1
  watershed[watershed$watershed == j, ]$n <- sum(hold$tally)
  hold.y <- hold[hold$AIpcr_susneg == "positive", ]
  watershed[watershed$watershed == j, ]$y <- sum(hold.y$tally)
}

watershed$apparent.prevalence <- watershed$y/watershed$n

#merge dataframes
temp <- data.frame(watershed = watershed$huc, 
                   app.prev.huc = watershed$apparent.prevalence)
temp.2 <- data.frame(watershed = sampling.events$watershed, 
                     sampling.event = sampling.events$sampling.event,
                     app.prev.se = sampling.events$apparent.prevalence)
plot <- inner_join(temp, temp.2, by="watershed")

n <- length(unique(plot$watershed))
ball <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pal <- sample(ball, n)

ggplot(data=plot, aes(x=app.prev.huc, y=app.prev.se)) + 
  geom_point(aes(color = watershed)) +
  scale_colour_manual(values = pal) +
  ylab("Sampling Event Apparent Prevalence") +
  xlab("Watershed Apparent Prevalence") +
  theme(panel.background=element_blank(),
        axis.line = element_line(),
        legend.position = "none")

ggplot(data=sampling.events, aes(x=n, y=apparent.prevalence)) +
  geom_point(aes(color = sampling.event)) +
  scale_colour_distiller(type = "div", palette = "Spectral") +
  xlab("Total Number of Samples") +
  ylab("Sampling Event Apparent Prevalence")+
  theme(panel.background = element_blank(), axis.line = element_line(), 
        legend.position = "none")
