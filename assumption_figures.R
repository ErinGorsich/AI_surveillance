##################################################################################
##################################################################################
#Figures to check assumptions of the models
##################################################################################
##################################################################################

library(ggplot2)

setwd("~/Github")
data <- readRDS("AVHS_samplingevent.rds")

##################################################################################
#assumption 1
#all sampling events within a single huc have true prevalences from the a huc-wide
#common Beta distribution
##################################################################################

#create a data set with sampling event, number of samples, number of positive
#samples, and the apparent prevalence
sampling.events <- data.frame(sampling.event = unique(data$event.number.week),
                              n = NA, y = NA, apparent.prevalence = NA)

for (i in 1:length(sampling.events$sampling.event)) {
  hold <- data[data$event.number.week == i, ]
  hold$tally <- 1
  sampling.events[sampling.events$sampling.event == i, ]$n <- sum(hold$tally)
  hold.y <- hold[hold$AIpcr_susneg == "positive", ]
  sampling.events[sampling.events$sampling.event == i, ]$y <- sum(hold.y$tally)
}

sampling.events$apparent.prevalence <- sampling.events$y/sampling.events$n

# data$apparent.prevalence <- data$y/data$n
# data[is.na(data)] <- 0
# 
# ggplot(data = data[data$watershed=='0316', ], aes(x=month, y=apparent.prevalence)) +
#   geom_point(aes(color = sample.event))
