#######################################################################
#######################################################################
#Sensitivity and Specificity results (ESA Poster)
#######################################################################
#######################################################################

library(ggplot2)

#read in fit data from selected model
fit <- readRDS("~/Github/queens_icar_fit.rds")

dabbling.se.median <- quantile.mcmc.list(fit[,1])[3]
dabbling.se.95 <- quantile.mcmc.list(fit[,1])[5]
dabbling.se.min <- quantile.mcmc.list(fit[,1])[1]
dabbling.sp.median <- quanitle.mcmc.list(fit[,4])[3]
dabbling.sp.95 <- quantile.mcmc.list(fit[,4])[5]
dabbling.sp.min <- quantile.mcmc.list(fit[,4])[1]

diving.se.median <- quantile.mcmc.list(fit[,2])[3]
diving.se.95 <- quantile.mcmc.list(fit[,2])[5]
diving.se.min <- quantile.mcmc.list(fit[,2])[1]
diving.sp.median <- quantile.mcmc.list(fit[,5])[3]
diving.sp.95 <- quantile.mcmc.list(fit[,5])[5]
diving.sp.min <- quantile.mcmc.list(fit[,5])[1]

geese.se.median <- quantile.mcmc.list(fit[,3])[3]
geese.se.95 <- quantile.mcmc.list(fit[,3])[5]
geese.se.min <- quantile.mcmc.list(fit[,3])[1]
geese.sp.median <- quantile.mcmc.list(fit[,6])[3]
geese.sp.95 <- quantile.mcmc.list(fit[,6])[5]
geese.sp.mi <- quantile.mcmc.list(fit[,6])[1]

plot <- data.frame(dabbling.sensitivity = c("dabbling.se.min", "dabbling.se.median", "dabbling.se.95"),
                   dabbling.specificity = c("dabbling.sp.min", "dabblin.sp.median", "dabbling.sp.95"),
                   diving.sensitivity = c("diving.se.min", "diving.se.median", "diving.se.95"), 
                   diving.specificity = c("diving.sp.min", "diving.sp.median", "diving.sp.95"),
                   geese.sensitivity = c("geese.se.min", "geese.se.median", "geese.se.95"), 
                   geese.specificity = c("geese.sp.min", "geese.sp.median", "geese.sp.95"))
plot.se <- data.frame(species = c("1", "2", "3"), 
                      minimum = c("dabbling.se.min", "diving.se.min", "geese.se.min"),
                      median = c("dabbling.se.median", "diving.se.median", "geese.se.median"),
                      ninetyfifth = c("dabbling.se.95", "diving.se.95", "geese.se.95"))
plot.sp <- data.frame(species = c("1", "2", "3"),
                      minimum = c("dabbling.sp.min", "diving.sp.min", "geese.sp.min"),
                      median = c("dabbling.sp.median", "dabbling.sp.median", "geese.sp.median"),
                      ninetyfifth = c("dabbling.sp.95", "dabbling.sp.95", "geese.sp.95"))

ggplot(data = plot.se) +
  geom_points(aes(x = median, y = species)) +
  error_bar(aes(xmin = minimum, xmax = ninetyfifth))

