#######################################################################
#######################################################################
#Sensitivity and Specificity results (ESA Poster)
#######################################################################
#######################################################################

library(ggplot2)

#read in fit data from selected model
fit <- readRDS("~/Github/short_queens_fit.rds")

dabbling.se.median <- quantile(fit[,1])[3]
dabbling.se.95 <- quantile(fit[,1], probs = 0.95)
dabbling.se.min <- quantile(fit[,1])[1]

dabbling.sp.median <- quantile(fit[,4])[3]
dabbling.sp.95 <- quantile(fit[,4], probs = 0.95)
dabbling.sp.min <- quantile(fit[,4])[1]

diving.se.median <- quantile(fit[,2])[3]
diving.se.95 <- quantile(fit[,2], probs = 0.95)
diving.se.min <- quantile(fit[,2])[1]

diving.sp.median <- quantile(fit[,5])[3]
diving.sp.95 <- quantile(fit[,5], probs = 0.95)
diving.sp.min <- quantile(fit[,5])[1]

geese.se.median <- quantile(fit[,3])[3]
geese.se.95 <- quantile(fit[,3], probs = 0.95)
geese.se.min <- quantile(fit[,3])[1]

geese.sp.median <- quantile(fit[,6])[3]
geese.sp.95 <- quantile(fit[,6], probs = 0.95)
geese.sp.min <- quantile(fit[,6])[1]

plot.se <- data.frame(species = c("1", "2", "3"), 
                      minimum = c(dabbling.se.min, diving.se.min, geese.se.min),
                      median = c(dabbling.se.median, diving.se.median, geese.se.median),
                      ninetyfifth = c(dabbling.se.95, diving.se.95, geese.se.95))
plot.sp <- data.frame(species = c("1", "2", "3"),
                      minimum = c(dabbling.sp.min, diving.sp.min, geese.sp.min),
                      median = c(dabbling.sp.median, dabbling.sp.median, geese.sp.median),
                      ninetyfifth = c(dabbling.sp.95, dabbling.sp.95, geese.sp.95))

jpeg("~/Honors Thesis/Thesis Work/Figures/sensitivityplot.jpeg")
ggplot(data = plot.se) +
  geom_point(aes(x = median, y = species), color = "darkgreen") +
  geom_errorbarh(aes(x = median, y = species, xmin = minimum, xmax = ninetyfifth), height = .1, color = "darkgreen") +
  ylab("Species Group") +
  xlab("Estimated Diagnostic Test Sensitivity") +
  scale_x_continuous(breaks = breaks, limits = c(0, 0.7)) +
  theme(panel.background=element_blank(), axis.ticks = element_line(), axis.line = element_line())
dev.off()

jpeg("~/Honors Thesis/Thesis Work/Figures/specificityplot.jpeg")
ggplot(data = plot.sp) +
  geom_point(aes(x = median, y = species), color = "darkgreen") +
  geom_errorbarh(aes(x = median, y = species, xmin = minimum, xmax = ninetyfifth), height = .1, color = 'darkgreen') +
  ylab("Species Group") +
  xlab("Estimated Diagnostic Test Specificity") +
  scale_x_continuous(breaks = breaks, limits = c(0.9, 1)) +
  theme(panel.background=element_blank(), axis.ticks = element_line(), axis.line = element_line())
dev.off()
