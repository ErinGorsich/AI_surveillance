#######################################################################
#######################################################################
#Sensitivity and Specificity results (ESA Poster)
#######################################################################
#######################################################################

library(ggplot2)

#read in fit data from selected model
fit.queens <- readRDS("~/Github/short_queens_fit.rds")
fit.base <- readRDS("~/Github/short_base_fit.rds")
fit.network <- readRDS("~/Github/short_network_fit.rds")

#queens quantile calculations
qdabbling.se.median <- quantile(fit.queens[,1])[3]
qdabbling.se.95 <- quantile(fit.queens[,1], probs = 0.95)
qdabbling.se.min <- quantile(fit.queens[,1])[1]

qdabbling.sp.median <- quantile(fit.queens[,4])[3]
qdabbling.sp.95 <- quantile(fit.queens[,4], probs = 0.95)
qdabbling.sp.min <- quantile(fit.queens[,4])[1]

qdiving.se.median <- quantile(fit.queens[,2])[3]
qdiving.se.95 <- quantile(fit.queens[,2], probs = 0.95)
qdiving.se.min <- quantile(fit.queens[,2])[1]

qdiving.sp.median <- quantile(fit.queens[,5])[3]
qdiving.sp.95 <- quantile(fit.queens[,5], probs = 0.95)
qdiving.sp.min <- quantile(fit.queens[,5])[1]

qgeese.se.median <- quantile(fit.queens[,3])[3]
qgeese.se.95 <- quantile(fit.queens[,3], probs = 0.95)
qgeese.se.min <- quantile(fit.queens[,3])[1]

qgeese.sp.median <- quantile(fit.queens[,6])[3]
qgeese.sp.95 <- quantile(fit.queens[,6], probs = 0.95)
qgeese.sp.min <- quantile(fit.queens[,6])[1]

#base quantile calculations
bdabbling.se.median <- quantile(fit.base[,1])[3]
bdabbling.se.95 <- quantile(fit.base[,1], probs = 0.95)
bdabbling.se.min <- quantile(fit.base[,1])[1]

bdabbling.sp.median <- quantile(fit.base[,4])[3]
bdabbling.sp.95 <- quantile(fit.base[,4], probs = 0.95)
bdabbling.sp.min <- quantile(fit.base[,4])[1]

bdiving.se.median <- quantile(fit.base[,2])[3]
bdiving.se.95 <- quantile(fit.base[,2], probs = 0.95)
bdiving.se.min <- quantile(fit.base[,2])[1]

bdiving.sp.median <- quantile(fit.base[,5])[3]
bdiving.sp.95 <- quantile(fit.base[,5], probs = 0.95)
bdiving.sp.min <- quantile(fit.base[,5])[1]

bgeese.se.median <- quantile(fit.base[,3])[3]
bgeese.se.95 <- quantile(fit.base[,3], probs = 0.95)
bgeese.se.min <- quantile(fit.base[,3])[1]

bgeese.sp.median <- quantile(fit.base[,6])[3]
bgeese.sp.95 <- quantile(fit.base[,6], probs = 0.95)
bgeese.sp.min <- quantile(fit.base[,6])[1]

#network quantile calculations
ndabbling.se.median <- quantile(fit.network[,1])[3]
ndabbling.se.95 <- quantile(fit.network[,1], probs = 0.95)
ndabbling.se.min <- quantile(fit.network[,1])[1]

ndabbling.sp.median <- quantile(fit.network[,4])[3]
ndabbling.sp.95 <- quantile(fit.network[,4], probs = 0.95)
ndabbling.sp.min <- quantile(fit.network[,4])[1]

ndiving.se.median <- quantile(fit.network[,2])[3]
ndiving.se.95 <- quantile(fit.network[,2], probs = 0.95)
ndiving.se.min <- quantile(fit.network[,2])[1]

ndiving.sp.median <- quantile(fit.network[,5])[3]
ndiving.sp.95 <- quantile(fit.network[,5], probs = 0.95)
ndiving.sp.min <- quantile(fit.network[,5])[1]

ngeese.se.median <- quantile(fit.network[,3])[3]
ngeese.se.95 <- quantile(fit.network[,3], probs = 0.95)
ngeese.se.min <- quantile(fit.network[,3])[1]

ngeese.sp.median <- quantile(fit.network[,6])[3]
ngeese.sp.95 <- quantile(fit.network[,6], probs = 0.95)
ngeese.sp.min <- quantile(fit.network[,6])[1]

plot.se <- data.frame(species = seq(1,3), 
                      minimum.base = c(bdabbling.se.min, bdiving.se.min, bgeese.se.min),
                      minimum.queens = c(qdabbling.se.min, qdiving.se.min, qgeese.se.min),
                      minimum.network = c(ndabbling.se.min, ndiving.se.min, ngeese.se.min),
                      median.base = c(bdabbling.se.median, bdiving.se.median, bgeese.se.median),
                      median.queens = c(qdabbling.se.median, qdiving.se.median, qgeese.se.median),
                      median.network = c(ndabbling.se.median, ndiving.se.median, ngeese.se.median),
                      ninetyfifth.base = c(bdabbling.se.95, bdiving.se.95, bgeese.se.95),
                      ninetyfifth.queens = c(qdabbling.se.95, qdiving.se.95, qgeese.se.95),
                      ninetyfifth.network = c(ndabbling.se.95, ndiving.se.95, ngeese.se.95))
plot.sp <- data.frame(species = seq(1,3),
                      minimum.base = c(bdabbling.sp.min, bdiving.sp.min, bgeese.sp.min),
                      minimum.queens = c(qdabbling.sp.min, qdiving.sp.min, qgeese.sp.min),
                      minimum.network = c(ndabbling.sp.min, ndiving.sp.min, ngeese.sp.min),
                      median.base = c(bdabbling.sp.median, bdabbling.sp.median, bgeese.sp.median),
                      median.queens = c(qdabbling.sp.median, qdabbling.sp.median, qgeese.sp.median),
                      median.network = c(ndabbling.sp.median, ndabbling.sp.median, ngeese.sp.median),
                      ninetyfifth.base = c(bdabbling.sp.95, bdabbling.sp.95, bgeese.sp.95),
                      ninetyfifth.queens = c(qdabbling.sp.95, qdabbling.sp.95, qgeese.sp.95),
                      ninetyfifth.network  = c(ndabbling.sp.95, ndabbling.sp.95, ngeese.sp.95))

breaks <- seq(0,1, 0.1)

setwd("~/HP/Figures")
jpeg("sensitivityplot.jpeg")
ggplot(data = plot.se) +
  geom_point(aes(x = median.base, y = species), color = "darkblue") +
  geom_errorbarh(aes(x = median.base, y = species, xmin = minimum.base, xmax = ninetyfifth.base), 
                 height = .1, color = "darkblue") +
  geom_point(aes(x = median.queens, y = species), color = "skyblue") +
  geom_errorbarh(aes(x = median.queens, y = species, xmin = minimum.queens, xmax = ninetyfifth.queens),
                 height = .1, color = "skyblue") +
  geom_point(aes(x = median.network, y = species), color = "cornflowerblue") +
  geom_errorbarh(aes(x = median.network, y = species, xmin = minimum.network, xmax = ninetyfifth.network),
                 height = .1, color = "cornflowerblue") +
  ylab("Species Group") +
  xlab("Estimated Diagnostic Test Sensitivity") +
  scale_x_continuous(breaks = breaks, limits = c(0, 0.7)) +
  scale_y_continuous(breaks = c(1, 2, 3), labels = c("1", "2", "3")) +
  theme(panel.background=element_blank(), axis.ticks = element_line(), axis.line = element_line())
dev.off()

jpeg("specificityplot.jpeg")
ggplot(data = plot.sp) +
  geom_point(aes(x = median.base, y = species), color = "darkblue") +
  geom_errorbarh(aes(x = median.base, y = species, xmin = minimum.base, xmax = ninetyfifth.base), 
                 height = .1, color = 'darkblue') +
  geom_point(aes(x = median.queens, y = species), color = "skyblue") +
  geom_errorbarh(aes(x = median.queens, y = species, xmin = minimum.queens, xmax = ninetyfifth.queens),
                 height = .1, color = "skyblue") +
  geom_point(aes(x = median.network, y = species), color = "cornflowerblue") +
  geom_errorbarh(aes(x = median.network, y = species, xmin = minimum.network, xmax = ninetyfifth.network),
                 height = .1, color = "cornflowerblue") +
  ylab("Species Group") +
  xlab("Estimated Diagnostic Test Specificity") +
  scale_x_continuous(breaks = breaks, limits = c(0.9, 1)) +
  scale_y_continuous(breaks = c(1, 2, 3), labels = c("1", "2", "3")) +
  theme(panel.background=element_blank(), axis.ticks = element_line(), axis.line = element_line())
dev.off()
