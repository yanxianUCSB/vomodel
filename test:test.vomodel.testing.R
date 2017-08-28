rm(list = ls())
source('vomodel.testing.R')
source('para.proteinRNA.R')
chipq <- -1
chipw <- 0
system.properties$Chi <- matrix(c(
    0, chipq, 0,0,chipw,
    chipq,0,0,0,0,
    0,0,0,0,0,
    0,0,0,0,0,
    chipw,0,0,0,0
), 5, 5)

# system.properties$polymer.num[1] <- 200
# system.properties$polymer.num[2] <- 3000
# system.properties$sigma[1] <- 0.05
# system.properties$sigma[2] <- 1
# system.properties$MW[1] <- 22e3
# system.properties$MW[2] <- 900e3
system.properties$size.ratio[1:4] <- c(3, 3, 1, 1)
# system.properties$molar.ratio[1:2] <- c(3000, 10)
# fitting.para$sampling.start <- 1e-10
fitting.para$sampling.end <- 0.05
# fitting.para$sampling.gap <- 5e-6
fitting.para$condensation <- T
fitting.para$counterion.release <- T
# fitting.para$binodal.guess <- c(1e-2, 1e-2)
# fitting.para$binodal.guess <- c(1e-3, 1e-4)


# alpha
# print(get.alpha(293, system.properties$size.ratio[1] * k.water.size))

# Binodal curve at 20 C
# ds <- get.binodal.curve(44, sysprop = system.properties, fitting.para = fitting.para,
# condensation = F, counterion.release = F)

# stop()

# phase diagram at 4, 14, 24, 34, 44 C
ds <- get.phase.diagram(system.properties, fitting.para, temp.range = seq(4, 64, 20))

saveRDS(ds, 'out.test.ds.data')
ds <- readRDS('out.test.ds.data')

ds2 <- get.phase.diagram.temp.conc(ds, system.properties, k.conc.salt = 0.030)
ds3 <- get.phase.diagram.temp.nacl(ds, system.properties, k.conc.polymer = 0.125)

phase.diagram.exp <- get.phase.diagram.exp(dataset.file = '~/Box/anywhere/dataset.csv')


# plot(ds2)
# plot(ds3)

library(yxplot)
library(ggplot2)
g <- yxplot.quick(ds$phi.polymer, ds$phi.salt)
g <- yxplot.quick(ds$conc.mass.polymer, ds$conc.salt)
print(g)
