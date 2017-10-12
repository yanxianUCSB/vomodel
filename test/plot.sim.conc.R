rm(list = ls())
source('vomodel.testing.R')
source('para.proteinRNA.R')
library(dplyr)
library(nleqslv)
library(ggplot2)
library(yxplot)

DEBUG <<- T
SAVE <<- T

k.conc.salt       <<- 0.030
k.conc.polymer    <<-  5E-6 * system.properties$MW[1] + 15E-3
# k.sim.temp.range  <<- seq(4, 44, 10)

# Water size
# watersize <<- 3.1E-10
# k.amino.acid.length       <<- watersize * 1.0
# k.dna.contour.unit.length <<- watersize * 1.0
# k.na.size                 <<- watersize * 1.1
# k.cl.size                 <<- watersize * 1.1
# k.water.size              <<- watersize
# system.properties$size.ratio <-  c(k.amino.acid.length, k.dna.contour.unit.length, k.na.size, k.cl.size, k.water.size) / k.water.size

# a <- 2.22
# b <- 5
# c <- 2
system.properties$size.ratio[1:4] <- c(2.818831, 13.912664,  1.288291,  1.288291)
# system.properties$size.ratio[1:4] <- c(a, a*b, c, c)

# Chi
Chi <-  matrix(rep(0, 25), 5, 5)
Chi[1,2] <- -0
Chi[1,5] <- 0

# Chi[1,5] <- 0
Chi <- Chi + t(Chi)
system.properties$Chi <- Chi

fitting.para$counterion.release <- F



# # #
fitting.para$sampling.start <- 1e-10
fitting.para$sampling.gap <- 1e-5
fitting.para$binodal.guess <- c(1e-3, 5e-3)

phase.diagram.exp <- get.phase.diagram.exp(dataset.file = '~/Box/anywhere/dataset.csv')

ds <- get.phase.diagram(system.properties, fitting.para, temp.range = seq(-160, -120, 5))
saveRDS(ds, 'out.test.ds.data')
ds <- readRDS('out.test.ds.data')

ds2 <- get.phase.diagram.temp.conc(ds, system.properties, k.conc.salt = k.conc.salt)
ds3 <- get.phase.diagram.temp.nacl(ds, system.properties, k.conc.polymer = k.conc.polymer)

assertthat::assert_that(!is.null(ds))

g <- ggplot(ds, aes(x = conc.polymer, y = conc.salt, group = tempC)) +
    geom_point(aes(col = tempC)) +
    geom_line(aes(x = conc.polymer, y = conc.salt, col = tempC)) +
    scale_color_continuous(guide = 'legend') +
    labs(x = 'Conc.polymer[mg/mL]',
         y= 'NaCl [M]', 
         col = 'Temp. [C]')
    

g2.transform <- function(x){-log(50-x)}
g2 <- ggplot() + 
    geom_point(data = ds2 %>% filter(conc.polymer < 20), 
               aes(x = conc.polymer, y = g2.transform(tempC), pch = 'model'), 
               size = 3) +
    geom_point(data = phase.diagram.exp %>% filter(conc.salt == 0.03),
               aes(x = conc.polymer, y = g2.transform(tempC.cpm), pch = 'experiment'),
               size = 3) +
    geom_errorbar(aes(x = conc.polymer,
                      ymax = g2.transform(tempC.cpm + tempC.cpsd),
                      ymin = g2.transform(tempC.cpm - tempC.cpsd)
                      ), data = phase.diagram.exp %>% filter(conc.salt == 0.03),
                  lwd = 1) +
    # geom_point(data = phase.diagram.exp %>% filter(conc.salt == 0.03), aes(x = conc.polymer, y = tempC.on, col = 'exp.on')) +
    scale_y_continuous(limits = g2.transform(c(-160, 40)),
                     breaks = g2.transform(seq(-160, 40, 40)), 
                     labels = seq(-160, 40, 40)) +
    labs(x = 'Conc.polymer[mg/mL]', 
         y = 'Temperature [C]', 
         pch = '') 

g3.transform <- function(x){-log(50-x)}
g3 <- ggplot() + 
    geom_point(data = ds3, 
               aes(x = conc.salt, y = g3.transform(tempC), pch = 'model'), 
               size = 3) +
    geom_point(data = phase.diagram.exp %>% filter(abs(conc.polymer - 0.125) < 1e-3),
               aes(x = conc.salt, y = g3.transform(tempC.cpm), pch = 'experiment'),
               size = 3) +
    geom_errorbar(aes(x = conc.salt,
                      ymax = g2.transform(tempC.cpm + tempC.cpsd),
                      ymin = g2.transform(tempC.cpm - tempC.cpsd)
                      ), 
                  data = phase.diagram.exp %>% filter(abs(conc.polymer - 0.125) < 1e-3)) +
    scale_y_continuous(limits = g3.transform(c(-160, 40)),
                     breaks = g3.transform(seq(-160, 40, 40)), 
                     labels = seq(-160, 40, 40)) +
    # geom_point(data = phase.diagram.exp %>% filter(abs(conc.polymer - 0.125) < 1e-3), aes(x = conc.salt, y = tempC.on, col = 'exp.on')) +
    labs(x = 'NaCl [M]', 
         y = 'Temperature [C]', 
         pch = '') 
 
g <- theme.background.1(g)
g <- theme.title.text.1(g)
g2 <- theme.background.1(g2)
g2 <- theme.title.text.1(g2)
g3 <- theme.background.1(g3)
g3 <- theme.title.text.1(g3)
if (!SAVE) {
multiplot(g, g2, g3, cols = 3)
}

if(SAVE) {
    ggsave('plot.sim.conc.bncurve.png', g, width = 5, height = 5)
    ggsave( 'plot.sim.conc.conc.png', g2, width = 5, height = 5)
    ggsave( 'plot.sim.conc.nacl.png', g3, width = 5, height = 5)
    ggsave('plot.sim.conc.bncurve.pdf', g, width = 5, height = 5)
    ggsave( 'plot.sim.conc.conc.pdf', g2, width = 5, height = 5)
    ggsave( 'plot.sim.conc.nacl.pdf', g3, width = 5, height = 5)
    saveRDS(object = ds, file = 'plot.sim.conc.ds', ascii = T)
    saveRDS(object = ds2, file = 'plot.sim.conc.conc.ds', ascii = T)
    saveRDS(object = ds3, file = 'plot.sim.conc.nacl.ds', ascii = T)
    saveRDS(system.properties, 'plot.sim.conc.sysprop', ascii = T)
    saveRDS(fitting.para, 'plot.sim.conc.fitpara', ascii = T)
}

