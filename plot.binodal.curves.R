rm(list = ls())
source('vomodel.testing.R')
source('para.proteinRNA.R')
library(nleqslv)
library(dplyr)
library(numDeriv)

# system.properties$polymer.num[1] <- 200
# system.properties$polymer.num[2] <- 3000
# system.properties$sigma[1] <- 0.05
# system.properties$sigma[2] <- 1
# system.properties$MW[1] <- 22e3
# system.properties$MW[2] <- 900e3
# system.properties$size.ratio[1:4] <- c(6, 6, 1, 1)
# system.properties$molar.ratio[1:2] <- c(3000, 10)
# fitting.para$sampling.start <- 1e-10
# fitting.para$sampling.end <- 0.1
# fitting.para$sampling.gap <- 5e-6
fitting.para$condensation <- F
fitting.para$counterion.release <- T
# fitting.para$binodal.guess <- c(1e-4, 1e-4)
fitting.para$binodal.guess <- c(1e-4, 1e-5)


# alpha
# print(get.alpha(293, system.properties$size.ratio[1] * k.water.size))

# Binodal curve at 20 C
# ds <- get.binodal.curve(44, sysprop = system.properties, fitting.para = fitting.para,
# condensation = F, counterion.release = F)

# stop()

# phase diagram at 4, 14, 24, 34, 44 C
ds <- bind_rows(lapply(seq(1, 1, 10), function(size.ratio12){
    bind_rows(lapply(seq(0, 60, 30), function(tempC){
        bind_rows(lapply(seq(0, 0, 1), function(chipq){
            bind_rows(lapply(seq(-0, -0, 0.1), function(chipw){
                
                # system.properties$size.ratio[1:4] <- c(size.ratio12, size.ratio12, 2, 2)
                system.properties$Chi <- matrix(c(
                    0, chipq, 0,0,chipw,
                    chipq,0,0,0,0,
                    0,0,0,0,0,
                    0,0,0,0,0,
                    chipw,0,0,0,0
                ), 5, 5)
                d <- get.phase.diagram(system.properties, fitting.para, temp.range = tempC) 
                if(!is.null(d))
                    d <- d %>% mutate(chipq = chipq,
                                      chipw = chipw,
                                      size.ratio12 = size.ratio12)
            }))
        }))
    }))
}))

# saveRDS(ds, 'out.test.ds.data')
# ds <- readRDS('out.test.ds.data')

ds2 <- get.phase.diagram.temp.conc(ds, system.properties, k.conc.salt = 0.030)
ds3 <- get.phase.diagram.temp.nacl(ds, system.properties, k.conc.polymer = 0.125)

phase.diagram.exp <- get.phase.diagram.exp(dataset.file = '~/Box/anywhere/dataset.csv')


# plot(ds2)
# plot(ds3)

library(yxplot)
library(ggplot2)
g <- ggplot(ds, aes(x = conc.polymer, y = (conc.salt))) +
    geom_point(aes(col = tempC)) +
    scale_color_continuous(guide = 'legend', breaks = unique(ds$tempC)) +
    # scale_y_log10() +
    facet_wrap(size.ratio12 + chipw~chipq, scales = 'free')
    
# g <- yxplot.quick(ds$phi.polymer, ds$phi.salt)
# g <- yxplot.quick(ds$conc.mass.polymer, ds$conc.salt)
print(g)

stop()

g2 <- ggplot(ds2, aes(x = conc.polymer, y = tempC)) +
    geom_point(aes(col = 'sim'), size = 2) +
    geom_line(aes(col = 'sim'), lwd = 1.5) +
    geom_point(data = phase.diagram.exp %>% filter(abs(conc.salt-0.03) < 1e-2), 
               aes(x = conc.polymer, y = tempC.cpm, col = 'exp.cp')) +
    geom_point(data = phase.diagram.exp %>% filter(abs(conc.salt-0.03) < 1e-2), 
               aes(x = conc.polymer, y = tempC.onm, col = 'exp.on')) +
    labs(x = 'Polymer [mg/mL]',
         y = 'Temp. [C]')
g3 <- ggplot(ds3, aes(x = conc.salt, y = tempC)) +
    geom_point(aes(col = 'sim'), size = 2) +
    geom_line(aes(col = 'sim'), lwd = 1.5) +
    geom_point(data = phase.diagram.exp %>% filter(abs(conc.polymer-0.125) < 1e-2), 
               aes(x = conc.salt, y = tempC.cpm, col = 'exp.cp')) +
    geom_point(data = phase.diagram.exp %>% filter(abs(conc.polymer-0.125) < 1e-2), 
               aes(x = conc.salt, y = tempC.onm, col = 'exp.on')) +
    labs(x = 'NaCl [M]',
         y = 'Temp. [C]')
g3
g2
