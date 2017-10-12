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
# a <- 2.22
# b <- 5
# c <- 2
# system.properties$size.ratio[1:4] <- c(1.5,  4,  1,  1)
# system.properties$size.ratio[1:4] <- c(a, a*b, c, c)
# system.properties$molar.ratio[1:2] <- c(3000, 10)
# fitting.para$sampling.start <- 1e-6
# fitting.para$sampling.end <- 0.005
# fitting.para$sampling.gap <- 5e-6
fitting.para$condensation <- F
fitting.para$counterion.release <- T
# fitting.para$binodal.guess <- c(1e-4, 1e-4)
fitting.para$binodal.guess <- c(1e-3, 5e-5)


# alpha
# print(get.alpha(293, system.properties$size.ratio[1] * k.water.size))

# Binodal curve at 20 C
# ds <- get.binodal.curve(44, sysprop = system.properties, fitting.para = fitting.para,
# condensation = F, counterion.release = F)

# stop()

# phase diagram at 4, 14, 24, 34, 44 C
# ds <- bind_rows(lapply(seq(1, 1, 10), function(size.ratio12){
#     bind_rows(lapply(seq(10, 10, 30), function(tempC){
#         bind_rows(lapply(seq(0, 0, 1), function(chipq){
#             bind_rows(lapply(seq(-0, -0, 0.1), function(chipw){
#                 
#                 chipp <- 0
#                 chipq <- 0
#                 chipw <- 0
#                 chiqw <- 0
#                 # system.properties$size.ratio[1:4] <- c(size.ratio12, size.ratio12, 2, 2)
#                 Chi <-  matrix(rep(0, 25), 5, 5)
#                 Chi[1,2] <- chipq
#                 Chi[1,5] <- chipw
#                 Chi[2,5] <- chiqw
#                 Chi[1,1] <- chipp
#                 # Chi[1,5] <- 0
#                 Chi <- Chi + t(Chi)
#                 system.properties$Chi <- Chi
#                 d <- get.phase.diagram(system.properties, fitting.para, temp.range = tempC) 
#                 if(!is.null(d))
#                     d <- d %>% mutate(chipq = chipq,
#                                       chipw = chipw,
#                                       size.ratio12 = size.ratio12)
#             }))
#         }))
#     }))
# }))

# saveRDS(ds, 'out.test.ds.data')
# ds <- readRDS('out.test.ds.data')

ds <- get.phase.diagram(system.properties, fitting.para, seq(0, 200, 100))
# ds3 <- get.phase.diagram.temp.nacl(ds, system.properties, k.conc.polymer = 0.125)

phase.diagram.exp <- get.phase.diagram.exp(dataset.file = '~/Box/anywhere/dataset.csv')


# plot(ds2)
# plot(ds3)

library(yxplot)
library(ggplot2)
g <- ggplot(ds, aes(x = phi.polymer, y = (phi.salt), group = tempC)) +
    geom_label(data = ds %>%
                   group_by(tempC) %>%
                   mutate(max.phi = phi.polymer[which.min(abs(phi.polymer - quantile(phi.polymer, 0.9)))]) %>%
                   mutate(max.phis = mean(phi.salt)) %>%
                   ungroup() %>%
                   mutate(note = paste0(
                       'alpha = ', round(get.alpha(tempC + 273.15, k.water.size), 2), '\n',
                       'sigma = ', round(sigma.q, 2)
                   )),
               aes(x = max.phi,
                   y = max.phis,
                   label = note,
                   col = tempC),
               show.legend = F) +
    geom_point(aes(col = tempC)) +
    scale_color_continuous(guide = 'legend', breaks = unique(ds$tempC)) +
    labs(x = 'Total Polymer Frac. [%]',
         y = 'Salt Frac. [%]',
         col = 'Temp. [ºC]')
g <- theme.background.1(g)
g <- theme.title.text.1(g)
print(g)

ggsave('~/Desktop/ParaProteinRNA_Commit_8253548_cirel.png', width = 6, height = 4)


ds2 <- get.phase.diagram.temp.conc(ds, system.properties, k.phi.salt = 0.0028)
g.temp.phi.polymer <- ggplot(ds2, aes(x = phi.polymer, y = tempC)) +
    geom_point( size = 2) +
    # geom_line(aes(col = 'sim'), lwd = 1.5) +
    labs(x = 'Total Polymer Frac. [%]',
         y = 'Salt Frac. [%]',
         col = '')
print(g.temp.phi.polymer)

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
