# test cloudpoint concentration 
rm(list = ls())
source('vomodel.R')
library(yxplot)
library(ggplot2)
library(dplyr)
library(rootSolve)
library(pracma)
DEBUG <- F

system.properties <- list(
    polymer.num = c(1000, 1000, 1, 1, 1),
    sigma = c(0.30, 0.30, 1, 1, 0),
    size.ratio = c(1, 1, 1, 1, 1),
    water.size = k.water.size
)
fitting.para <- list()
fitting.para$epsilon <- 1E-8
fitting.para$sampling.gap <- 1e-4
fitting.para$critical.point.guess <- c(phi.polymer=0.01, phi.salt=0.15)

tempCs <- 4

fitting.para$binodal.guess <- c(0.03, 0.008)
p2 <- get.binodal.curves(tempCs, 0, system.properties, fitting.para)

head(p2)

g <- ggplot(p2, aes(y = conc.salt, group = tempC)) +
    geom_line(aes(x = conc.p, col = tempC), lwd = 2) +
    scale_colour_continuous(guide = 'colorbar') +
    labs(x = 'Polymer [mol/L]',
         y = 'Salt [mol/L]',
         col = expression('Temperature'~degree~C))
g <- theme.title.text.1(g)
print(g)

