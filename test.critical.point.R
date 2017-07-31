
rm(list = ls())
source('vomodel.R')
source('proteinRNA.para.R')
library(yxplot)
library(ggplot2)
library(dplyr)
library(rootSolve)
library(pracma)
library(nleqslv)
DEBUG <- T
SAVE <- F


# system.properties <- list(
#     polymer.num   = c(207, 900e3/306.2,   1, 1, 1),
#     sigma         = c(11/207, 0.5,   1, 1, 0),
#     size.ratio    = c(   1,    1,   1, 1, 1),
#     molar.ratio   = c(  900e3/306.2,  11, 0.5, 0.5, 0),  # the polycation:polyanion and cation:anion molar ratio
#     water.size    = k.water.size
# )
# fitting.para <- list()
# fitting.para$epsilon <- 1E-8
# fitting.para$sampling.gap <- 1e-4
# fitting.para$critical.point.guess <- c(phi.polymer=0.005, phi.salt=0.005)
# fitting.para$binodal.guess <- c( 0.1,  0.005)  # phi.polymer.2, phi.salt = 0.9 * critical salt



c.point.temp.ds <- c.point.temp(system.properties, fitting.para)
c.point.temp.f <- c.point.temp.fun(c.point.temp.ds)
c.point.temp.f(seq(273.15, 333.15, 0.1))

g <- ggplot(c.point.temp.ds, aes(x = phi.polymer, y = phi.salt, col = temp)) +
    geom_point(size = 1) + 
    labs(x = 'critical polymer fraction',
         y = 'critical salt fraction',
         col = 'temperature [K]') +
    scale_color_continuous(breaks = c(283, 308, 333))
    
g <- theme.background.1(g)
g <- theme.title.text.1(g)
print(g)
if(SAVE) ggsave('test.critical.point.png', width = 5, height = 5)
