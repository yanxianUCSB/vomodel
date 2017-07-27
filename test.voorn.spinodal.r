# test voorn.salt  //todo:
rm(list = ls())
source('vomodel.R')
library(yxplot)
library(yxhelper)
library(ggplot2)
library(rootSolve)

phi.polymer <- c(seq(0.000001, 0.03, 0.0001), seq(0.03, 0.4, 0.01))
phi.salt <- seq(0.0001, 0.2, 0.001)
temp <- 300
polymer.num <- c(1000, 1000, 1, 1, 1)
alpha <- 3.655
sigma <- c(0.44, 0.44, 1, 1, 0)
size.ratio <- c(1, 1, 1, 1, 1)


ds <- spinodal.curve(phi.polymer, phi.salt, x.axis = 'phi.polymer', 
                     temp, alpha, sigma, Chi = 0, polymer.num, size.ratio, 
                  curve.type = 'spinodal')

g <- yxplot.quick(ds$root, ds$para) +
  labs(x = 'Polymer Conc [% v.v]',
       y = 'NaCl Conc [% v.v]')
ggsave('test.voorn.spinodal.png', width = 5, height = 5)

