# test binodal curve
rm(list = ls())
source('vomodel.R')
library(yxplot)
library(ggplot2)
library(rootSolve)

phi.polymer <- seq(0.0001, 0.4, 0.001)
phi.salt <- seq(0.0001, 0.2, 0.01)
temp <- 300
polymer.num <- c(1000, 1000, 1, 1, 1)
alpha <- 3.655
sigma <- c(0.44, 0.44, 1, 1, 0)
size.ratio <- c(1, 1, 1, 1, 1)

y <-  free.energy.funs(phi.polymer, phi.salt = 0.012,
                                                       temp = temp,
                                                       alpha = alpha,
                                                       sigma = sigma,
                                                       Chi = 0,
                                                       polymer.num = polymer.num,
                                                       size.ratio = size.ratio)

p <- binodal.curve(phi.polymer.seq = phi.polymer)
print(p)

plot(y$phi, y$f)
abline(a = p$par[1], b = p$par[2], col = 'red')
