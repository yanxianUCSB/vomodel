# test.freeEnergy.grad free energy
# using pracama::grad() to generate derivative

rm(list = ls())
source('vomodel.R')
phi.polymer <- seq(0.01, 0.4, 0.0001)
phi.salt <- 0.155
temp <- 300
polymer.num <- c(1000, 1000, 1, 1, 1)
alpha <- 3.655
sigma <- c(0.44, 0.44, 1, 1, 0)
size.ratio <- c(1, 1, 1, 1, 1)
Chi <- 0

ds <- gibbs.funs(phi.polymer, phi.salt = phi.salt, temp = temp, polymer.num = polymer.num,
                       alpha = alpha, sigma = sigma, size.ratio = size.ratio, Chi = Chi)
plot(ds[-2])
head(ds)
