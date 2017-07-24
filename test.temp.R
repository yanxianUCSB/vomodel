# test.numericDeriv free energy

rm(list = ls())
source('vomodel.R')
phi.polymer <- seq(0.01, 0.4, 0.0001)
phi.salt <- 0.01
temp <- 300
polymer.num <- c(1000, 1000, 1, 1, 1)
alpha <- 3.655
sigma <- c(0.44, 0.44, 1, 1, 0)
size.ratio <- c(1, 1, 1, 1, 1)
Chi <- 0

ds <- free.energy.funs(phi.polymer, phi.salt = 0.2, temp = temp, polymer.num = polymer.num,
alpha = alpha, sigma = sigma, size.ratio = size.ratio, Chi = Chi)
# ds <- critical.point()
plot(ds)

ds <- sapply(phi.polymer, function(x) {
  
  critical.point.fun(
    c(x, 0.12),
    alpha = alpha,
    sigma = sigma,
    size.ratio = size.ratio,
    Chi = Chi,
    polymer.num = polymer.num,
    temp = temp
  )[[1]]
})

plot(phi.polymer, ds)

ds <- critical.point(alpha = alpha,
               sigma = sigma,
               size.ratio = size.ratio,
               Chi = Chi,
               polymer.num = polymer.num,
               temp = temp)
print(ds)
# critical.point(alpha, sigma, Chi)
