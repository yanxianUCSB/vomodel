# test phi
rm(list = ls())
source('vomodel.R')
phi.polymer <- seq(0.001, 0.4, 0.0001)
phi.salt <- 0.155
temp <- 300
polymer.num <- c(1000, 1000, 1, 1, 1)
alpha <- 3.655
sigma <- c(0.44, 0.44, 1, 1, 0)
size.ratio <- c(1, 1, 1, 1, 1)

g <- sapply(phi.polymer, function(p.polymer) {
  free.energy(phi.polymer = p.polymer,
                 phi.salt = phi.salt,
                 temp = temp,
                 alpha = alpha,
                 sigma = sigma,
                 Chi = 0,
                 polymer.num = polymer.num,
                 size.ratio = size.ratio
                 )
})

dg <- diff(g) / diff(phi.polymer)
ddg <- diff(dg) / diff(phi.polymer)[-1]
dddg <- diff(ddg) / diff(phi.polymer)[-2:-1]
par(mfrow = c(2, 2), mar = rep(2, 4))
plot(phi.polymer, g)
plot(phi.polymer[-1], dg)
plot(phi.polymer[-2:-1], ddg)
plot(phi.polymer[-3:-1], dddg)

# 
# ddg.root <- function(ddg, phi.polymer) {
#   ddg.sp <- splinefun(x = phi.polymer[-2:-1], y = ddg)
# }
