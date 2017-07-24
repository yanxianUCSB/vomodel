# test phi
rm(list = ls())
source('vomodel.R')
source('test.peek.para.R')
phi.polymer <- phi.polymer.seq
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
