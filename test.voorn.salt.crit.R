# test phi
rm(list = ls())
source('vomodel.R')
phi.polymer <- seq(0.001, 0.05, 0.001)
phi.salt <- 0.15
temp <- 300
polymer.num <- c(1000, 1000, 1, 1, 1)
alpha <- 3.655
sigma <- c(0.44, 0.44, 1, 1, 0)
size.ratio <- c(1, 1, 1, 1, 1)

phi.salt.crit.get <-
  function(phi.polymer,
           phi.salt = seq(1E-6, 0.5, 1E-3),
           ...) {
    sapply(phi.polymer, function(phi.polymer) {
      g <- sapply(phi.salt, function(phi.salt) {
        free.energy(
          phi.polymer = phi.polymer,
          phi.salt = phi.salt,
          temp = temp,
          alpha = alpha,
          sigma = sigma,
          Chi = 0,
          polymer.num = polymer.num,
          size.ratio = size.ratio
        )
      })
      dg <- diff(g) / diff(phi.salt)
      ddg <- diff(dg) / diff(phi.salt)[-1]
      dddg <- diff(ddg) / diff(phi.salt)[-2:-1]
      dddg.sp.fun <- function(x) {
        spline(x = phi.salt[-3:-1],
               y = dddg,
               xout = x)$y
      }
      root <-
        uniroot(f = dddg.sp.fun, interval = range(phi.salt[-3:-1]))$root
      return(root)
    })
  }

g <- yxplot.quick(phi.polymer, phi.salt.crit.get(phi.polymer)) +
  labs(x = 'Polymer Conc [% v.v]',
       y = 'Critical NaCl Conc [% v.v]') 
ggsave('test.voorn.salt.crit.png', width = 5, height = 5)

