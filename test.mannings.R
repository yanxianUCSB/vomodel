# manning's condensation
rm(list = ls())
source('vomodel.R')
system.properties <- list(
    polymer.num   = c(207, 900e3/306.2,   1, 1, 1),
    sigma         = c(11/207, 0.5,   1, 1, 0),
    size.ratio    = c(   1,    1,   1, 1, 1),
    molar.ratio   = c(  900e3/306.2,  11, 0.5, 0.5, 0),  # the polycation:polyanion and cation:anion molar ratio
    water.size    = k.water.size
)
fitting.para <- list()
fitting.para$epsilon <- 1E-8
fitting.para$sampling.gap <- 1e-4
fitting.para$critical.point.guess <- c(phi.polymer=0.005, phi.salt=0.01)
fitting.para$binodal.guess <- c( 0.1,  0.005)  # phi.polymer.2, phi.salt = 0.9 * critical salt
