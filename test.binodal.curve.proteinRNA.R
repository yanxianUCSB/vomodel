# Binodal curve at different temperature
rm(list = ls())
source('vomodel.R')
library(yxplot)
library(ggplot2)
library(dplyr)
library(rootSolve)
library(pracma)
library(nleqslv)
DEBUG <- T
SAVE <- F

get.binodal.curve.proteinRNA <- function() {
    
    system.properties <- list(
        polymer.num   = c(207, 900e3 / 306.2,   1, 1, 1),
        sigma         = c(11 / 207, 0.5,   1, 1, 0),
        size.ratio    = c(k.water.size * 2.1, k.dna.contour.unit.length, k.na.size, k.cl.size, k.water.size) / k.water.size,
        MW            = c(22e3, 900e3),
        charge.ratio  = c(1, 1),
        molar.ratio   = c(900e3 / 306.2,  11, 0.5, 0.5, 0),
        # the polycation:polyanion and cation:anion molar ratio
        water.size    = k.water.size,
        Chi = matrix(c(.0,.0,0,0,0,
                       .0,.0,0,0,0,
                       0,0,0,0,0,
                       0,0,0,0,0,
                       0,0,0,0,0), 5,5)
    )
    fitting.para <- list(
        epsilon = 1E-8 , 
        sampling.start = 6e-5,
        sampling.gap = 6e-8 ,
        critical.point.guess = c(phi.polymer = 0.004, phi.salt = 0.00035) ,
        c.point.temp.fun = c.point.temp.fun(c.point.temp(system.properties, fitting.para)) ,
        binodal.guess = 0.021  # phi.polymer.2
    )
    sampling <- list(
        tempC = 20
    )
    
    p <- lapply(sampling$tempC, function(tempC) {
        # update critical.point.guess
        fitting.para$critical.point.guess <- as.numeric(fitting.para$c.point.temp.fun(tempC + 273))
        get.binodal.curve(tempC, Chi = system.properties$Chi, system.properties, fitting.para, unit = 'mol')
    })
    if (DEBUG) {
        print(head(p[[1]]))
        par(mfrow = c(1, 2))
        plot(p[[1]]$phi.polymer, p[[1]]$phi.salt)
        plot(p[[1]]$conc.p, p[[1]]$conc.salt)
    }
    
    return(do.call(rbind, p))
}

# for(to.test in seq(2.1, 6, 0.2)) tryCatch(
p <- get.binodal.curve.proteinRNA()
# )
if (SAVE)
    saveRDS(p, 'binodal.curve.4C.40C.Chi0.dataset', ascii = T)
# test Chi

# test assymmetric