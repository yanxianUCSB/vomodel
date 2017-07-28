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
        size.ratio    = c(1,    1,   1, 1, 1),
        molar.ratio   = c(900e3 / 306.2,  11, 0.5, 0.5, 0),
        # the polycation:polyanion and cation:anion molar ratio
        water.size    = k.water.size,
        Chi = matrix(c(0,0,0,0,0,
                       0,0,0,0,0,
                       0,0,0,0,0,
                       0,0,0,0,0,
                       0,0,0,0,0), 5,5)
    )
    fitting.para <- list(
        epsilon = 1E-8 , 
        sampling.gap = 1e-5 ,
        critical.point.guess = c(phi.polymer = 0.005, phi.salt = 0.005) ,
        c.point.temp.fun = c.point.temp.fun(c.point.temp(system.properties, fitting.para)) ,
        binodal.guess = 0.1  # phi.polymer.2
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
        plot(p[[1]]$phi.polymer, p[[1]]$phi.salt)
    }
    
    return(do.call(rbind, p))
}

p <- get.binodal.curve.proteinRNA()
if (SAVE)
    saveRDS(p, 'binodal.curve.4C.60C.Chi0.dataset', ascii = T)
# test Chi

# test assymmetric