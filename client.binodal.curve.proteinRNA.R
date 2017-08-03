# Binodal curve at different temperature
rm(list = ls())
source('vomodel.R')
library(yxplot)
library(ggplot2)
library(dplyr)
library(rootSolve)
library(pracma)
library(nleqslv)
DEBUG <- F
SAVE <- T

get.binodal.curve.proteinRNA <- function() {
    
    system.properties <- list(
        polymer.num   = c(207, 900e3 / 306.2,   1, 1, 1),
        sigma         = c(11 / 207, 0.5,   1, 1, 0),
        size.ratio    = c(k.amino.acid.length, k.dna.contour.unit.length, k.na.size, k.cl.size, k.water.size) / k.water.size,
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
        sampling.start = 1e-5,
        sampling.gap = 1e-7 ,
        critical.point.guess = c(phi.polymer = 0.01, phi.salt = 0.0001) ,
        c.point.temp.fun = c.point.temp.fun(c.point.temp(system.properties, fitting.para)) ,
        binodal.guess = 0.1  # phi.polymer.2
    )
    sampling <- list(
        tempC = seq(4, 60, 0.1)
    )
    
    p <- lapply(sampling$tempC, function(tempC) {
        # update Bjerrum length and the effective charge density of RNA
        lB <- ke^2 / (kEr*kkB*(tempC+273.15))
        system.properties$sigma[2] <- system.properties$size.ratio[2]*k.water.size / lB
        # update critical.point.guess
        fitting.para$critical.point.guess <- as.numeric(fitting.para$c.point.temp.fun(tempC + 273))
        cat(paste0('Temp [C]: ', tempC, '\n'))
        cat(paste0('Bjerrum length [m]: ', lB, '\n'))
        cat(paste0('Sigma RNA: ', system.properties$sigma[2], '\n'))
        cat('fitting >>>\n')
        out <- get.binodal.curve(tempC, Chi = system.properties$Chi, system.properties, fitting.para, unit = 'mol')
        cat('succeeded!\n')
        return(out)
    })
    
    if (DEBUG) {
        print(head(p[[1]]))
        par(mfrow = c(1, 2))
        plot(p[[1]]$phi.polymer, p[[1]]$phi.salt)
        plot(p[[1]]$conc.p, p[[1]]$conc.salt)
    }
    
    return(do.call(rbind, p))
}

p <- get.binodal.curve.proteinRNA()

if (SAVE) {
    cat('saving data >>>\n')
    saveRDS(p, 'binodal.curve.4C.60C.Chi0.dataset', ascii = T)
    cat('saved!')
}

