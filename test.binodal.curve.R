# Binodal curve at different temperature
rm(list = ls())
source('vomodel.R')
source('para.voorn.R')
library(yxplot)
library(yxhelper)
library(ggplot2)
library(dplyr)
library(rootSolve)
library(pracma)
library(nleqslv)
DEBUG <- T
SAVE <- F

sysprop.sigma.fix <- function(tempC, system.properties, fitting.para) {
    # Sigma
    if (fitting.para$condensation) {
        lB <- 0.7E-9  # Bjerrum length at 298K
        system.properties$sigma[2] <- system.properties$size.ratio[2]*k.water.size / lB
    } 
    if (fitting.para$counterion.release) {
        # update Bjerrum length and the effective charge density of RNA
        lB <- ke^2 / (kEr*kkB*(tempC+273.15))
        system.properties$sigma[2] <- system.properties$size.ratio[2]*k.water.size / lB
    } 
    return(system.properties$sigma)
}

test.binodal.curve.3temps <- function(system.properties, fitting.para, ...) {
    
    system.properties$sigma <- sysprop.sigma.fix(tempC = 10, system.properties, fitting.para)
    
    p2 <- get.binodal.curves(10, system.properties$Chi, system.properties, fitting.para) %>% 
        mutate(conc.polymer = conc.p + conc.q)
    
    g4 <- ggplot(p2, aes(y = conc.salt, group = tempC)) +
        geom_point(aes(x = conc.p + conc.q, col = tempC), lwd = 2) +
        scale_colour_continuous(breaks = seq(20, 40, 10), guide = 'legend') +
        labs(x = 'Polymer [mol/L]',
             y = 'Salt [mol/L]',
             col = 'Temperature C' )
    g4 <- theme.title.text.1(g4)
    # g <- ggplot(p2) +
        # geom_line(aes(x = seqr(range(conc.polymer), 1e-4), y = spline(conc.polymer, conc.salt, xout = seqr(range(conc.polymer), 1e-4))$y, col = tempC))
    # print(g)
    print(g4)
}
test.binodal.curve.phi.3temps <- function(system.properties, fitting.para, 
                                          tempC = 10, ...) {
    
    system.properties$sigma <- sysprop.sigma.fix(tempC, system.properties, fitting.para)
    
    p2 <- get.binodal.curves(tempC, system.properties$Chi, system.properties, fitting.para)
    g4 <- ggplot(p2, aes(y = phi.salt, group = tempC)) +
        geom_point(aes(x = phi.polymer, col = tempC), lwd = 2) +
        scale_colour_continuous(breaks = seq(20, 40, 10), guide = 'legend') +
        labs(x = 'Polymer ',
             y = 'Salt ',
             col = 'Temperature C' )
    g4 <- theme.title.text.1(g4)
    print(g4)
}

DEBUG.kpq.fac              <- 1
fitting.para$binodal.guess <- 0.1
system.properties$Chi[1, 2] <- -0.
system.properties$Chi[2, 1] <- -0.


temp <- nleqslv(x = 273, function(x) get.alpha(x, k.water.size)-3.622 )$x
tempC <- temp - 273

g <- test.binodal.curve.phi.3temps(system.properties, fitting.para, tempC = tempC)
ggsave('test.binodal.curve.chin0.png', width = 5, height = 5)


if(SAVE) {
    ggsave('test.binodal.curve.temp.png', g2, width = 5, height = 5)
    ggsave('test.binodal.curve.sigma.png', g3, width = 5, height = 5)
    ggsave('test.binodal.curve.polymernum.png', g4, width = 5, height = 5)
} 