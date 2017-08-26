# manning's condensation
rm(list = ls())
source('vomodel.R')
source('para.proteinRNA.R')
library(yxplot)
library(ggplot2)
library(dplyr)
library(rootSolve)
library(pracma)
library(nleqslv)
DEBUG <- T
SAVE <- T

test.mannings <- function(system.properties, fitting.para) {
    # system.properties <- list(
    #     polymer.num   = c(207, 900e3 / 306.2,   1, 1, 1),
    #     sigma         = c(11 / 207, 0.5,   1, 1, 0),
    #     size.ratio    = c(1,    1,   1, 1, 1),
    #     molar.ratio   = c(900e3 / 306.2,  11, 0.5, 0.5, 0),
    #     # the polycation:polyanion and cation:anion molar ratio
    #     water.size    = k.water.size
    # )
    # fitting.para <- list()
    # fitting.para$epsilon <- 1E-8
    # fitting.para$sampling.gap <- 1e-4
    # fitting.para$critical.point.guess <-
    #     c(phi.polymer = 0.005, phi.salt = 0.005)
    # fitting.para$c.point.temp.fun <-
    #     c.point.temp.fun(c.point.temp(system.properties, fitting.para))
    # fitting.para$binodal.guess <-
    #     c(0.1,  0.005)  # phi.polymer.2, phi.salt = 0.9 * critical salt
    
    temp <- 298
    bjerrum.length <- ke ^ 2 / (kEr * kkB * temp)
    sigma.dna <- k.dna.contour.unit.length / bjerrum.length
    
    # mannings
    system.properties$sigma <- c(11 / 207, sigma.dna,   1, 1, 0)
    ds.mannings <-
        get.binodal.curve(
            20,
            Chi = system.properties$Chi,
            sysprop = system.properties,
            fitting.para = fitting.para,
            unit = 'mol'
        ) %>% mutate(mannings = T)
    
    # not mannings
    system.properties$sigma <- c(11 / 207, 1,   1, 1, 0)
    fitting.para$critical.point.guess <-
        c(phi.polymer = 0.005, phi.salt = 0.2)
    ds.flory <-
        get.binodal.curve(
            20,
            Chi = system.properties$Chi,
            sysprop = system.properties,
            fitting.para = fitting.para,
            unit = 'mol'
        ) %>% mutate(mannings = F)
    
    return(rbind(ds.mannings, ds.flory))
}

ds <- test.mannings(system.properties, fitting.para)

g <- ggplot(ds, aes(x = conc.p, y = conc.salt, group = mannings)) +
    geom_line(aes(col = mannings), lwd = 2) +
    labs(x = 'Protein Conc [mol/L]',
         y = 'Salt Conc [mol/L]',
         col = 'Mannings') +
    scale_y_log10()
g <- theme.background.1(g)
g <- theme.title.text.1(g)
print(g)

if (SAVE)
    ggsave('test.mannings.png', width = 5, height = 5)