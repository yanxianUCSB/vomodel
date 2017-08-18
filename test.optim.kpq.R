# optimize DEBUG
rm(list = ls())
source('vomodel.R')
source('para.proteinRNA.R')
library(yxplot)
library(yxhelper)
library(ggplot2)
library(dplyr)
library(rootSolve)
library(pracma)
library(nleqslv)
DEBUG <- T
SAVE <- F
# DEBUG.good.size.p <- 1
DEBUG.good.kpqfac <- 15
DEBUG.good.guess <- 0.008


test.binodal.curve <- function(system.properties, fitting.para, DEBUG.guess, tempC = 10, ...) {
    
    
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
    
    # Chi
    Chi <-  matrix(rep(0, 25), 5, 5)
    Chi[1,2] <- 0.
    Chi[1,5] <- -0.
    Chi[2,1] <- Chi[1,2]
    Chi[5,1] <- Chi[1,5]
    system.properties$Chi <- Chi
    fitting.para$binodal.guess <- DEBUG.guess
    
    # size
    size <- system.properties$size.ratio
    # system.properties$size.ratio[1] <- 2
    # system.properties$size.ratio[2] <- 1
    # system.properties$size.ratio[3] <- 1
    # system.properties$size.ratio[4] <- 1
    
    out <- get.binodal.curves(tempC, system.properties$Chi, system.properties, fitting.para)
    
}


# Find guess
for (k in 1:1000) {
    
    
    DEBUG.guess.seq <- seq(1e-3, 1e-2, 1e-4)
    y <- sapply(DEBUG.guess.seq, function(i) {
        DEBUG.kpq.fac <<- DEBUG.good.kpqfac
        # system.properties$size.ratio[1] <- DEBUG.good.size.p
        # system.properties$size.ratio[2] <- DEBUG.good.size.p
        # system.properties$size.ratio[3] <- DEBUG.good.size.p
        # system.properties$size.ratio[4] <- DEBUG.good.size.p
        out <- test.binodal.curve(system.properties, fitting.para, DEBUG.guess = i)
        if (is.null(out)) return(NA)
        else return(i)
    })
    # print(DEBUG.guess.seq)
    print(y)
    DEBUG.good.guess <- min(y[!is.na(y)])
    
    
    # Find kpq.fac
    
    DEBUG.guess.seq <- seq(1, 20, 1e-1)
    y <- sapply(DEBUG.guess.seq, function(i) {
        DEBUG.kpq.fac <<- i
        # system.properties$size.ratio[1] <- DEBUG.good.size.p
        # system.properties$size.ratio[2] <- DEBUG.good.size.p
        # system.properties$size.ratio[3] <- DEBUG.good.size.p
        # system.properties$size.ratio[4] <- DEBUG.good.size.p
        out <- test.binodal.curve(system.properties, fitting.para, DEBUG.guess = DEBUG.good.guess)
        if (is.null(out)) return(NA)
        else return(i)
    })
    print(y)
    DEBUG.good.kpqfac <- min(y[!is.na(y)])
    
    # Find size
    
    # DEBUG.guess.seq2 <- seq(1, 2, 0.1)
    # y2 <- sapply(DEBUG.guess.seq2, function(j){
    #     DEBUG.kpq.fac <<- DEBUG.good.kpqfac
    #     system.properties$size.ratio[1] <- j
    #     system.properties$size.ratio[2] <- j
    #     system.properties$size.ratio[3] <- j
    #     system.properties$size.ratio[4] <- j
    #     out <- test.binodal.curve(system.properties, fitting.para, DEBUG.guess = DEBUG.good.guess)
    #     if (is.null(out)) return(NA)
    #     else return(j)
    # })
    # print(y2)
    # DEBUG.good.size.p <- max(y2[!is.na(y2)])
    
}

DEBUG.kpq.fac <<- DEBUG.good.kpqfac

out <- test.binodal.curve(system.properties, fitting.para, DEBUG.guess = DEBUG.good.guess)

plot(out$conc.p + out$conc.q, out$conc.salt)

plot(out$phi.polymer, out$phi.salt)

# fitting.para$binodal.guess <- c(0.05, 0.05)
# p1 <- get.binodal.curve(20, system.properties$Chi, system.properties, fitting.para)
# fitting.para$binodal.guess <- c(0.1, 0.1)
# p3 <- get.binodal.curve(40, system.properties$Chi, system.properties, fitting.para)
# 
# p3 <- do.call(rbind, lapply(seq(0.34, 0.44, 0.05), function(sigma){
#     system.properties$sigma <- c(sigma, sigma, 1, 1, 0)
#     get.binodal.curve(30, system.properties$Chi, system.properties, fitting.para)
# }))
# p4 <- do.call(rbind, lapply(seq(500, 1000, 250), function(polymer.num){
#     system.properties$polymer.num <- c(polymer.num, polymer.num, 1, 1, 0)
#     get.binodal.curve(30, system.properties$Chi, system.properties, fitting.para)
# }))


# g2 <- ggplot(p2, aes(y = conc.salt, group = tempC)) +
#   geom_line(aes(x = conc.p + conc.q, col = tempC), lwd = 2) +
#     scale_colour_continuous(breaks = seq(20, 40, 10), guide = 'legend') +
#   labs(x = 'Polymer [mol/L]',
#        y = 'Salt [mol/L]',
#        col = expression('Temperature'~degree~C))
# g2 <- theme.title.text.1(g2)
# 
# g3 <- ggplot(p3, aes(y = conc.salt, group = sigma.p)) +
#   geom_line(aes(x = conc.p + conc.q, col = sigma.p), lwd = 2) +
#     scale_colour_continuous(breaks = seq(0.34, 0.44, 0.05), guide = 'legend') +
#   labs(x = 'Polymer [mol/L]',
#        y = 'Salt [mol/L]',
#        col = expression(sigma) )
# g3 <- theme.title.text.1(g3)
