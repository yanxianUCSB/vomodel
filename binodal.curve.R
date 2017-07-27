# EXPLORE binodal curve
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

system.properties <- list(
  polymer.num = c(1000, 1000, 1, 1, 1),
  sigma = c(0.34, 0.34, 1, 1, 0),
  size.ratio = c(1, 1, 1, 1, 1),
  water.size = k.water.size
)
fitting.para <- list()
fitting.para$epsilon <- 1E-8
fitting.para$sampling.gap <- 1e-4
fitting.para$critical.point.guess <- c(phi.polymer=0.01, phi.salt=0.15)
fitting.para$binodal.guess <- c(1E-1, 1E-1)


# fitting.para$binodal.guess <- c(0.05, 0.05)
# p1 <- get.binodal.curve(20, 0, system.properties, fitting.para)
# fitting.para$binodal.guess <- c(0.1, 0.1)
# p3 <- get.binodal.curve(40, 0, system.properties, fitting.para)

p2 <- get.binodal.curves(seq(20, 40, 10), 0, system.properties, fitting.para)
p3 <- do.call(rbind, lapply(seq(0.34, 0.44, 0.05), function(sigma){
    system.properties$sigma <- c(sigma, sigma, 1, 1, 0)
    get.binodal.curve(30, 0, system.properties, fitting.para)
}))
p4 <- do.call(rbind, lapply(seq(500, 1000, 250), function(polymer.num){
    system.properties$polymer.num <- c(polymer.num, polymer.num, 1, 1, 0)
    get.binodal.curve(30, 0, system.properties, fitting.para)
}))


g2 <- ggplot(p2, aes(y = conc.salt, group = tempC)) +
  geom_line(aes(x = conc.p + conc.q, col = tempC), lwd = 2) +
    scale_colour_continuous(breaks = seq(20, 40, 10), guide = 'legend') +
  labs(x = 'Polymer [mol/L]',
       y = 'Salt [mol/L]',
       col = expression('Temperature'~degree~C))
g2 <- theme.title.text.1(g2)

g3 <- ggplot(p3, aes(y = conc.salt, group = sigma.p)) +
  geom_line(aes(x = conc.p + conc.q, col = sigma.p), lwd = 2) +
    scale_colour_continuous(breaks = seq(0.34, 0.44, 0.05), guide = 'legend') +
  labs(x = 'Polymer [mol/L]',
       y = 'Salt [mol/L]',
       col = expression(sigma) )
g3 <- theme.title.text.1(g3)

g4 <- ggplot(p4, aes(y = conc.salt, group = length.p)) +
  geom_line(aes(x = conc.p + conc.q, col = length.q), lwd = 2) +
    scale_colour_continuous(breaks = seq(500, 1000, 250), guide = 'legend') +
  labs(x = 'Polymer [mol/L]',
       y = 'Salt [mol/L]',
       col = 'Poly. Num.' )
g4 <- theme.title.text.1(g4)


if(SAVE) {
    ggsave('binodal.curve.temp.png', g2, width = 5, height = 5)
    ggsave('binodal.curve.sigma.png', g3, width = 5, height = 5)
    ggsave('binodal.curve.polymernum.png', g4, width = 5, height = 5)
} 