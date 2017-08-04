# Binodal curve at different temperature
rm(list = ls())
source('vomodel.R')
source('proteinRNA.para.R')
library(yxplot)
library(ggplot2)
library(dplyr)
library(rootSolve)
library(pracma)
library(nleqslv)
DEBUG <- F
SAVE <- T

# system.properties <- list(
#   polymer.num = c(1000, 1000, 1, 1, 1),
#   sigma = c(0.34, 0.34, 1, 1, 0),
#   size.ratio = c(1, 1, 1, 1, 1),
#   water.size = k.water.size,
#   molar.ratio = rep(1, 5)
# )
# fitting.para <- list()
# fitting.para$epsilon <- 1E-8
# fitting.para$sampling.gap <- 1e-4
# fitting.para$critical.point.guess <- c(phi.polymer=0.01, phi.salt=0.15)
# fitting.para$binodal.guess <- c(1E-1, 1E-1)


# fitting.para$binodal.guess <- c(0.05, 0.05)
ds <- do.call(rbind, lapply(c(20, 30, 40), function(tempC){
    lB <- ke^2 / (kEr*kkB*(tempC+273.15))
    system.properties$sigma[2] <- system.properties$size.ratio[2]*k.water.size / lB
    # update critical.point.guess
    fitting.para$critical.point.guess <- as.numeric(fitting.para$c.point.temp.fun(tempC + 273))
    return(get.binodal.curve(tempC, system.properties$Chi, system.properties, fitting.para))
}))
# p1 <- get.binodal.curve(20, system.properties$Chi, system.properties, fitting.para)
# p2 <- get.binodal.curve(30, system.properties$Chi, system.properties, fitting.para)
# p3 <- get.binodal.curve(40, system.properties$Chi, system.properties, fitting.para)
# ds <- rbind(p1, p2, p3)

g2 <- ggplot(ds, aes(y = conc.salt, group = tempC)) +
  geom_line(aes(x = conc.p + conc.q, col = tempC), lwd = 2) +
    scale_colour_continuous(breaks = seq(20, 40, 10), guide = 'legend') +
  labs(x = 'Polymer [mol/L]',
       y = 'Salt [mol/L]',
       col = expression('Temperature'~degree~C))
g2 <- theme.title.text.1(g2)

print(g2)

if(SAVE) {
    ggsave('test.binodal.curve.proteinRNA.temp.png', g2, width = 5, height = 5)
} 