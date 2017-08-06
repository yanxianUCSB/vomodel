# Binodal curve at different temperature
rm(list = ls())
source('vomodel.R')
source('para.voorn.R')
library(yxplot)
library(ggplot2)
library(dplyr)
library(rootSolve)
library(pracma)
library(nleqslv)
DEBUG <- F
SAVE <- F

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


fitting.para$binodal.guess <- c(0.05, 0.00001)
p20 <- get.binodal.curve(75, system.properties$Chi, system.properties, fitting.para)




g2 <- ggplot(p20, aes(y = conc.salt)) +
  geom_line(aes(x = phi.polymer), lwd = 2) +
  labs(x = 'Polymer [%]',
       y = 'Salt [mol/L]',
       col = expression('Temperature'~degree~C))
g2 <- theme.title.text.1(g2)


print(g2)
if(SAVE) {
    ggsave('test.binodal.curve.voorn.png', g2, width = 6, height = 3)
} 