# compare with published data binodal curve
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
  sigma = c(0.44, 0.44, 1, 1, 0),
  size.ratio = c(1, 1, 1, 1, 1),
  water.size = k.water.size,
  molar.ratio = rep(1, 5)
)
fitting.para <- list()
fitting.para$epsilon <- 1E-8
fitting.para$sampling.gap <- 1e-4
fitting.para$critical.point.guess <- c(phi.polymer=0.01, phi.salt=0.15)
fitting.para$binodal.guess <- c(1E-1, 1E-1)

temp <- ke ^ 2 / ( (3.655 /  (2 / 3 * sqrt(pi) ) ) ^ (2/3)  * k.water.size)  / (kEr * kkB)


p2 <- get.binodal.curve(temp - 273.15, 0, system.properties, fitting.para) %>%
    mutate(ID = 'Yanxian') %>%
    select(phi.polymer, phi.salt, ID)


# CHIMAD standard
# CHIMAD of sigma 0.44 alpha 3.655
chimad <- read.delim('ref.chimad.txt') %>%
    rename(phi.polymer = VolumeFraction, phi.salt = Binodal) %>%
    mutate(ID = 'CHIMAD') %>%
    select(phi.polymer, phi.salt, ID)

# plot and save 
ds <- rbind(p2, chimad) %>% 
    ungroup()
head(ds)

g <- ggplot(ds, aes(y = phi.salt, group = ID)) +
  geom_line(aes(x = phi.polymer, col = ID), lwd = 2, lty = 2, alpha = 0.8) +
  labs(x = 'Polymer phi',
       y = 'Salt psi',
       col = paste0('RMSE = ', format(rmserr(x = spline(ds$phi.polymer[which(ds$ID == 'Yanxian')], 
                                      ds$phi.salt[which(ds$ID == 'Yanxian')], 
                                      xout = ds$phi.polymer[which(ds$ID == 'CHIMAD')])$y,
                           y = ds$phi.salt[which(ds$ID == 'CHIMAD')])$rmse, digits = 3)))
g <- theme.title.text.1(g)
print(g)
if(SAVE) ggsave('test.binodal.curve.verify.png', width = 5, height = 5)
