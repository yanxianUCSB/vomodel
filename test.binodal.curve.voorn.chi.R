# Binodal curve at different temperature
rm(list = ls())
source('vomodel.R')
source('para.voorn.chi.R')
library(yxplot)
library(ggplot2)
library(dplyr)
library(rootSolve)
library(pracma)
library(nleqslv)
DEBUG <- T
SAVE <- F


fitting.para$binodal.guess <- c(0.05, 0.03)
# fitting.para$sampling.start <- 0.0007875

ds <- NULL
# Chi1
for (chipw in c(0.01, -0.02)) {
system.properties$Chi <- matrix(c(.0,.0,0,0, chipw,
                                .0,.0,0,0,0,
                                0,0,0,0,0,
                                0,0,0,0,0,
                                chipw,0,0,0,0), 5,5)
p20 <- get.binodal.curve(40, system.properties$Chi, system.properties, fitting.para) %>% 
    mutate(Chi = chipw)
ds <- rbind(ds, p20)
}

# critical points
print(unique(ds$critic.salt))
# Plots
g2 <- ggplot(ds, aes(y = phi.salt*100)) +
    geom_point(aes(x = phi.polymer*100, col = as.factor(Chi)), lwd = 2) +
    # geom_line(aes(x = phi.polymer*100, group = pairing, col = as.factor(Chi))) +
    labs(x = 'Polymer [%]',
         y = 'Salt [%]',
         col = expression('Chi'))
g2 <- theme.title.text.1(g2)


print(g2)
if(SAVE) {
    ggsave('test.binodal.curve.voorn.png', g2, width = 6, height = 3)
} 