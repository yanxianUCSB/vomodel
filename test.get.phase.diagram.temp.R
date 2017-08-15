rm(list = ls())
source('vomodel.R')
source('para.proteinRNA.R')
library(dplyr)
library(nleqslv)
library(ggplot2)

DEBUG <<- T
SAVE <<- F

k.conc.salt  <<- 0.030
k.conc.polymer  <<-  5E-6 * system.properties$MW[1] + 15E-3
k.sim.temp.range  <<- seq(4, 40, 1)

get.phase.diagram.exp <- function(dataset.file = 'dataset.csv') {
    dataset.file <- '~/Box/anywhere/dataset.csv'
    dataset <- read.csv(dataset.file) %>% 
        mutate(conc.polymer = protein * 1E-6 * system.properties$MW[1] + rna * 1E-3) %>% 
        mutate(conc.salt = nacl * 1E-3) %>% 
        mutate(tempC.cp = cloudpoint,
               tempC.on = onset) %>% 
        select(conc.polymer, conc.salt, tempC.cp, tempC.on)
    return(dataset)
    
}

Chi <-  matrix(rep(0, 25), 5, 5)
Chi[1,2] <- 0
Chi[1,5] <- 0
Chi[2,1] <- Chi[1,2]
Chi[5,1] <- Chi[1,5]

system.properties$Chi <- Chi

phase.diagram.exp <- get.phase.diagram.exp()

ds <- get.phase.diagram(system.properties, fitting.para)

ds2 <- get.phase.diagram.temp.conc(ds, system.properties)

ds3 <- get.phase.diagram.temp.nacl(ds, system.properties)


g <- ggplot(ds %>% filter(conc.salt < 0.1), aes(x = conc.p, y = conc.salt, group = tempC)) +
    geom_point(aes(col = tempC))

g2 <- ggplot(ds2, aes(x = conc.polymer, y = tempC)) + geom_point(aes(col = 'sim')) +
    geom_point(data = phase.diagram.exp %>% filter(conc.salt == 0.03), aes(x = conc.polymer, y = tempC.cp, col = 'exp.cp')) +
    geom_point(data = phase.diagram.exp %>% filter(conc.salt == 0.03), aes(x = conc.polymer, y = tempC.on, col = 'exp.on')) +
    labs(x = 'Conc.polymer[mg/mL]', 
         y = 'Temperature [C]', 
         col = '')

g3 <- ggplot(ds3, aes(x = conc.salt, y = tempC)) + geom_point(aes(col = 'sim')) +
    geom_point(data = phase.diagram.exp %>% filter(abs(conc.polymer - 0.125) < 1e-3), aes(x = conc.salt, y = tempC.cp, col = 'exp.cp')) +
    geom_point(data = phase.diagram.exp %>% filter(abs(conc.polymer - 0.125) < 1e-3), aes(x = conc.salt, y = tempC.on, col = 'exp.on')) +
    labs(x = 'NaCl [M]', 
         y = 'Temperature [C]', 
         col = '')

print(ds2)
print(g)
print(g2)
print(g3)

if(SAVE) {
    g2 <- theme.background.1(g2)
    g2 <- theme.title.text.2(g2)
    ggsave( 'get.phase.diagram.temp.conc.png', g2, width = 5, height = 5)
    g3 <- theme.background.1(g3)
    g3 <- theme.title.text.2(g3)
    ggsave( 'get.phase.diagram.temp.nacl.png', g3, width = 5, height = 5)
}