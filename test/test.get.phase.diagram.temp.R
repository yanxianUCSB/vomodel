rm(list = ls())
source('vomodel.R')
source('para.proteinRNA.R')
library(dplyr)
library(nleqslv)
library(ggplot2)
library(yxplot)

DEBUG <<- T
SAVE <<- T

k.conc.salt       <<- 0.030
k.conc.polymer    <<-  5E-6 * system.properties$MW[1] + 15E-3
k.sim.temp.range  <<- seq(4, 16, 2)

# Water size
watersize <<- 3.1E-10
# k.amino.acid.length       <<- watersize * 1.0
# k.dna.contour.unit.length <<- watersize * 1.0
# k.na.size                 <<- watersize * 1.1
# k.cl.size                 <<- watersize * 1.1
# k.water.size              <<- watersize
system.properties$size.ratio <-  c(k.amino.acid.length, k.dna.contour.unit.length, k.na.size, k.cl.size, k.water.size) / k.water.size


# Chi
Chi <-  matrix(rep(0, 25), 5, 5)
Chi[1,2] <- 0.0875303
Chi[1,5] <- 0.070017
Chi[2,1] <- Chi[1,2]
Chi[5,1] <- Chi[1,5]
system.properties$Chi <- Chi



# # #

fitting.para$binodal.guess <- 0.015

phase.diagram.exp <- get.phase.diagram.exp(dataset.file = 'C:/Users/Yetsun/Box/anywhere/dataset.csv')

ds <- get.phase.diagram(system.properties, fitting.para)
saveRDS(ds, 'out.test.ds.data')
ds <- readRDS('out.test.ds.data')

ds2 <- get.phase.diagram.temp.conc(ds, system.properties, k.conc.salt = k.conc.salt)
ds3 <- get.phase.diagram.temp.nacl(ds, system.properties, k.conc.polymer = k.conc.polymer)

assertthat::assert_that(!is.null(ds))

g <- ggplot(ds, aes(x = conc.polymer, y = conc.salt, group = tempC)) +
    geom_point(aes(col = tempC)) +
    geom_line(aes(x = conc.polymer, y = conc.salt, col = tempC))

g2 <- ggplot(ds2 %>% filter(conc.polymer < 20), aes(x = conc.polymer, y = tempC)) + 
    geom_line(aes(col = 'sim'), lwd = 2) +
    geom_point(data = phase.diagram.exp %>% filter(conc.salt == 0.03), aes(x = conc.polymer, y = tempC.cp, col = 'exp.cp')) +
    geom_point(data = phase.diagram.exp %>% filter(conc.salt == 0.03), aes(x = conc.polymer, y = tempC.on, col = 'exp.on')) +
    labs(x = 'Conc.polymer[mg/mL]', 
         y = 'Temperature [C]', 
         col = '')

g3 <- ggplot(ds3, aes(x = conc.salt, y = tempC)) + 
    geom_line(aes(col = 'sim'), lwd = 2) +
    geom_point(data = phase.diagram.exp %>% filter(abs(conc.polymer - 0.125) < 1e-3), aes(x = conc.salt, y = tempC.cp, col = 'exp.cp')) +
    geom_point(data = phase.diagram.exp %>% filter(abs(conc.polymer - 0.125) < 1e-3), aes(x = conc.salt, y = tempC.on, col = 'exp.on')) +
    labs(x = 'NaCl [M]', 
         y = 'Temperature [C]', 
         col = '')

print(ds2)
print(g)
stop()
# readline('>>> ')
print(g2)
readline('>>> ')
print(g3)
readline('>>> ')

if(SAVE) {
    g2 <- theme.background.1(g2)
    g2 <- theme.title.text.2(g2)
    ggsave( 'get.phase.diagram.temp.conc.png', g2, width = 5, height = 5)
    g3 <- theme.background.1(g3)
    g3 <- theme.title.text.2(g3)
    ggsave( 'get.phase.diagram.temp.nacl.png', g3, width = 5, height = 5)
}

