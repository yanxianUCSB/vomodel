rm(list = ls()[which(ls() != 'DEBUG.k.pd.sim')])
source('vomodel.R')
source('para.proteinRNA.R')
library(dplyr)
library(nleqslv)
library(ggplot2)
library(pracma)
library(mcGlobaloptim)

DEBUG <<- F
SAVE <<- T

k.conc.salt  <<- 0.030
k.conc.polymer  <<-  5E-6 * system.properties$MW[1] + 15E-3
k.sim.temp.range  <<- seq(4, 40, 1)
k.NULLPUNISH  <<- 1E10
# DEBUG.k.pd.sim <<- pd.sim

# sp <- system.properties
# k.phi.salt <<- k.conc.salt * 1000 * kNa * sp$polymer.num[3] * sp$size.ratio[3] * sp$water.size ^ 3 
# k.phi.polymer <<- 

get.phase.diagram.exp <- function(dataset.file = 'dataset.csv') {
    dataset <- read.csv(dataset.file) %>% 
        mutate(conc.polymer = protein * 1E-6 * system.properties$MW[1] + rna * 1E-3) %>% 
        mutate(conc.salt = nacl * 1E-3) %>% 
        mutate(tempC.cp = cloudpoint,
               tempC.on = onset) %>% 
        select(conc.polymer, conc.salt, tempC.cp, tempC.on)
    return(dataset)
        
}

chi.fn <- function(Chi.vec, ds.exp, system.properties, fitting.para) {
    Chi      <- matrix(rep(0, 25), 5, 5)
    Chi[1,2] <- Chi.vec[1]
    Chi[1,5] <- Chi.vec[2]
    Chi      <- Chi + t(Chi)
    system.properties$Chi <- Chi
    
    cat('current guess Chi\n')
    print(Chi)
    
    ds.sim      <- get.phase.diagram(system.properties, fitting.para)
    if (is.null(ds.sim)) return(k.NULLPUNISH)
    
    cat('current binodal curve \n')
    print(head(ds.sim[c('phi.polymer', 'phi.salt')]))
    
    ds.sim.conc <- get.phase.diagram.temp.conc(ds.sim, system.properties)
    ds.sim.nacl <- get.phase.diagram.temp.nacl(ds.sim, system.properties)
    predict.conc <- predict(lm(tempC ~ conc.polymer, ds.sim.conc), newdata = ds.exp)
    predict.nacl <- predict(lm(tempC ~ conc.salt, ds.sim.nacl), newdata = ds.exp)
    rmse.conc    <- rmserr(ds.exp$tempC.cp, predict.conc)$rmse
    rmse.conc[2] <- rmserr(ds.exp$tempC.on, predict.conc)$rmse
    rmse.nacl    <- rmserr(ds.exp$tempC.cp, predict.nacl)$rmse
    rmse.nacl[2] <- rmserr(ds.exp$tempC.on, predict.nacl)$rmse
    print(rmse.conc)
    print(rmse.nacl)
    # stop()
    # rmse.conc    <- rmserr(ds.exp$tempC.cp, spline(ds.sim.conc$conc.polymer, ds.sim.conc$tempC, xout = ds.exp$conc.polymer)$y )$rmse
    # rmse.conc[2] <- rmserr(ds.exp$tempC.on, spline(ds.sim.conc$conc.polymer, ds.sim.conc$tempC, xout = ds.exp$conc.polymer)$y )$rmse
    # rmse.nacl    <- rmserr(ds.exp$tempC.cp, spline(ds.sim.nacl$conc.salt,    ds.sim.nacl$tempC, xout = ds.exp$conc.salt   )$y )$rmse
    # rmse.nacl[2] <- rmserr(ds.exp$tempC.on, spline(ds.sim.nacl$conc.salt,    ds.sim.nacl$tempC, xout = ds.exp$conc.salt   )$y )$rmse
    return(sum(c(rmse.conc, rmse.nacl)))
}

ds.exp <- get.phase.diagram.exp(dataset.file=choose.files(caption='Select Exp Dataset', multi = F))

out <- multiStartoptim(c(0.1, -0.03), chi.fn, ds.exp=ds.exp, system.properties=system.properties, fitting.para=fitting.para)

if (SAVE) {
    cat('saving data >>>\n')
    saveRDS(out, 'out.fit.dat', ascii = T)
    cat('saved!')
}
