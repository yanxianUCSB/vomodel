rm(list = ls())
source('vomodel.R')
source('para.proteinRNA.R')

DEBUG <<- F
SAVE <<- F

k.conc.salt  <<- 0.030
k.conc.polymer  <<-  5E-6 * system.properties$MW[1] + 15E-3
# sp <- system.properties
# k.phi.salt <<- k.conc.salt * 1000 * kNa * sp$polymer.num[3] * sp$size.ratio[3] * sp$water.size ^ 3 
# k.phi.polymer <<- 

get.phase.diagram <- function(system.properties, fitting.para) {
    sampling <- list(
        tempC = seq(4, 60, 1)
    )
    p <- lapply(sampling$tempC, function(tempC) {
        if (fitting.para$condensation == T) {
            lB <- 0.7E-9  # Bjerrum length at 298K
            system.properties$sigma[2] <- system.properties$size.ratio[2]*k.water.size / lB
        }
        if (fitting.para$counterion.release == T) {
            # update Bjerrum length and the effective charge density of RNA
            lB <- ke^2 / (kEr*kkB*(tempC+273.15))
            system.properties$sigma[2] <- system.properties$size.ratio[2]*k.water.size / lB
        }
        # update critical.point.guess
        fitting.para$critical.point.guess <- as.numeric(fitting.para$c.point.temp.fun(tempC + 273))
        cat(paste0('Temp [C]: ', tempC, '\n'))
        cat(paste0('Bjerrum length [m]: ', lB, '\n'))
        cat(paste0('Sigma RNA: ', system.properties$sigma[2], '\n'))
        cat('fitting >>>\n')
        out <- get.binodal.curve(tempC, Chi = system.properties$Chi, system.properties, fitting.para, unit = 'mol')
        cat('succeeded!\n')
        return(out)
    })
    
    out <- do.call(rbind, p) %>% 
        filter(phase.separated)
    return(out) 
}
get.phase.diagram.temp.conc <- function(phase.diagram.ds, system.properties) {
    precision <- 1e-5
    
    ds <- phase.diagram.ds
    ds <- 
        as.data.frame(ds) %>% 
        mutate(conc.polymer = conc.p * system.properties$MW[1] + conc.q * system.properties$MW[2])
    
    ds2 <- lapply(unique(ds$tempC), function(tempc) {
        ds.t1 <- ds %>% filter(tempC == tempc, phase == 'dilute')
        ds.t2 <- ds %>% filter(tempC == tempc, phase == 'dense')
        ds3 <- rbind(
            data.frame(conc.salt = 0.030) %>%
                mutate(
                    conc.polymer = spline(ds.t1$conc.salt, ds.t1$conc.polymer, xout = conc.salt)$y,
                    phase = 'dilute'
                ),
            data.frame(conc.salt = 0.030) %>%
                mutate(
                    conc.polymer = spline(ds.t2$conc.salt, ds.t2$conc.polymer, xout = conc.salt)$y,
                    phase = 'dense'
                )
        ) %>%
            mutate(
                conc.salt = 0.030,
                tempC = tempc
            ) 
    })
    
    ds3 <- do.call(rbind, ds2) 
    
    ds4 <- ds3 %>% filter(conc.polymer < 20, conc.polymer > 0) %>% select(conc.polymer, tempC)
    return(ds4)
    
}
get.phase.diagram.temp.nacl <- function(phase.diagram.ds, system.properties, nacl.range = NULL) {
    precision <- 1e-5
    
    ds <- phase.diagram.ds
    ds <- 
        as.data.frame(ds) %>% 
        mutate(conc.polymer = conc.p * system.properties$MW[1] + conc.q * system.properties$MW[2])
    
    ds2 <- lapply(unique(ds$tempC), function(tempc) {
        ds.t1 <- ds %>% filter(tempC == tempc, phase == 'dilute') 
        
        if (min(ds.t1$conc.salt) > k.conc.salt) return()
        
        # if (DEBUG){
        #     plot(ds.t1$conc.polymer, ds.t1$conc.salt)
        #     readline('>>>')
        # }
        ds3 <- rbind(
            data.frame(conc.polymer = k.conc.polymer) %>%
                mutate(
                    conc.salt = spline(ds.t1$conc.polymer, ds.t1$conc.salt, xout = conc.polymer)$y,
                    phase = 'dilute'
                )
        ) %>%
            mutate(
                conc.polymer = k.conc.polymer,
                tempC = tempc
            ) 
    })
    
    ds3 <- do.call(rbind, ds2) 
    
    ds4 <- ds3 %>% filter(conc.polymer < 20, conc.polymer > 0) %>% select(conc.salt, tempC)
    return(ds4)
    
}

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

get.base.chi.fun <- function(Chi.vec, pd.exp, system.properties, fitting.para) {
    # Chi.vec = c(chi.pp, chi.pq, chi.pw, chi.qw)
    chi.pp <- Chi.vec[1]
    chi.pq <- Chi.vec[2]
    chi.pw <- Chi.vec[3]
    chi.qw <- Chi.vec[4]
    system.properties$Chi <- matrix(c(
        chi.pp,    chi.pq,    0, 0, chi.pw,
        chi.pq,    0,         0, 0, chi.qw,
        0,         0,         0, 0,      0,
        0,         0,         0, 0,      0,
        chi.pw,    chi.qw,    0, 0,      0
    ))
    pd.sim <- get.phase.diagram(system.properties, fitting.para)
    pd.sim.con <- get.phase.diagram.temp.conc(pd.sim, system.properties)
    pd.sim.nacl <- get.phase.diagram.temp.nacl(pd.sim, system.properties)
    rmse.conc <- rmserr(pd.exp$tempC.cp, spline(pd.sim.con$conc.polymer, pd.sim.con$tempC, xout = pd.exp$conc.polymer)$y)$rmse
    rmse.conc.2 <- rmserr(pd.exp$tempC.on, spline(pd.sim.con$conc.polymer, pd.sim.con$tempC, xout = pd.exp$conc.polymer)$y)$rmse
    rmse.nacl <- rmserr(pd.exp$tempC.cp, spline(pd.sim.con$conc.salt, pd.sim.con$tempC, xout = pd.exp$conc.salt)$y)$rmse
    rmse.nacl.2 <- rmserr(pd.exp$tempC.on, spline(pd.sim.con$conc.salt, pd.sim.con$tempC, xout = pd.exp$conc.salt)$y)$rmse
    
    return(c(rmse.conc, rmse.conc.2, rmse.nacl, rmse.nacl.2))
}




system.properties$Chi = matrix(c(
    0.0, -0.08, 0, 0, 0.01,
    -0.08, 0.0, 0, 0, 0.001,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0.01, 0.001, 0, 0, 0
), 5, 5)

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