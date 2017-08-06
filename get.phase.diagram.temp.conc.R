rm(list = ls())
source('vomodel.R')
source('para.proteinRNA.R')

DEBUG <<- T
SAVE <<- F

get.phase.diagram <- function(system.properties, fitting.para) {
    sampling <- list(
        tempC = seq(22, 60, 1)
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
    return(do.call(rbind, p))
}
get.phase.diagram.temp.conc <- function(phase.diagram.ds, system.properties) {
    precision <- 1e-5
    
    ds <- phase.diagram.ds
    ds <- 
        as.data.frame(ds) %>% mutate(conc.polymer = conc.p * system.properties$MW[1] + conc.q * system.properties$MW[2])
    
    ds2 <- lapply(unique(ds$tempC), function(tempc) {
        ds.t1 <- ds %>% filter(tempC == tempc, phase == 'dilute')
        ds.t2 <- ds %>% filter(tempC == tempc, phase == 'dense')
        ds3 <- rbind(
            data.frame(phi.salt = seq(1e-8, ds.t1$critic.salt[1], precision)) %>%
                mutate(
                    phi.polymer = spline(ds.t1$phi.salt, ds.t1$phi.polymer, xout = phi.salt)$y,
                    phase = 'dilute'
                ),
            data.frame(phi.salt = seq(ds.t2$critic.salt[1], 1e-8, - precision)) %>%
                mutate(
                    phi.polymer = spline(ds.t2$phi.salt, ds.t2$phi.polymer, xout = phi.salt)$y,
                    phase = 'dense'
                )
        ) %>%
            mutate(
                conc.polymer = phi.polymer * ds.t2$conc.polymer[1] / ds.t2$phi.polymer[1],
                conc.salt = phi.salt * ds.t2$conc.salt[1] / ds.t2$phi.salt[1],
                tempC = tempc
            ) %>%
            filter(abs(phi.polymer) < 1)
    })
    
    ds3 <- do.call(rbind, ds2) 
    return(ds3)
    
}

test <- function(s, f){
    ds <- get.phase.diagram(s, f)
    print(head(ds))
    return(ds)
}

ds <- test(system.properties, fitting.para)

ds2 <- get.phase.diagram.temp.conc(ds, system.properties)
print(ds2)
    g <- ggplot(ds2 %>% filter(conc.salt > 0.029, conc.salt < 0.031, conc.polymer < 20, conc.polymer > 0), aes(x = conc.polymer, y = tempC)) + geom_point()
    print(g)
    