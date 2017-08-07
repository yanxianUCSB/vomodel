rm(list = ls())
source('vomodel.R')
source('para.proteinRNA.R')

DEBUG <<- T
SAVE <<- F

# sp <- system.properties
# k.conc.salt  <<- 0.030
# k.phi.salt <<- k.conc.salt * 1000 * kNa * sp$polymer.num[3] * sp$size.ratio[3] * sp$water.size ^ 3 
# k.conc.polymer  <<-  5E-6 * system.properties$MW[1] + 15E-3
# k.phi.polymer <<- 
get.binodal.curve_.fun <- function(x, phi.polymer, system.properties, fitting.para) {
    # x[1]: phi.salt
    # x[2]: temp
    phi.salt <- x[1]
    temp <- x[2]
    
    alpha <- get.alpha(temp, size = system.properties$water.size
    
    
    phi.polymer.1 <- x[1]
    phi.salt <- x[2]
    phi.polymer <- c(phi.polymer.1, phi.polymer.2)
    return(
        c(
            gibbs.d(phi.polymer[1], phi.salt, ...) - gibbs.d(phi.polymer[2], phi.salt, ...),
            (phi.polymer[2] - phi.polymer[1]) * gibbs.d(phi.polymer[1], phi.salt, ...) - 
                (gibbs(phi.polymer[2], phi.salt, ...) - gibbs(phi.polymer[1], phi.salt, ...))
        )
    )
}

get.binodal.curve__ <- function( sysprop = NULL, fitting.para = NULL, ...) {
    #' generate binodal curve
    if (is.null(sysprop)) arg <- list(...)
    else arg <- sysprop
    
    c.point <- critical.point_(fitting.para$critical.point.guess, ...)
    assertthat::assert_that(c.point$phi.polymer > 0)
    # if (DEBUG) print(c('critical point', c.point))
    
    # search binodal point return c(phi.polymer, phi.salt)
    # phi.polymer.2.seq <- c(seq( arg$sampling.gap, arg$sampling.gap*1e3, arg$sampling.gap),
    #                        seq( arg$sampling.gap * 1e3, c.point$phi.polymer, arg$sampling.gap * 1e4))
    n <- floor(0.5 * (1 + sqrt(1 + 8 * c.point$phi.polymer / fitting.para$sampling.gap)))
    phi.polymer.2.seq <-  fitting.para$sampling.gap * seq(1, n, 1) * seq(2, n + 1, 1) * 0.5
    phi.polymer.2.seq <-  phi.polymer.2.seq[which(phi.polymer.2.seq >= fitting.para$sampling.start)]
    
    # Binodal Guess
    binodal.guess <- fitting.para$binodal.guess
    binodal.guess[2] <- c.point$phi.salt * 0.9
    
    output <- list()
    
    for (phi.polymer.2 in phi.polymer.2.seq) {
        
        test <- binodal.curve.fun_(x = binodal.guess, phi.polymer.2 = phi.polymer.2, ...)
        if (anyNA(test) ||
            is.nan(test[1]) ||
            is.nan(test[2])) {
            next
        }
        roots <- nleqslv (
            binodal.guess,
            binodal.curve.fun_,
            binodal.curve.jacobian_,
            phi.polymer.2 = phi.polymer.2, 
            ...
        )
        if (
            roots$x[1] > 0 &&
            roots$x[1] > phi.polymer.2 &&
            roots$x[2] > 0 &&
            roots$x[2] < c.point$phi.salt &&
            max(abs(roots$fvec)) < 1 &&
            roots$termcd == 2) {
            
            binodal.guess <- roots$x
            
            output[[length(output) + 1]] <-
                c(
                    phi.polymer.1 = roots$x[1],
                    phi.polymer.2 = phi.polymer.2,
                    phi.salt = roots$x[2],
                    f1 = roots$fvec[1],
                    f2 = roots$fval[2],
                    pair = phi.polymer.2
                )
        } else {
            # if(DEBUG) print(binodal.guess)
            # if(DEBUG) print(c.point$phi.polymer)
            # if(DEBUG) print(roots)
            binodal.guess <- binodal.guess * 1.01
            next
        }
        # if(DEBUG) print(roots$x)
        # if(DEBUG) print(roots$termcd)
        # if(DEBUG) print(binodal.guess)
    }
    
    if (length(output) == 0) {
        return()
    }
    assertthat::assert_that(length(output) > 0, msg = 'binodal.curve_ failed')
    
    p2 <- as.data.frame.matrix(do.call(rbind, output)) %>%
        filter(phi.polymer.1 > max(phi.polymer.2))  # requiring the dense phase should be larger than dilute phase
    
    assertthat::assert_that(nrow(p2) > 0)
    
    ds <- data.frame(
        phi.polymer = c(p2$phi.polymer.2, rev(p2$phi.polymer.1)),
        phi.salt = c(p2$phi.salt, rev(p2$phi.salt)),
        phase = factor(c(rep('dilute', nrow(p2)), rep('dense', nrow(p2))),
                       levels = c('dilute', 'dense')),
        critic.polymer = c.point$phi.polymer,
        critic.salt = c.point$phi.salt,
        pairing = c(p2$phi.polymer.2, rev(p2$phi.polymer.2)) * (p2$phi.salt[1] + 1),
        phase.separated = ifelse(nrow(p2) < 3, F, T)  # if the roots only has few rows then we say there is no binodal curve
    )
    return(ds)
}

get.binodal.curve_ <- function(tempC,
         Chi = 0,
         sysprop,
         fitting.para,
         unit = 'mol') {
    #' get binodal curve by Voorn Overbeek Model 3D component system
    #' temp in degree C
    #' Chi = 0 or 5 * 5 matrix
    temp <- tempC + 273.15
    
    kpq <- pkpq(sysprop = sysprop)
    
    ds <- binodal.curve_(
        sysprop = sysprop,
        fitting.para = fitting.para,
        temp = temp,
        Chi = Chi,
        alpha = get.alpha(temp, sysprop$water.size),
        sigma = sysprop$sigma,
        polymer.num = sysprop$polymer.num,
        size.ratio = sysprop$size.ratio,
        molar.ratio = sysprop$molar.ratio
    ) 
    if (is.null(ds)) return()
    ds <- ds %>%
        mutate(
            temp = temp,
            kpq = kpq,
            sigma.p = sysprop$sigma[1],
            sigma.q = sysprop$sigma[2],
            length.p = sysprop$polymer.num[1],
            length.q = sysprop$polymer.num[2]
        )
    
    if (unit == 'mol') {
        ds <- ds %>%
            rowwise() %>%
            mutate(
                tempC = tempC,
                conc.salt = phi.salt / (
                    kNa * sysprop$polymer.num[3] *
                        sysprop$size.ratio[3] *
                        sysprop$water.size ^ 3
                ) / 1000,
                # mol/L
                conc.p = phi.polymer / (1 + kpq) / (
                    kNa * sysprop$polymer.num[1] *
                        sysprop$size.ratio[1] *
                        sysprop$water.size ^ 3
                ) / 1000,
                conc.q = phi.polymer * kpq / (1 + kpq) / (
                    kNa * sysprop$polymer.num[2] *
                        sysprop$size.ratio[2] *
                        sysprop$water.size ^ 3
                ) / 1000
            ) %>%
            ungroup()
        
    }
    
    return(ds)
}

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
get.phase.diagram.temp.nacl <- function(phase.diagram.ds, system.properties) {
    precision <- 1e-5
    
    ds <- phase.diagram.ds
    ds <- 
        as.data.frame(ds) %>% 
        mutate(conc.polymer = conc.p * system.properties$MW[1] + conc.q * system.properties$MW[2])
    
    ds2 <- lapply(unique(ds$tempC), function(tempc) {
        ds.t1 <- ds %>% filter(tempC == tempc, phase == 'dilute')
        
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



g <- ggplot(ds %>% filter(conc.salt < 0.1, tempC == 25), aes(x = conc.p, y = conc.salt, group = tempC)) +
    geom_point(aes(col = tempC))
print(g)

ds <- get.phase.diagram(system.properties, fitting.para)
ds2 <- get.phase.diagram.temp.conc(ds, system.properties)
ds3 <- get.phase.diagram.temp.nacl(ds, system.properties)

print(ds2)
g2 <- ggplot(ds2, aes(x = conc.polymer, y = tempC)) + geom_point()
print(g2)
g3 <- ggplot(ds3, aes(x = conc.salt, y = tempC)) + geom_point()
print(g3)

