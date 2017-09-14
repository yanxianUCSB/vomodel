# Voorn-Overbeek Modeling
rm(list = ls())
kNa                         <<- 6.02E23
# kkB                         <<- 1.38064852E-23         # Boltzmann Constant
kkB                         <<- 1.38E-23         # Boltzmann Constant
# ke                          <<- 1.60217662E-19         # Elementary charge
ke                          <<- 1.60E-19         # Elementary charge
kEr                         <<- 4 * pi * 80 * 8.85E-12 # F/m; vaccuum permitivity = 8.85E-12

# k.water.size                <<- 0.31E-9     # calculated from 18 cm^3 / mol
k.water.size                <<- (18E-6 / kNa)^(1/3)     # calculated from 18 cm^3 / mol
k.na.size                   <<- 0.235E-9    # (Marcus 1988)
k.cl.size                   <<- 0.318E-9    # (Marcus 1988)
k.water.conc                <<- 1000 / 18.01528 * 1000  # water concentration

k.vol                       <<- 120E-9

k.binodal.guess.offset      <<- 1.05

pkpq                       <- function(..., sysprop = NULL) {
    # kpq = phi_q / pho_p
    # arg <- ifelse(missing(arg), list(...), arg)
    if (is.null(sysprop))
        arg <- list(...)
    else
        arg <- sysprop
    size.1 <- arg$size.ratio[1]
    size.2 <- arg$size.ratio[2]
    sigma.1 <- arg$sigma[1]
    sigma.2 <- arg$sigma[2]
    M.1 <- arg$polymer.num[1]
    M.2 <- arg$polymer.num[2]
    N.1 <- arg$molar.ratio[1]
    N.2 <- arg$molar.ratio[2]
    
    out1 <- (N.2 * M.2 * size.2^3 ) / (N.1 * M.1 * size.1^3 )
    
    # out2 <- (N.2 * M.2 * size.2 * sigma.2) / (N.1 * M.1 * size.1 * sigma.1)
    
    return( out1 )
}
pS                         <- function(sigma.p, sigma.q, kpq) {
    return((sigma.p + sigma.q * kpq) / (1 + kpq))
}
pA                         <- function(r.p, r.q, kpq, size.ratio.p, size.ratio.q) {
    return(
        (1 / (r.p*size.ratio.p) + 
             kpq / (r.q*size.ratio.q)) * (1 / (1 + kpq))
    )
}
pB                         <- function(r.p, r.q, kpq, size.ratio.p, size.ratio.q) {
    return(
        -log(1 + kpq) / (r.p * size.ratio.p * (1 + kpq)) -
            kpq * (log(1 + kpq) - log(kpq)) / ((r.q * size.ratio.q) * (1 + kpq))
    )
}
pp                         <- function(..., sysprop = NULL) {
    if (is.null(sysprop)) {
        arg <- list(...)
        
    } else {
        arg <- sysprop
    }
    out <- list()
    kpq <- pkpq(sysprop = arg)
    out$kpq <- kpq
    out$S <- pS(arg$sigma[1], arg$sigma[2], kpq)
    out$A <- pA(arg$polymer.num[1], arg$polymer.num[2], kpq, arg$size.ratio[1], arg$size.ratio[2])
    out$B <- pB(arg$polymer.num[1], arg$polymer.num[2], kpq, arg$size.ratio[1], arg$size.ratio[2])
    return(out)
}
get.alpha                   <- function(temp, size) {
    # units: K, m
    alpha <-
        2 / 3 * sqrt(pi) * ((ke ^ 2 / (kEr * kkB * temp)) / size) ^ (3 / 2)
}
get.temp                    <- function(alpha, size) {
    a <- 
        alpha / (2/3*sqrt(pi)) ^ (-3/2) * size
    temp <- ke^2/(kEr*kkB) / a
}

gibbs                      <- function(phi.polymer, phi.salt, ...) {
    #' ASSUMPTIONS
    #' 1. the solutions are in lattice with lattice size equal to the monomer size and ion sizes.
    #' 2. lattice size, monomer size, ion and water size are equal to 0.31 nm.
    #' 3. 
    # ... = alpha, sigma, polymer.num, kpq
    # return F/NkT
    arg <- list(...)
    alpha <- arg$alpha
    para <- pp(sysprop = arg)
    kpq <- para$kpq
    phi <- c(phi.polymer/(1+kpq), phi.polymer*kpq/(1+kpq), phi.salt*0.5, phi.salt*0.5, 1-phi.salt-phi.polymer)
    # temp <- get.temp(alpha, k.water.size)
    return(
        # (k.vol * kkB * temp) / k.water.size ^ 3 * (
        -alpha * (para$S * phi.polymer + phi.salt) ^ 1.5 +
            para$A * phi.polymer * log(phi.polymer) +
            para$B * phi.polymer +
            phi.salt * log(0.5 * phi.salt) +
            (1 - phi.polymer - phi.salt) * log(1 - phi.polymer - phi.salt) +
            t(phi) %*% arg$Chi %*% phi
        # )
    )
}
gibbs.d                    <- function(phi.polymer, phi.salt, ...) {
    # ... = alpha, sigma, polymer.num, kpq
    arg <- list(...)
    alfa <- arg$alpha
    Chi <- arg$Chi
    para <- pp(sysprop = arg)
    kpq <- para$kpq
    phi <- c(phi.polymer/(1+kpq), phi.polymer*kpq/(1+kpq), phi.salt*0.5, phi.salt*0.5, 1-phi.salt-phi.polymer)
    # temp <- get.temp(alpha, k.water.size)
    return(numDeriv::grad(func = gibbs, x = phi.polymer, method = 'simple', phi.salt = phi.salt, ...))
    return(
        # (k.vol * kkB * temp) / k.water.size ^ 3 * (
        -alfa * 1.5 * (para$S * phi.polymer + phi.salt) ^ 0.5 * para$S +
            para$A * log(phi.polymer) + para$A +
            para$B -
            log(1 - phi.polymer - phi.salt) - 1 +
            # 4 * Chi[1, 2] * kpq / (1+kpq)^2 * phi.polymer
            # 2 * Chi[1, 5] * 1   / (1+kpq)   * phi[5]
            2*phi.polymer*(1/(1+kpq)^2) * (Chi[1,1] + kpq*Chi[1,2] + kpq*Chi[2,1] + kpq^2*Chi[2,2]) +
              1/(1+kpq) * 2*(phi[3]*Chi[1,3] + phi[4]*Chi[1,4]) +
            kpq/(1+kpq) * 2*(phi[3]*Chi[2,3] + phi[4]*Chi[2,4]) +
            2 * 1/(1+kpq) * (1-phi.polymer-phi.salt) * Chi[1,5] - 2 * 1/(1+kpq) * phi.polymer * Chi[1,5] +
            2 * kpq/(1+kpq) * (1-phi.polymer-phi.salt) * Chi[2,5] - 2 * kpq/(1+kpq) * phi.polymer * Chi[2,5] +
            -2 * (phi[3] + phi[4]) +
            -4 * (1-phi.polymer-phi.salt) * Chi[5,5]
            
            # (1/(1+kpq))*(( Chi[1,3] + Chi[3,1] + kpq*Chi[2,3] + kpq*Chi[3,2])*phi[3] +
            #                  (Chi[1,4] + Chi[4,1] + kpq*Chi[2,4] + kpq*Chi[4,2])*phi[4] +
            #                  (Chi[1,5] + Chi[5,1] + kpq*Chi[2,5] + kpq*Chi[5,2])*phi[5])
        # 2*Chi[1,1]*kpq/(1+kpq)^2 * phi.polymer + 2/(1+kpq)*Chi[1,1:5]%*%phi - 2/(1+kpq)*Chi[1,1]*phi[1]
        # )
    )
}
# gibbs.dd                   <- function(phi.polymer, phi.salt, ...) {
#     # ... = alpha, sigma, polymer.num, kpq
#     arg <- list(...)
#     Chi <- arg$Chi
#     alpha <- arg$alpha
#     para <- pp(sysprop = arg)
#     kpq <- para$kpq
#     phi <- c(phi.polymer/(1+kpq), phi.polymer*kpq/(1+kpq), phi.salt*0.5, phi.salt*0.5, 1-phi.salt-phi.polymer)
#     return(
#         # (k.vol * kkB * arg$temp) / k.water.size ^ 3 * (
#             -alpha * 0.75 * (para$S * phi.polymer + phi.salt) ^ (-0.5) * para$S ^ 2 +
#                 para$A * 1 / phi.polymer +
#                 1 / (1 - phi.polymer - phi.salt) +
#                 2*(1/(1+kpq)^2)*(Chi[1,1]+Chi[1,2]+Chi[2,1]+Chi[2,2])
#             # 2*arg$Chi[1,1]*kpq/(1+kpq)^2
#         # )
#     )
# }
# gibbs.ddd                  <- function(phi.polymer, phi.salt, ...) {
#     # ... = alpha, sigma, polymer.num, kpq
#     arg <- list(...)
#     Chi <- arg$Chi
#     alpha <- arg$alpha
#     para <- pp(sysprop = arg)
#     kpq <- para$kpq
#     phi <- c(phi.polymer/(1+kpq), phi.polymer*kpq/(1+kpq), phi.salt*0.5, phi.salt*0.5, 1-phi.salt-phi.polymer)
#     return(
#         # (k.vol * kkB * arg$temp) / k.water.size ^ 3 * (
#             alpha * 0.375 * (para$S * phi.polymer + phi.salt) ^ (-1.5) * para$S ^ 3 -
#                 para$A * 1 / phi.polymer ^ 2 +
#                 1 / (1 - phi.polymer - phi.salt) ^ 2
#         # )
#     )
# }
# 
# critical.point_            <- function(guess,  sysprop, fitting.para, default = NULL) {
#     # output list of critical points
#     
#     i <- critical.point.fun_(guess, 
#                              alpha = sysprop$alpha,
#                              sigma = sysprop$sigma,
#                              Chi = sysprop$Chi,
#                              polymer.num = sysprop$polymer.num,
#                              size.ratio = sysprop$size.ratio,
#                              molar.ratio = sysprop$molar.ratio)
#     if (any(is.nan(i))) return(default)
#     
#     out <- nleqslv(
#         x = guess,
#         fn = critical.point.fun_,
#         # jac = critical.point.jac,
#         control = list(allowSingular = T, xtol = 1e-10),
#         global = 'cline',
#         method = 'Newton',
#         alpha = sysprop$alpha,
#         sigma = sysprop$sigma,
#         Chi = sysprop$Chi,
#         polymer.num = sysprop$polymer.num,
#         size.ratio = sysprop$size.ratio,
#         molar.ratio = sysprop$molar.ratio)
#     
#     result <- list(phi.polymer = out$x[1], phi.salt = out$x[2])
#     
#     if (result$phi.polymer < 0) return(default)
#     
#     return(result)
# }
# critical.point.fun_        <- function(phis, ...) {
#     arg <- list(...)
#     ddg <- gibbs.dd(phi.polymer = phis[1], phi.salt = phis[2], ...)
#     dddg <- gibbs.ddd(phi.polymer = phis[1], phi.salt = phis[2], ...)
#     return(c(
#         ddg,
#         dddg
#     ))
# }
binodal.curve.fun_              <- function(x, phi.polymer.2, ...) {
    
    phi.polymer.1 <- x[1]
    phi.salt <- x[2]
    phi.polymer <- c(phi.polymer.1, phi.polymer.2)
    
    phi1 <- phi.polymer[1]
    phi2 <- phi.polymer[2]
    gphi1 <- gibbs(phi.polymer[1], phi.salt, ...)
    gphi2 <- gibbs(phi.polymer[2], phi.salt, ...)
    dgdphi1 <- gibbs.d(phi.polymer[1], phi.salt, ...)
    dgdphi2 <- gibbs.d(phi.polymer[2], phi.salt, ...)
    out <- c(
        dgdphi1 - dgdphi2,
        (phi2 - phi1) * dgdphi1 - (gphi2 - gphi1)
    )
    
    return(out)
}
binodal.curve_                  <- function( sysprop = NULL, fitting.para = NULL) {
    #' generate binodal curve
    arg <- sysprop
    
    c.point <- fitting.para$default.critical.point
    
    n <- floor(0.5 * (1 + sqrt(1 + 8 * fitting.para$sampling.end / fitting.para$sampling.gap)))
    phi2.seq <-  fitting.para$sampling.gap * seq(1, n, 1) * seq(2, n + 1, 1) * 0.5
    phi2.seq <-  phi2.seq[which(phi2.seq >= fitting.para$sampling.start)]
    
    
    # Binodal Guess
    binodal.guess0 <- fitting.para$binodal.guess
    # binodal.guess0[2] <- c.point$phi.salt
    
    for (i in seq(1, 4, 1)) {
        
        # i <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)[i]
        i <- 10^i
        
        binodal.guess <- binodal.guess0*i
        cat('current binodal guess: ')
        cat(i)
        cat('/10000\n')
        
        guess.i <- 1
        phi.salt.old <- 0
        phi.1.old <- 1
        
        output <- list()
        
        for (phi.polymer.2 in phi2.seq) {
            
            # protect nleqslv from receiving NaN
            test <- binodal.curve.fun_(x = binodal.guess, phi.polymer.2 = phi.polymer.2,
                                       alpha = sysprop$alpha,
                                       sigma = sysprop$sigma,
                                       Chi = sysprop$Chi,
                                       polymer.num = sysprop$polymer.num,
                                       size.ratio = sysprop$size.ratio,
                                       molar.ratio = sysprop$molar.ratio)
            if (anyNA(test) || is.nan(test[1]) || is.nan(test[2])) {
                next
            }
            # nleqslv solving binodal curve points
            roots <- nleqslv (
                binodal.guess,
                binodal.curve.fun_,
                # binodal.curve.jacobian_,
                phi.polymer.2 = phi.polymer.2, 
                control = list(xtol = fitting.para$epsilon),
                method = 'Newton',
                global = 'cline',
                alpha = sysprop$alpha,
                sigma = sysprop$sigma,
                Chi = sysprop$Chi,
                polymer.num = sysprop$polymer.num,
                size.ratio = sysprop$size.ratio,
                molar.ratio = sysprop$molar.ratio
                )
            # print(sysprop)
            # stop()
            # filter 
            # if (roots$x[2] < phi.salt.old || 
            #     roots$x[1] > phi.1.old
            #     ) break()
            
            if (
                roots$x[2] > phi.salt.old &&
                roots$x[1] < phi.1.old &&
                roots$x[1] > 0 &&
                roots$x[1] > phi.polymer.2 + 1e-8 &&
                roots$x[2] > 0 &&
                # roots$x[2] < c.point$phi.salt &&
                max(abs(roots$fvec)) < 1 &&
                roots$termcd %in% c(1)
                ) {
                
                
                binodal.guess <- roots$x
                # print(phi.polymer.2)
                # print(roots)
                
                output[[length(output) + 1]] <- 
                    c(
                        phi.polymer.1 = roots$x[1],
                        phi.polymer.2 = phi.polymer.2,
                        phi.salt = roots$x[2],
                        f1 = roots$fvec[1],
                        f2 = roots$fval[2],
                        pair = phi.polymer.2
                    )
                
                guess.i <- 1
                phi.salt.old <- roots$x[2]
                phi.1.old <- roots$x[1]
                
            # modify guess and re compute
            } else {
                # print(phi.polymer.2)
                # print(c.point)
                # print(roots$termcd)
                # print(binodal.guess)
                # if(roots$termcd == 1) print(roots)
                # if(DEBUG) k.binodal.guess.offset <- 1.05
                if (guess.i == 1) 
                    binodal.guess[2] <- binodal.guess[2] * 1.05
                else 
                    binodal.guess[2] <- binodal.guess[2] * 1.05
                guess.i <- guess.i + 1
                next
            }
        }
        
        if (length(output) > 4) {
            cat('binodal guess good: ')
            cat(binodal.guess)
            cat(' nleqslv completed\n')
            break
        }
    }
    
    if (length(output) < 1) return()
    # assertthat::assert_that(length(output) > 0, msg = 'binodal.curve_ failed')
    
    p2 <- as.data.frame.matrix(do.call(rbind, output)) %>%
        filter(phi.polymer.1 > max(phi.polymer.2))  # requiring the dense phase should be larger than dilute phase
    
    if (nrow(p2) < 1) return()
    # assertthat::assert_that(nrow(p2) > 0)
    
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
    
    cat(' nleqslv makes sense\n')
    return(ds)
}
get.binodal.curve           <- function(tempC, 
                                        sysprop, fitting.para, 
                                        condensation = T,
                                        counterion.release = T,
                                        unit = 'mol') {
    #' get binodal curve by Voorn Overbeek Model 3D component system
    #' temp in degree C
    #' Chi = 0 or 5 * 5 matrix
    temp <- tempC + 273.15
    
    sysprop$alpha <- get.alpha(temp, k.water.size)
    cat('alpha: ')
    cat(sysprop$alpha)
    cat('\n')
    
    
    # Condensation
    if (condensation) {
        lB <- 0.7E-9  # Bjerrum length at 298K
        # sysprop$sigma[2] <- sysprop$size.ratio[2]*k.water.size / lB
        sysprop$sigma[2] <- k.water.size / lB
        cat('condensation: current sigma[2] = ')
        cat(sysprop$sigma[2])
        cat('\n')
    } 
    if (counterion.release) {
        # update Bjerrum length and the effective charge density of RNA
        lB <- ke^2 / (kEr*kkB*(tempC+273.15))
        # sysprop$sigma[2] <- sysprop$size.ratio[2]*k.water.size / lB
        sysprop$sigma[2] <- k.water.size / lB
        cat('counterion release: current sigma[2] = ')
        cat(sysprop$sigma[2])
        cat('\n')
    } 
    
    
    # # Chi
    Chi298 <- sysprop$Chi
    sysprop$Chi <- Chi298 * 298 / temp
    # cat('Chipq, Chipw: ')
    # cat(sysprop$Chi[1,2])
    # cat(' ')
    # cat(sysprop$Chi[1,5])
    # cat('\n')
    
    kpq <- pkpq(sysprop = sysprop)
    
    ds <- binodal.curve_(sysprop = sysprop, fitting.para = fitting.para) 
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
                        (sysprop$size.ratio[3] *
                        sysprop$water.size) ^ 3
                ) / 1000,
                # mol/L
                conc.p = phi.polymer / (1 + kpq) / (
                    kNa * sysprop$polymer.num[1] *
                        (sysprop$size.ratio[1] *
                        sysprop$water.size) ^ 3
                ) / 1000,
                conc.q = phi.polymer * kpq / (1 + kpq) / (
                    kNa * sysprop$polymer.num[2] *
                        (sysprop$size.ratio[2] *
                        sysprop$water.size) ^ 3
                ) / 1000
            ) %>%
            ungroup() %>% 
            mutate(
                conc.mass.polymer = conc.p * sysprop$MW[1] + conc.q * sysprop$MW[2]
            )
    }
    return(ds)
}
get.phase.diagram           <- function(system.properties, fitting.para, temp.range = NULL) {
    

    
    p <- lapply(temp.range, function(tempC) {
        
        
        # update critical.point.guess
        # fitting.para$critical.point.guess <- as.numeric(fitting.para$c.point.temp.fun(tempC + 273))
        
            cat(paste0('Cur. Temp [C]: ', tempC, '\n'))
            cat('fitting >>>\n')
            
            out <- get.binodal.curve(tempC, sysprop = system.properties, 
                                   fitting.para = fitting.para, 
                                   condensation = fitting.para$condensation, 
                                   counterion.release = fitting.para$counterion.release)
        
        if (is.null(out) || nrow(out) < 2) {
            cat('get.binodal.curve failed!\n\n')
            return(NULL)
        } else {
            cat('get.binodal.curve succeeded!\n\n')
            return(out)
        }
    })
    
    out <- do.call(rbind, p) 
    
    if (is.null(out)) return()
    
    out <- out %>% 
        filter(phase.separated) %>% 
        mutate(conc.polymer = conc.p * system.properties$MW[1] + conc.q * system.properties$MW[2])
    
    return(out) 
}

get.phase.diagram.temp.conc <- function(phase.diagram.ds, system.properties, k.conc.salt = 0.030) {
    if (is.null(phase.diagram.ds)) return()
    # //TODO:
    # precision <- 1e-5
    
    ds <- as.data.frame(phase.diagram.ds)
    
    ds2 <- lapply(unique(ds$tempC), function(tempc) {
        ds.t1 <- ds %>% filter(tempC == tempc, phase == 'dilute')
        ds.t2 <- ds %>% filter(tempC == tempc, phase == 'dense')
        
        if(k.conc.salt > max(c(ds.t1$conc.salt, ds.t2$conc.salt))) return()
        
        ds3 <- rbind(
            data.frame(conc.salt = k.conc.salt) %>%
                mutate(
                    conc.polymer = spline(ds.t1$conc.salt, ds.t1$conc.polymer, xout = conc.salt)$y,
                    conc.salt = k.conc.salt,
                    phase = 'dilute'
                ),
            data.frame(conc.salt = k.conc.salt) %>% 
                mutate(
                    conc.polymer = spline(ds.t2$conc.salt, ds.t2$conc.polymer, xout = conc.salt)$y,
                    conc.salt = k.conc.salt,
                    phase = 'dense'
                )
        ) %>%
            mutate(
                tempC = tempc
            ) 
        
        return(ds3)
    })
    
    ds3 <- do.call(rbind, ds2) 
    if(is.null(ds3)) return()
    ds4 <- ds3 %>% filter(conc.polymer > 0) %>% select(conc.polymer, tempC)
    return(ds4)
    
}
get.phase.diagram.temp.nacl <- function(phase.diagram.ds, system.properties, k.conc.polymer = 0.125) {
    if (is.null(phase.diagram.ds)) return()
    
    precision <- 1e-5
    
    ds <- phase.diagram.ds
    ds <- 
        as.data.frame(ds) %>% 
        mutate(conc.polymer = conc.p * system.properties$MW[1] + conc.q * system.properties$MW[2])
    
    ds2 <- lapply(unique(ds$tempC), function(tempc) {
        ds.t1 <- ds %>% filter(tempC == tempc, phase == 'dilute') 
        
        if (min(ds.t1$conc.salt) > 30) return()
        
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
        return(ds3)
    })
    
    ds3 <- do.call(rbind, ds2) 
    if (is.null(ds3)) return()
    
    ds4 <- ds3 %>% filter(conc.polymer < 20, conc.polymer > 0) %>% select(conc.salt, tempC)
    return(ds4)
    
}
get.phase.diagram.exp       <- function(dataset.file = 'dataset.csv') {
    dataset <- read.csv(dataset.file) %>% 
        mutate(conc.polymer = protein * 1E-6 * system.properties$MW[1] + rna * 1E-3) %>% 
        mutate(conc.salt = nacl * 1E-3) %>% 
        mutate(tempC.cp = cloudpoint,
               tempC.on = onset) %>% 
        select(conc.polymer, conc.salt, tempC.cp, tempC.on) %>% 
        group_by(conc.polymer, conc.salt) %>% 
        mutate(tempC.cpm = mean(tempC.cp),
               tempC.onm = mean(tempC.on)) %>% 
        ungroup()
    dataset <- dataset[!duplicated(dataset[c('conc.polymer', 'conc.salt')]), ]
    return(dataset)
}

