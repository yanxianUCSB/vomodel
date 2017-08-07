# Analytical
SysProp <- R6Class(
    'Binodal System',
    private = list(),
    public = list(
        initialize = function(){},
        phi.polymer = NULL,
        phi.salt = NULL,
        temp = NULL,
        alpha = NULL,
        sigma = NULL,
        lattice.size = NULL,
        phi.polymer.pair = NULL
    )
)
Gibbs.Funs <- R6Class(
    'Gibbs.Funs',
    private = list(
        system.properties = list(),
        pkpq = function(..., sysprop = NULL) {
            # kpq = phi_q / pho_p
            # arg <- ifelse(missing(arg), list(...), arg)
            if (is.null(sysprop))
                arg <- list(...)
            else
                arg <- sysprop
            
            return( (arg$size.ratio[2] * arg$sigma[2] * arg$polymer.num[2] * arg$molar.ratio[2]) /
                        (arg$size.ratio[1] * arg$sigma[1] * arg$polymer.num[1] * arg$molar.ratio[1]) )
        },
        pS = function(sigma.p, sigma.q, kpq) {
            return((sigma.p + sigma.q * kpq) / (1 + kpq))
        },
        pA = function(r.p, r.q, kpq, size.ratio.p, size.ratio.q) {
            return(
                (1 / (r.p*size.ratio.p) + 
                     kpq / (r.q*size.ratio.q)) * (1 / (1 + kpq))
            )
        },
        pB = function(r.p, r.q, kpq, size.ratio.p, size.ratio.q) {
            return(
                -log(1 + kpq) / (r.p * size.ratio.p * (1 + kpq)) -
                    kpq * (log(1 + kpq) - log(kpq)) / ((r.q * size.ratio.q) * (1 + kpq))
            )
        },
        pp = function(sysprop) {
            arg <- sysprop
            
            kpq <- (arg$size.ratio[2] * arg$sigma[2] * arg$polymer.num[2] * arg$molar.ratio[2]) /
                (arg$size.ratio[1] * arg$sigma[1] * arg$polymer.num[1] * arg$molar.ratio[1])
            
            sigma.p <- arg$sigma[1]
            sigma.q <- arg$sigma[2]
            r.p <- arg$polymer.num[1]
            r.q <- arg$polymer.num[2]
            size.ratio.p <- arg$size.ratio[1]
            size.ratio.q <- arg$size.ratio[2]
            
            out <- list()
            out$S <- (sigma.p + sigma.q * kpq) / (1 + kpq)
            out$A <- (1 / (r.p*size.ratio.p) + kpq / (r.q*size.ratio.q)) / (1 + kpq)
            out$B <- 
                -log(1 + kpq) / (r.p * size.ratio.p * (1 + kpq)) -
                kpq * (log(1 + kpq) - log(kpq)) / ((r.q * size.ratio.q) * (1 + kpq))
            out$kpq <- kpq
            
            return(out)
        },
        get.alpha = function(temp, size) {
            # units: K, m
            alpha <-
                2 / 3 * sqrt(pi) * ((ke ^ 2 / (kEr * kkB * temp)) / size) ^ (3 / 2)
        },
        alpha.dT = function(temp, size) {
            lB <- ke ^ 2 / (kEr * kkB * temp)
            1.5 * sqrt(pi) * 1.5 * sqrt(lB / size) * 1 / size * lB * (-1 / temp)
        },
        gibbs = function(phi.polymer, phi.salt, alpha, Chi, para, temp) {
            #' ASSUMPTIONS
            #' 1. the solutions are in lattice with lattice size equal to the monomer size and ion sizes.
            #' 2. lattice size, monomer size, ion and water size are equal to 0.31 nm.
            #' 3. 
            # ... = alpha, sigma, polymer.num, kpq
            # return F/NkT
            kpq <- para$kpq
            phi <- c(phi.polymer/(1+kpq), phi.polymer*kpq/(1+kpq), phi.salt*0.5, phi.salt*0.5, 1-phi.salt-phi.polymer)
            return(
                (k.vol * kkB * temp) / k.water.size ^ 3 * (
                    -alpha * (para$S * phi.polymer + phi.salt) ^ 1.5 +
                        para$A * phi.polymer * log(phi.polymer) +
                        para$B * phi.polymer +
                        phi.salt * log(0.5 * phi.salt) +
                        (1 - phi.polymer - phi.salt) * log(1 - phi.polymer - phi.salt) +
                        t(phi) %*% Chi %*% phi
                )
            )
        },
        gibbs.d = function(phi.polymer, phi.salt, alpha, Chi, para, temp) {
            # ... = alpha, sigma, polymer.num, kpq
            kpq <- para$kpq
            phi <- c(phi.polymer/(1+kpq), phi.polymer*kpq/(1+kpq), phi.salt*0.5, phi.salt*0.5, 1-phi.salt-phi.polymer)
            return(
                (k.vol * kkB * temp) / k.water.size ^ 3 * (
                    -alfa * 1.5 * (para$S * phi.polymer + phi.salt) ^ 0.5 * para$S +
                        para$A * log(phi.polymer) + para$A +
                        para$B -
                        log(1 - phi.polymer - phi.salt) - 1 +
                        2*phi.polymer*(1/(1+kpq)^2)*(Chi[1,1]+Chi[1,2]+Chi[2,1]+Chi[2,2]) +
                        (1/(1+kpq))*((Chi[1,3]+Chi[3,1]+Chi[2,3]+Chi[3,2])*phi[3] + (Chi[1,4]+Chi[4,1]+Chi[2,4]+Chi[4,2])*phi[4] + (Chi[1,5]+Chi[5,1]+Chi[2,5]+Chi[5,2])*phi[5])
                    # 2*Chi[1,1]*kpq/(1+kpq)^2 * phi.polymer + 2/(1+kpq)*Chi[1,1:5]%*%phi - 2/(1+kpq)*Chi[1,1]*phi[1]
                )
            )
        },
        gibbs.dd = function(phi.polymer, phi.salt, alpha, Chi, para, temp) {
            # ... = alpha, sigma, polymer.num, kpq
            kpq <- para$kpq
            phi <- c(phi.polymer/(1+kpq), phi.polymer*kpq/(1+kpq), phi.salt*0.5, phi.salt*0.5, 1-phi.salt-phi.polymer)
            return(
                (k.vol * kkB * temp) / k.water.size ^ 3 * (
                    -alpha * 0.75 * (para$S * phi.polymer + phi.salt) ^ (-0.5) * para$S ^ 2 +
                        para$A * 1 / phi.polymer +
                        1 / (1 - phi.polymer - phi.salt) +
                        2*(1/(1+kpq)^2)*(Chi[1,1]+Chi[1,2]+Chi[2,1]+Chi[2,2])
                    # 2*arg$Chi[1,1]*kpq/(1+kpq)^2
                )
            )
        },
        gibbs.ddd = function(phi.polymer, phi.salt, alpha, Chi, para, temp) {
            # ... = alpha, sigma, polymer.num, kpq
            kpq <- para$kpq
            phi <- c(phi.polymer/(1+kpq), phi.polymer*kpq/(1+kpq), phi.salt*0.5, phi.salt*0.5, 1-phi.salt-phi.polymer)
            return(
                (k.vol * kkB * temp) / k.water.size ^ 3 * (
                    alpha * 0.375 * (para$S * phi.polymer + phi.salt) ^ (-1.5) * para$S ^ 3 -
                        para$A * 1 / phi.polymer ^ 2 +
                        1 / (1 - phi.polymer - phi.salt) ^ 2
                )
            )
        },
        gibbs.pfps = function(phi.polymer, phi.salt, alpha, Chi, para, temp) {
            # -para.alpha * 1.5 * (para.S * phip + phis) ** 0.5 +
            # log(0.5 * phis) + 1 +
            # -log(1 - phip - phis) - 1
            kpq <- para$kpq
            phi <- c(phi.polymer/(1+kpq), phi.polymer*kpq/(1+kpq), phi.salt*0.5, phi.salt*0.5, 1-phi.salt-phi.polymer)
            return((k.vol * kkB * temp) / k.water.size ^ 3 * (
                -alpha * 1.5 * (para$S * phi.polymer + phi.salt) ^ 0.5 +
                    log(0.5 * phi.salt) + 1 +
                    -log(1 - phi.polymer - phi.salt) - 1 +
                    2*phi.salt*(1/4)*(Chi[3,3]+Chi[3,4]+Chi[4,3]+Chi[4,4]) +
                    (1/2)*((Chi[1,3]+Chi[3,1]+Chi[4,1]+Chi[1,4])*phi[1] + (Chi[3,2]+Chi[2,3]+Chi[4,2]+Chi[2,4])*phi[2] + (Chi[3,5]+Chi[5,3]+Chi[4,5]+Chi[5,4])*phi[5])
                # (Chi[1,3]+Chi[1,4]-2*Chi[1,5])/(1+kpq) * phi.polymer
            ))
        },
        gibbs.pdfps = function(phi.polymer, phi.salt, alpha, Chi, para, temp) {
            #         -para.alpha * 0.75 * (para.S * phip + phis) ** (-0.5) * para.S +
            # 1/(1 - phip - phis)
            kpq <- para$kpq
            phi <- c(phi.polymer/(1+kpq), phi.polymer*kpq/(1+kpq), phi.salt*0.5, phi.salt*0.5, 1-phi.salt-phi.polymer)
            return((k.vol * kkB * temp) / k.water.size ^ 3 * (
                -alpha * 0.75 * (para$S * phi.polymer + phi.salt) ^ (-0.5) * para$S +
                    1 / (1 - phi.polymer - phi.salt) +
                    (1/(1+kpq))*((Chi[1,3]+Chi[3,1]+Chi[2,3]+Chi[3,2])*(1/2) + (Chi[1,4]+Chi[4,1]+Chi[2,4]+Chi[4,2])*(1/2) - (Chi[1,5]+Chi[5,1]+Chi[2,5]+Chi[5,2]))
                # (Chi[1,3]+Chi[1,4]-2*Chi[1,5])/(1+kpq)
            ))
        }
    ),
    public = list(
        initialize = function() {
        },
        get.gibbs = function() {
            sysprop <-  self$sysprop
            
            gibbs <- private$gibbs
            gibbs.d <- private$gibbs.d
            gibbs.dd <- private$gibbs.dd
            gibbs.ddd <- private$gibbs.ddd
            
            phi.polymer <- sysprop$phi.polymer.pair
            phi.salt <- sysprop$phi.salt
            temp <- sysprop$temp
            Chi <- sysprop$Chi
            
            para <- private$pp(sysprop)
            alpha <- private$get.alpha(temp, sysprop$lattice.size)
            
            out <-
                expand.grid(phi.polymer = phi.polymer, phi.salt = phi.salt) %>%
                group_by(phi.polymer) %>% 
                mutate(
                    f = gibbs(phi.polymer, phi.salt, alpha, Chi, para, temp),
                    df = gibbs.d(phi.polymer, phi.salt, alpha, Chi, para, temp),
                    ddf = gibbs.dd(phi.polymer, phi.salt, alpha, Chi, para, temp),
                    dddf = gibbs.ddd(phi.polymer, phi.salt, alpha, Chi, para, temp)
                ) %>% 
                ungroup()
            return(out)
        },
        get.binodal.curve.fun = function(phi.polymer.pair) {
            sysprop = self$sysprop
            
            gibbs <- private$gibbs
            gibbs.d <- private$gibbs.d
            phi.polymer <- phi.polymer.pair
            phi.salt <- sysprop$phi.salt
            temp <- sysprop$temp
            Chi <- sysprop$Chi
            
            para <- private$pp(sysprop)
            alpha <- private$get.alpha(temp, sysprop$lattice.size)
            
            return(
                c(
                    gibbs.d(phi.polymer[1], phi.salt, alpha, Chi, para, temp) - 
                        gibbs.d(phi.polymer[2], phi.salt, alpha, Chi, para, temp),
                    (phi.polymer[2] - phi.polymer[1]) * gibbs.d(phi.polymer[1], phi.salt, alpha, Chi, para, temp) - 
                        (gibbs(phi.polymer[2], phi.salt, alpha, Chi, para, temp) - 
                             gibbs(phi.polymer[1], phi.salt, alpha, Chi, para, temp))
                )
            )
        }
        # binodal.curve.jacobian_ = function(x, phi.polymer.2, sysprop){
        #     phi.polymer.1 <- x[1]
        #     phi.salt <- x[2]
        #     
        # },
    )
)
Binodal.Generator <- R6Class(
    'Binodal.Generator',
    private = list(
        binodal.curve.jacobian.phi.polymer.phi.salt = function(x, phi.polymer.2) {
            sysprop <-  self$sysprop
            
            phi.polymer.1 <- x[1]
            phi.polymer.2 <- phi.polymer.2
            phi.salt <- x[2]
            
            gibbs <- private$gibbs
            gibbs.d <- private$gibbs.d
            
            temp <- sysprop$temp
            Chi <- sysprop$Chi
            
            para <- private$pp(sysprop)
            alpha <- private$get.alpha(temp, sysprop$lattice.size)
            
            return(matrix(
                c(
                    gibbs.dd(phi.polymer.1, phi.salt, alpha, Chi, para, temp),
                    
                    gibbs.d(phi.polymer.1, phi.salt, alpha, Chi, para, temp) - gibbs.d(phi.polymer.2, phi.salt, alpha, Chi, para, temp),
                    
                    gibbs.pdfps(phi.polymer.1, phi.salt, alpha, Chi, para, temp) - gibbs.pdfps(phi.polymer.2, phi.salt, alpha, Chi, para, temp),
                    
                    (phi.polymer.2 - phi.polymer.1) * gibbs.pdfps(phi.polymer.2, phi.salt, alpha, Chi, para, temp) -
                        (
                            gibbs.pfps(phi.polymer.2, phi.salt, alpha, Chi, para, temp) - gibbs.pfps(phi.polymer.1, phi.salt, alpha, Chi, para, temp)
                        )
                ),
                nrow = 2,
                ncol = 2
            ))
        }
        
    ),
    public = list(
        sysprop = NULL,
        gibbs.funs = NULL,
        get.binodal.curve = function( sysprop, fitting.para, critical.generator ) {
            #' generate binodal curve
            arg <- sysprop
            
            c.point <- critical.generator$get.critical.point()
            
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
        },
        get.temp.conc.curve = function(){},
        get.temp.nacl.curve = function(){}
    )
)
Critical.Generator <- R6Class(
    'Critical.Generator',
    private = list(
        ## Critical point
        c.point.temp = function(sysprop, fitting.para) {
            temp <- seq(273.15, 333.15, 1)
            do.call(rbind, lapply(temp, function(temp) {
                alpha <- get.alpha(temp, sysprop$water.size)
                polymer.num <- sysprop$polymer.num
                sigma <- sysprop$sigma
                size.ratio <- sysprop$size.ratio
                Chi <- sysprop$Chi
                molar.ratio <- sysprop$molar.ratio
                ds <- critical.point_(guess = fitting.para$critical.point.guess, temp = temp, polymer.num = polymer.num,
                                      alpha = alpha, sigma = sigma, size.ratio = size.ratio, Chi = Chi, molar.ratio = molar.ratio)
                # if (DEBUG) print(ds)
                ds$temp <-  temp
                return(as.data.frame(ds))
            }))
        },
        c.point.temp.fun = function(c.point.temp.ds) {
            return(function(temp) {
                list(
                    phi.polymer = spline(c.point.temp.ds$temp, c.point.temp.ds$phi.polymer, xout = temp)$y,
                    phi.salt = spline(c.point.temp.ds$temp, c.point.temp.ds$phi.salt, xout = temp)$y)
            })
        },
        critical.point.fun_ = function(phis, ...) {
            arg <- list(...)
            return(c(
                gibbs.dd(phi.polymer = phis[1], phi.salt = phis[2], ...),
                gibbs.ddd(phi.polymer = phis[1], phi.salt = phis[2], ...)
            ))
        },
        critical.point_ = function(guess, ...) {
            # output list of critical points
            arg <- list(...)
            guess <- critical.point.get.guess(guess, ...)
            out <- nleqslv(
                x = guess,
                fn = critical.point.fun_,
                ...)
            return(list(phi.polymer = out$x[1], phi.salt = out$x[2]))
        },
        critical.point.get.guess = function(guess, ...) {
            for(i in 1:10) {
                test <- critical.point.fun_(guess, ...)
                if(anyNA(test) ||
                   any(is.nan(test)) ||
                   abs(test[1]) > 1e6 ||
                   abs(test[2]) > 1e6
                ) {
                    guess[1] <- guess[1]
                    guess[2] <- guess[2] * 1.01
                    next
                } else 
                    break
            }
            return(guess)
        }
    ),
    public = list(
        initialize = function(sysprop, fitting.para){
            self$sysprop <- sysprop
            self$fitting.para <- fitting.para
        },
        sysprop = NULL,
        fitting.para = NULL,
        get.critical.point = function(){
            assertthat::assert_that(c.point$phi.polymer > 0)
            
        }
    )
)

