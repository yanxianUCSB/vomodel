# Voorn-Overbeek Modeling
kNa <<- 6.02E23
kkB <<- 1.38064852E-23
ke  <<- 1.60217662E-19
kEr <<- 4 * pi * 80 * 8.85E-12 # F/m; vaccuum permitivity = 8.85E-12
k.water.size <<- 0.31E-9  # 0.31nm as a size of a water molecule
k.vol <<- 120E-9
k.water.conc  <<- 1000 / 18.01528 * 1000  # water concentration
k.dna.contour.unit.length <<- 0.33E-9  # 0.33 nm


## USERS

get.binodal.curve <-
    function(tempC,
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
            temp = temp,
            Chi = Chi,
            alpha = get.alpha(temp, sysprop$water.size),
            sigma = sysprop$sigma,
            polymer.num = sysprop$polymer.num,
            size.ratio = sysprop$size.ratio,
            molar.ratio = sysprop$molar.ratio,
            guess.critical.point = fitting.para$critical.point.guess,
            binodal.guess = fitting.para$binodal.guess,
            epsilon = fitting.para$epsilon,
            sampling.gap = fitting.para$sampling.gap
        ) %>%
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
get.binodal.curves <-
    function(tempCs, Chi = 0, sysprop, fitting.para) {
        do.call(rbind, lapply(tempCs, function(tempC) {
            get.binodal.curve(tempC, Chi = 0, sysprop, fitting.para)
        }))
    }


## IMPLEMENTATION

get.alpha <- function(temp, size) {
    # units: K, m
    alpha <-
        2 / 3 * sqrt(pi) * ((ke ^ 2 / (kEr * kkB * temp)) / size) ^ (3 / 2)
}
alpha.dT <- function(temp, size) {
    lB <- ke ^ 2 / (kEr * kkB * temp)
    1.5 * sqrt(pi) * 1.5 * sqrt(lB / size) * 1 / size * lB * (-1 / temp)
}
Fen <- function(phis, rs, ws) {
    # units: 1, 1, 1
    Fen <- sum(phis * log(phis) / (ws * rs))
}
Fel <- function(alpha, sigma, phi) {
    # units: 1, 1, 1
    Fel <- 0 - alpha * sum(sigma * phi) ^ (1.5)
}
Fchi <- function(Chi, phi) {
    # units: 1, 1
    Fchi <- sum(Chi * (phi %*% t(phi)))
}
free.energy <-
    function(phi.polymer,
             phi.salt,
             alpha,
             sigma,
             Chi,
             temp,
             polymer.num,
             size.ratio) {
        # free.energy in Spruijt et al's work
        phi <-
            c(
                phi.polymer * 0.5,
                phi.polymer * 0.5,
                phi.salt * 0.5,
                phi.salt * 0.5,
                1 - phi.polymer - phi.salt
            )
        g <-
             (k.vol * kkB * temp)/k.water.size^3 *
            (Fen(phi, polymer.num, size.ratio) + Fel(alpha, sigma, phi) + Fchi(Chi, phi))
        return(g)
    }


# Numeircal Derivative

free.energy.f <- function(x, phi.salt, temp, ...) {
    l <- list(...)
    # free energy f
    free.energy(
        phi.polymer = x,
        phi.salt,
        alpha = l$alpha,
        sigma = l$sigma,
        Chi = l$Chi,
        temp,
        polymer.num = l$polymer.num,
        size.ratio = l$size.ratio
    )
}
free.energy.df <- function(x, phi.salt, temp, ...) {
    # df / dphi
    return(yx.numdiff(
        free.energy.f,
        x,
        phi.salt = phi.salt,
        temp = temp,
        ...
    ))
}
free.energy.ddf <- function(x, phi.salt, temp, ...) {
    # ddf / dphi
    return(yx.numdiff(
        free.energy.df,
        x,
        phi.salt = phi.salt,
        temp = temp,
        ...
    ))
}
free.energy.dddf <- function(x, phi.salt, temp, ...) {
    # dddf / dphi
    return(yx.numdiff(
        free.energy.ddf,
        x,
        phi.salt = phi.salt,
        temp = temp,
        ...
    ))
}
free.energy.funs <- function(phi, phi.salt, temp, ...) {
    # f df ddf dddf
    return(data.frame(
        phi = phi,
        f = sapply(phi, free.energy.f, phi.salt, temp, ...),
        df = sapply(phi, free.energy.df, phi.salt, temp, ...),
        ddf = sapply(phi, free.energy.ddf, phi.salt, temp, ...),
        dddf = sapply(phi, free.energy.dddf, phi.salt, temp, ...)
    ))
}


# Analytical

pkpq <- function(..., sysprop = NULL) {
    # kpq = phi_q / pho_p
    # arg <- ifelse(missing(arg), list(...), arg)
    if (is.null(sysprop))
        arg <- list(...)
    else
        arg <- sysprop
    
    return( (arg$size.ratio[2] * arg$sigma[2] * arg$polymer.num[2] * arg$molar.ratio[2]) /
        (arg$size.ratio[1] * arg$sigma[1] * arg$polymer.num[1] * arg$molar.ratio[1]) )
}
pS <- function(sigma.p, sigma.q, kpq) {
        return((sigma.p + sigma.q * kpq) / (1 + kpq))
    }
pA <- function(r.p, r.q, kpq) {
        return((1 / r.p + kpq / r.q) * (1 / (1 + kpq)))
    }
pB <- function(r.p, r.q, kpq) {
    return(-log(1 + kpq) / (r.p * (1 + kpq)) - kpq * (log(1 + kpq) - log(kpq)) /
               ((r.q) * (1 + kpq)))
}
pp <- function(..., sysprop = NULL) {
    if (is.null(sysprop)) {
        arg <- list(...)
        
    } else {
        arg <- sysprop
    }
    out <- list()
    kpq <- pkpq(sysprop = arg)
    out$kpq <- kpq
    out$S <- pS(arg$sigma[1], arg$sigma[2], kpq)
    out$A <- pA(arg$polymer.num[1], arg$polymer.num[2], kpq)
    out$B <- pB(arg$polymer.num[1], arg$polymer.num[2], kpq)
    return(out)
}

gibbs <- function(phi.polymer, phi.salt, ...) {
    # ... = alpha, sigma, polymer.num, kpq
    # return F/NkT
    arg <- list(...)
    alpha <- arg$alpha
    para <- pp(sysprop = arg)
    kpq <- para$kpq
    phi <- c(phi.polymer/(1+kpq), phi.polymer*kpq/(1+kpq), phi.salt*0.5, phi.salt*0.5, 1-phi.salt-phi.polymer)
    return(
        (k.vol * kkB * arg$temp) / k.water.size ^ 3 * (
            -alpha * (para$S * phi.polymer + phi.salt) ^ 1.5 +
                para$A * phi.polymer * log(phi.polymer) +
                para$B * phi.polymer +
                phi.salt * log(0.5 * phi.salt) +
                (1 - phi.polymer - phi.salt) * log(1 - phi.polymer - phi.salt) +
                t(phi) %*% arg$Chi %*% phi
        )
    )
}
gibbs.d <- function(phi.polymer, phi.salt, ...) {
    # ... = alpha, sigma, polymer.num, kpq
    arg <- list(...)
    alfa <- arg$alpha
    Chi <- arg$Chi
    para <- pp(sysprop = arg)
    kpq <- para$kpq
    phi <- c(phi.polymer/(1+kpq), phi.polymer*kpq/(1+kpq), phi.salt*0.5, phi.salt*0.5, 1-phi.salt-phi.polymer)
    return(
        (k.vol * kkB * arg$temp) / k.water.size ^ 3 * (
            -alfa * 1.5 * (para$S * phi.polymer + phi.salt) ^ 0.5 * para$S +
                para$A * log(phi.polymer) + para$A +
                para$B -
                log(1 - phi.polymer - phi.salt) - 1 +
                2*Chi[1,1]*kpq/(1+kpq)^2 * phi.polymer + 2/(1+kpq)*Chi[1,1:5]%*%phi - 2/(1+kpq)*Chi[1,1]*phi[1]
        )
    )
}
gibbs.dd <- function(phi.polymer, phi.salt, ...) {
    # ... = alpha, sigma, polymer.num, kpq
    arg <- list(...)
    alpha <- arg$alpha
    para <- pp(sysprop = arg)
    kpq <- para$kpq
    return(
        (k.vol * kkB * arg$temp) / k.water.size ^ 3 * (
            -alpha * 0.75 * (para$S * phi.polymer + phi.salt) ^ (-0.5) * para$S ^ 2 +
                para$A * 1 / phi.polymer +
                1 / (1 - phi.polymer - phi.salt) +
                2*arg$Chi[1,1]*kpq/(1+kpq)^2
        )
    )
}
gibbs.ddd <- function(phi.polymer, phi.salt, ...) {
    # ... = alpha, sigma, polymer.num, kpq
    arg <- list(...)
    alpha <- arg$alpha
    para <- pp(sysprop = arg)
    return(
        (k.vol * kkB * arg$temp) / k.water.size ^ 3 * (
            alpha * 0.375 * (para$S * phi.polymer + phi.salt) ^ (-1.5) * para$S ^ 3 -
                para$A * 1 / phi.polymer ^ 2 +
                1 / (1 - phi.polymer - phi.salt) ^ 2
        )
    )
}
gibbs.pfps <- function(phi.polymer, phi.salt, ...) {
    # -para.alpha * 1.5 * (para.S * phip + phis) ** 0.5 +
    # log(0.5 * phis) + 1 +
    # -log(1 - phip - phis) - 1
    arg <- list(...)
    Chi <- arg$Chi
    para <- pp(sysprop = arg)
    kpq <- para$kpq
    return((k.vol * kkB * arg$temp) / k.water.size ^ 3 * (
        -arg$alpha * 1.5 * (para$S * phi.polymer + phi.salt) ^ 0.5 +
            log(0.5 * phi.salt) + 1 +
            -log(1 - phi.polymer - phi.salt) - 1 +
            (Chi[1,3]+Chi[1,4]-2*Chi[1,5])/(1+kpq) * phi.polymer
    ))
}
gibbs.pdfps <- function(phi.polymer, phi.salt, ...) {
    #         -para.alpha * 0.75 * (para.S * phip + phis) ** (-0.5) * para.S +
    # 1/(1 - phip - phis)
    arg <- list(...)
    Chi <- arg$Chi
    para <- pp(sysprop = arg)
    kpq <- para$kpq
    return((k.vol * kkB * arg$temp) / k.water.size ^ 3 * (
        -arg$alpha * 0.75 * (para$S * phi.polymer + phi.salt) ^ (-0.5) * para$S +
            1 / (1 - phi.polymer - phi.salt) +
            (Chi[1,3]+Chi[1,4]-2*Chi[1,5])/(1+kpq)
    ))
}
gibbs.funs <- function(phi.polymer, phi.salt, ...) {
    out <-
        expand.grid(phi.polymer = phi.polymer, phi.salt = phi.salt) %>%
        mutate(
            f = gibbs(phi.polymer, phi.salt, ...),
            df = gibbs.d(phi.polymer, phi.salt, ...),
            ddf = gibbs.dd(phi.polymer, phi.salt, ...),
            dddf = gibbs.ddd(phi.polymer, phi.salt, ...)
        )
    return(out)
}
binodal.curve.fun <- function(phi.polymer, ...) {
    c(
        gibbs.d(phi.polymer[1], ...) - gibbs.d(phi.polymer[2], ...),
        (phi.polymer[2] - phi.polymer[1]) * gibbs.d(phi.polymer[1], ...) - (gibbs(phi.polymer[2], ...) - gibbs(phi.polymer[1], ...))
    )
}


## From concentration to Phi

phi <- function(conc, size.ratio, length.water, polym.num) {
    # units: mol/m^3, 1, m, 1
    phi <- kNa * conc * polym.num * size.ratio * length.water ^ 3
    return(phi)
}
get.phi <- function(conc, system.properties) {
    phis <-
        phi(
            conc[1:4],
            system.properties$size.ratio[1:4],
            system.properties$water.size,
            system.properties$polym.num[1:4]
        )
    return(c(phis, 1 - sum(phis)))
}
get.phis <- function(concs, system.properties) {
    return((lapply(concs, get.phi, system.properties)))
}


## From Concentration to Phi to Free energy

get.free.energy <- function(temp, conc, Chi, Para) {
    # units: K, mol/m^3, 1
    phi <- get.phi(conc, Para)
    Fen <- Fen(phi, Para$polym.num, Para$size.ratio)
    Fel <- Fel(alpha = alpha(temp, Para$size), Para$charge.den, phi)
    Fchi <- Fchi(Chi, phi)
    G <-
        k.water.size ^ 3 / (k.vol * kkB * temp) * (Fen + Fel + Fchi)
    return(G)
}
get.conc <- function(tot.conc, salt.conc, Para) {
    # units: mg/mL, mM
    protein.conc <-
        tot.conc * 1E3 * Para$mass.ratio[1] / Para$MW[1]  # mol/m^3
    rna.conc <-
        tot.conc * 1E3 * Para$mass.ratio[2] / Para$MW[2]  # mol/m^3
    salt.conc <- salt.conc * 1  # mol/m^3
    conc <-
        c(protein.conc, rna.conc, salt.conc, salt.conc, k.water.conc)
    return(conc)
}
get.free.energy_ <- function(temps, tot.concs, salt.concs, Chis, Paras) {
        # input: list of subjects
        # units: K, mg/mL, mL
        temps <- temps  # K
        ds <-
            expand.grid(
                temp = temps,
                tConc = tot.concs,
                sConc = salt.concs,
                Chi = Chis,
                Para = Paras
            ) %>%
            rowwise() %>%
            mutate(free.energy = get.free.energy(temp, get.conc(tConc, sConc, Para), Chi, Para)) %>%
            mutate(entropy = Fen(
                phis = get.phi(conc = get.conc(tConc, sConc, Para), Para),
                Para$polym.num,
                Para$size.ratio
            )) %>%
            mutate(enthalpy.el = Fel(
                alpha = alpha(temp, size = Para$size),
                sigma = Para$charge.den,
                get.phi(conc = get.conc(tConc, sConc, Para), Para)
            ))
        return(ds)
    }


## math

normed <- function(x) {
    return((x - x[1]))
}
yx.numdiff <- function(f, ...) {
    return(grad(f, ..., method = 'central'))
}
yx.fsolve <- function (f,
                  x,
                  J = NULL,
                  maxiter = 100,
                  tol = .Machine$double.eps ^ (0.5),
                  ...) {
            x0 <- x
            if (!is.numeric(x0))
                stop("Argument 'x0' must be a numeric vector.")
            x0 <- c(x0)
            fun <- match.fun(f)
            f <- function(x)
                fun(x, ...)
            n <- length(x0)
            m <- length(f(x0))
            if (!is.null(J)) {
                Jun <- match.fun(J)
                J <- function(x)
                    Jun(x, ...)
            }
            else {
                J <- function(x)
                    jacobian(f, x)
            }
            if (m == n) {
                sol = broyden(f,
                              x0,
                              J0 = J(x0),
                              maxiter = maxiter,
                              tol = tol)
                xs <- sol$zero
                fs <- f(xs)
            }
            else {
                sol <- gaussNewton(x0,
                                   f,
                                   Jfun = J,
                                   maxiter = maxiter,
                                   tol = tol)
                xs <- sol$xs
                fs <- sol$fs
                if (fs > tol)
                    warning("Minimum appears not to be a zero -- change starting point.")
            }
            return(list(x = xs, fval = fs))
        }
yx.nr <- function(f,
                  x,
                  J,
                  ...,
                  epsilon = 1E-10,
                  maxiter = 1E3) {
    # Newton Raphson method
    iter <- 1
    new.guess <- x  #TODO
    while (iter < maxiter) {
        iter <- iter + 1
        guess <- new.guess
        jacobian <-  J(guess, ...)
        invjac <- inv(((jacobian)))
        # new.guess <- guess + c(1, -1) * (guess[2]-guess[1])/maxiter
        f1 <- f(guess, ...)
        if (DEBUG)
            print(c(' f1', f1))
        new.guess <- guess - 0.1 * invjac %*% f1
        if (DEBUG)
            print(c('invjac %*% f1', invjac %*% f1))
        if (abs(new.guess[1] - guess[1]) < epsilon &&
            abs(new.guess[2] - guess[2]) < epsilon) {
            break
        } else {
            if (DEBUG)
                print('guess missed')
        }
    }
    if (iter >= maxiter) {
        print('maxiter reached')
        return(list(x = c(NA, NA)))
    } else {
        print('nr completed')
        return(list(x = new.guess))
    }
    return(list(x = ifelse(iter == maxiter, new.guess, c(NA, NA))))
}
stupid.fsolve <- function(f, x, x.critic, epsilon = 1E-10, ...) {
    # x[1] --- x.critic
    for (x1 in x[which(x < x.critic)]) {
        # x.critic --- x[2]
        for (x2 in x[which(x.critic <= x)]) {
            f.out <- f(x = c(x1, x2), ...)
            if (abs(f.out[1]) < epsilon && abs(f.out[2]) < epsilon) {
                return(list(
                    x = c(x1, x2),
                    fval = c(epsilon, epsilon)
                ))
            }
        }
    }
    return(list(x = c(NA, NA), fval = c(1E15, 1E15)))
}


## Critical point
c.point.temp <- function(sysprop, fitting.para) {
    temp <- seq(273.15, 333.15, 0.1)
    do.call(rbind, lapply(temp, function(temp) {
        alpha <- get.alpha(temp, sysprop$water.size)
        polymer.num <- sysprop$polymer.num
        sigma <- sysprop$sigma
        size.ratio <- sysprop$size.ratio
        Chi <- sysprop$Chi
        molar.ratio <- sysprop$molar.ratio
        ds <- critical.point_(guess.critical.point = fitting.para$critical.point.guess, temp = temp, polymer.num = polymer.num,
                              alpha = alpha, sigma = sigma, size.ratio = size.ratio, Chi = Chi, molar.ratio = molar.ratio)
        ds$temp <-  temp
        return(as.data.frame(ds))
    }))
}
c.point.temp.fun <- function(c.point.temp.ds) {
    return(function(temp) {
        list(
             phi.polymer = spline(c.point.temp.ds$temp, c.point.temp.ds$phi.polymer, xout = temp)$y,
            phi.salt = spline(c.point.temp.ds$temp, c.point.temp.ds$phi.salt, xout = temp)$y)
    })
}
critical.point.fun_ <- function(phis, ...) {
    arg <- list(...)
    return(c(
        gibbs.dd(phi.polymer = phis[1], phi.salt = phis[2], ...),
        gibbs.ddd(phi.polymer = phis[1], phi.salt = phis[2], ...)
    ))
}
critical.point_ <- function(...) {
    # output list of critical points
    arg <- list(...)
    guess <- arg$guess.critical.point
    out <- fsolve(f = critical.point.fun_,
                  x0 = guess,
                  ...)
    return(list(phi.polymer = out$x[1], phi.salt = out$x[2]))
}


## Binodal curve

binodal.curve.jacobian.temp.phi <-
    function(x, phi.polymer.1, phi.salt, sysprop) {
        phi1 <- phi.polymer.1
        phis <- phi.salt
        t <- x[1]
        phi2 <- x[2]
        
        sp <- sysprop
        sab <- pp(sysprop = sp)
        
        a <-
            get.alpha(t, sp$water.size)
        dadT <- alpha.dT(t, sp$water.size)
        S <- sab$S
        A <- sab$A
        B <- sab$B
        
        matrix((k.vol * kkB * t) / k.water.size ^ 3 * c(
            (-dadT * 1.5 * sqrt(S * phi1 + phis) * S) -
                (-dadT * 1.5 * sqrt(S * phi2 + phis) * S),
            
            (phi2 - phi1) * (-dadT * 1.5 * sqrt(S * phi2 + phis) * S) -
                ((-dadT * (
                    sqrt(S * phi2 + phis) ^ 3
                )) - (-dadT * (
                    sqrt(S * phi1 + phis) ^ 3
                ))),-(
                    -a * 0.75 * 1 / sqrt(S * phi2 + phis) * S ^ 2 +  A * 1 / phi2 +  1 / (1 -
                                                                                              phi2 - phis)
                ),
            
            (
                -a * 1.5 * sqrt(S * phi1 + phis) * S + A * log(phi1) +-log(1 - phi1 - phis)
            ) -
                (
                    -a * 1.5 * sqrt(S * phi2 + phis) * S + A * log(phi2) +-log(1 - phi2 - phis)
                )
        ), 2, 2)
    }
binodal.curve.fun.temp.phi <-
    function(x, phi.polymer.1, phi.salt, sysprop) {
        phi1 <- phi.polymer.1
        phis <- phi.salt
        t <- x[1]
        phi2 <- x[2]
        
        sp <- sysprop
        sab <- pp(sysprop = sp)
        
        a <-
            get.alpha(t, sp$water.size)
        dadT <- alpha.dT(t, sp$water.size)
        S <- sab$S
        A <- sab$A
        B <- sab$B
        
        return(c(
            (
                -a * 1.5 * sqrt(S * phi1 + phis) * S + A * log(phi1) +-log(1 - phi1 - phis)
            ) -
                (
                    -a * 1.5 * sqrt(S * phi2 + phis) * S + A * log(phi2) +-log(1 - phi2 - phis)
                ),
            (phi2 - phi1) * (
                -a * 1.5 * sqrt(S * phi2 + phis) * S + A * log(phi2) +-log(1 - phi2 - phis)
            ) -
                ((
                    -a * sqrt(S * phi2 + phis) ^ 3 + A * phi2 * log(phi2) + B * phi2 + (1 -
                                                                                            phi2 - phis) * log(1 - phi2 - phis)
                ) -
                    (
                        -a * sqrt(S * phi1 + phis) ^ 3 + A * phi1 * log(phi1) + B * phi1 + (1 -
                                                                                                phi1 - phis) * log(1 - phi1 - phis)
                    )
                )
        ))
    }
binodal.curve.jacobian_ <- function(x, phi.polymer.2, ...) {
    phi.polymer.1 <- x[1]
    phi.salt <- x[2]
    return(matrix(
        c(
            gibbs.dd(phi.polymer.1, phi.salt, ...),
            
            gibbs.d(phi.polymer.1, phi.salt, ...) - gibbs.d(phi.polymer.2, phi.salt, ...),
            
            gibbs.pdfps(phi.polymer.1, phi.salt, ...) - gibbs.pdfps(phi.polymer.2, phi.salt, ...),
            
            (phi.polymer.2 - phi.polymer.1) * gibbs.pdfps(phi.polymer.2, phi.salt, ...) -
                (
                    gibbs.pfps(phi.polymer.2, phi.salt, ...) - gibbs.pfps(phi.polymer.1, phi.salt, ...)
                )
        ),
        nrow = 2,
        ncol = 2
    ))
}
binodal.curve.fun_ <- function(x, phi.polymer.2, ...) {
    phi.polymer.1 <- x[1]
    phi.salt <- x[2]
    return(binodal.curve.fun(
        c(phi.polymer.1, phi.polymer.2),
        phi.salt = phi.salt,
        ...
    ))
}
binodal.curve_ <- function(..., sysprop = NULL, fitting.para = NULL) {
    # generate binodal curve
    if (is.null(sysprop)) arg <- list(...)
    else arg <- sysprop
    
    c.point <- critical.point_(...)
    if (DEBUG) print(c('critical point', c.point))
    
    binodal.guess <- arg$binodal.guess
    binodal.guess[2] <- c.point$phi.salt * 0.01
    
    # search binodal point return c(phi.polymer, phi.salt)
    # phi.polymer.2.seq <- c(seq( arg$sampling.gap, arg$sampling.gap*1e3, arg$sampling.gap),
    #                        seq( arg$sampling.gap * 1e3, c.point$phi.polymer, arg$sampling.gap * 1e4))
    n <- floor(0.5 * (1 + sqrt(1 + 8 * c.point$phi.polymer / arg$sampling.gap)))
    phi.polymer.2.seq <- arg$sampling.gap * seq(1, n, 1) * seq(2, n + 1, 1) * 0.5
    
    output <- list()
    for (phi.polymer.2 in phi.polymer.2.seq) {
        roots <- nleqslv::nleqslv (
            binodal.guess,
            binodal.curve.fun_,
            binodal.curve.jacobian_,
            phi.polymer.2 = phi.polymer.2,
            ...
        )
        output[[length(output) + 1]] <-
            c(
                phi.polymer.1 = roots$x[1],
                phi.polymer.2 = phi.polymer.2,
                phi.salt = roots$x[2],
                f1 = roots$fval[1],
                f2 = roots$fval[2]
            )
        binodal.guess <- roots$x
    }
    # output <- lapply((phi.polymer.2.seq), function(phi.polymer.2) {
    #     roots <- nleqslv::nleqslv (
    #         binodal.guess,
    #         binodal.curve.fun_,
    #         binodal.curve.jacobian_,
    #         phi.polymer.2 = phi.polymer.2,
    #         ...
    #     )
    #     # binodal.curve.guess <- roots$x
    #     c(
    #         phi.polymer.1 = roots$x[1],
    #         phi.polymer.2 = phi.polymer.2,
    #         phi.salt = roots$x[2],
    #         f1 = roots$fval[1],
    #         f2 = roots$fval[2]
    #     )
    # })
    
    p2 <- as.data.frame.matrix(do.call(rbind, output)) %>%
        # requiring the dense phase should be larger than dilute phase
        filter(phi.polymer.1 > max(phi.polymer.2))
    
    ds <- data.frame(
        phi.polymer = c(p2$phi.polymer.2, rev(p2$phi.polymer.1)),
        phi.salt = c(p2$phi.salt, rev(p2$phi.salt)),
        phase = factor(c(rep('dilute', nrow(p2)), rep('dense', nrow(p2))),
        levels = c('dilute', 'dense')),
        critic.polymer = c.point$phi.polymer,
        critic.salt = c.point$phi.salt
    )
    return(ds)
}


## Spinodal curve

spinodal.curve <-
    function(phi.polymer,
             phi.salt,
             x.axis = 'phi.polymer',
             temp,
             alpha,
             sigma,
             Chi = 0,
             polymer.num,
             size.ratio,
             curve.type = 'spinodal') {
        if (x.axis == 'phi.polymer') {
            paras <- phi.salt
            root.interval <- phi.polymer
        } else {
            paras <- phi.polymer
            root.interval <- phi.salt
        }
        
        # get roots at all the x
        root.and.x <- lapply(paras, function(para) {
            # g = free energy
            g <- sapply(root.interval, function(y) {
                free.energy(
                    phi.polymer = ifelse(x.axis == 'phi.polymer', y, para),
                    phi.salt = ifelse(x.axis == 'phi.polymer', para, y),
                    temp = temp,
                    alpha = alpha,
                    sigma = sigma,
                    Chi = Chi,
                    polymer.num = polymer.num,
                    size.ratio = size.ratio
                )
            })
            dg <- diff(g) / diff(root.interval)
            ddg <- diff(dg) / diff(root.interval)[-1]
            dddg <- diff(ddg) / diff(root.interval)[-2:-1]
            # function to find root of ddg == 0
            if (curve.type == 'spinodal') {
                func <- ddg
                func.r <- root.interval[-2:-1]
            } else if (curve.type == 'critic') {
                func <- dddg
                func.r <- root.interval[-3:-1]
            } else if (curve.type == 'binodal') {
                func <- dg
                func.r <- root.interval[-1]
            } else if (curve.type == 'deltaG') {
                func <- g
                func.r <- root.interval
            } else {
                return()
            }
            sp.fun <- function(x) {
                spline(x = func.r,
                       y = func,
                       xout = x)$y
            }
            # all the roots using N-R method
            root <- uniroot.all(f = sp.fun, interval = range(func.r))
            # output as a list of two vector
            ifelse(
                length(root) == 0,
                output <- data.frame(root = NA, para = NA),
                output <-
                    data.frame(root = root, para = rep(para, length(root)))
            )
            return(output)
        })
        return(do.call(rbind, root.and.x))
    }

## Deprecated
{
    # critical.point.fun <- function(x, alpha, sigma, Chi, temp,
    #                                polymer.num, size.ratio ) {
    #   # function for critical.point()
    #   phi.polymer <- x[1]
    #   phi.salt <- x[2]
    #
    #   output <- c(free.energy.ddf(x = x[1], phi.salt = x[2],
    #                               alpha = alpha, sigma = sigma, Chi = Chi,
    #                               temp = temp, polymer.num = polymer.num, size.ratio =size.ratio),
    #               free.energy.dddf(x = x[1], phi.salt = x[2],
    #                               alpha = alpha, sigma = sigma, Chi = Chi,
    #                               temp = temp, polymer.num = polymer.num, size.ratio =size.ratio))
    #   return(sum(output^2))
    # }
    # critical.point <- function(alpha, sigma, Chi, temp,
    #                            polymer.num, size.ratio) {
    #     # generate critical points //TODO:
    #
    #     # guess: GOOD LUCK!!!
    #   guess.phi.salt <- seq(0.02, 0.15, 1e-3)
    #   guess.phi.polymer <- seq(0.0001, 0.1, 1e-3)
    #   sp.curve <- spinodal.curve(guess.phi.polymer, guess.phi.salt, 'phi.polymer',
    #                              temp, alpha, sigma, 0, polymer.num, size.ratio)
    #     # call non-linear-equation-set solver, return c(phi.polymer, phi.salt)
    #   fsolve(f = critical.point.fun, x0 = guess,
    #          alpha = alpha, sigma = sigma, Chi = Chi,
    #          temp=temp, polymer.num = polymer.num, size.ratio = size.ratio)
    # }
    # binodal.curve.jacobian <- function(x, ...){
    #   # return(jacobian(binodal.curve.fun, x, ...))
    #   return(matrix(c(
    #     gibbs.dd(x[1], ...),  # dF1/dX1
    #     -gibbs.d(x[1], ...) + (x[2]-x[1])*gibbs.dd(x[1], ...) + gibbs.d(x[1], ...),  # dF2/dX1
    #     -gibbs.dd(x[2], ...),  # dF1/dX2
    #     gibbs.d(x[1], ...) - gibbs.d(x[2], ...)  # dF2/dX2
    #   ), nrow = 2, ncol = 2))
    # }
    # binodal.curve <- function(phi.polymer.seq, ...) {
    #   # generate binodal curve
    #   arg <- list(...)
    #   # binodal.curve.guess <- range(phi.polymer.seq)
    #   binodal.curve.guess <- c(1E-5, 0.5)
    #   # critical point
    #   c.point <- critical.point_(...)
    #   if (DEBUG) print(c('critical point', c.point))
    #   # range of phi.salt = seq(0, critical.salt)
    #   phi.salt.seq <- seq(1e-15, c.point$phi.salt, 1e-2)
    #   if (DEBUG) phi.salt.seq <- c(0.13)
    #   # search binodal point return c(phi.polymer, phi.salt)
    #   output <- list()
    #   for(i in seq_along(phi.salt.seq)) {
    #   # }
    #   # output <- lapply(seq_along(phi.salt.seq), function(i) {
    #     phi.salt <- phi.salt.seq[i]
    #     # roots <- stupid.fsolve(f = binodal.curve.fun, x = phi.polymer.seq, x.critic = c.point$phi.polymer,
    #     #                        phi.salt = phi.salt, ...)
    #     roots <-
    #         yx.nr (
    #          binodal.curve.fun,
    #          binodal.curve.guess,
    #         binodal.curve.jacobian,
    #         phi.salt = phi.salt,
    #         ...
    #       )
    #     # binodal.curve.guess <- roots$x
    #     print(roots)
    #     # return(
    #     output[[i]] <- c(
    #         phi.salt = phi.salt,
    #         phi.polymer1 = roots$x[1],
    #         phi.polymer2 = roots$x[2],
    #         f1 = roots$fval[1],
    #         f2 = roots$fval[2]
    #       )
    #     # )
    #   }
    #   # )
    #   return(do.call(rbind, output))
    # }
}