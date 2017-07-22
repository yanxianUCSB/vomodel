# Voorn-Overbeek Modeling
kNa <<- 6.02E23
kkB <<- 1.38064852E-23
ke  <<- 1.60217662E-19
kEr <<- 4 * pi * 80 * 8.85E-12 # F/m; vaccuum permitivity = 8.85E-12
k.water.size <<- 0.31E-9  # 0.31nm as a size of a water molecule
k.vol <<- 120E-9
k.water.conc  <<- 1000 / 18.01528 * 1000  # water concentration

# # c(+, -, sp, sn, w)
# Para <- list()
# Para$size <- k.water.size
# Para$charge.den <- c(11/207, 1, 1, 1, 0)
# Para$polym.num <- c(207, 900E3 / (324 - 18), 1, 1, 1)  # polym.num. REF[Yanxian's Notebook]
# Para$size.ratio <- c(1, 1, 1, 1, 1)

## From Phi to free.energy and phase diagram 

Fen <- function(phis, rs, ws) {
  # units: 1, 1, 1
  Fen <- sum(phis * log(phis) / (ws * rs))
}

alpha <- function(temp, size) {
  # units: K, m
  alpha <-
    2 / 3 * sqrt(pi) * ((ke ^ 2 / (kEr * kkB * temp)) / size) ^ (3 / 2)
}

Fel <- function(alpha, sigma, phi) {
  # units: 1, 1, 1
  Fel <- 0 - alpha * sum(sigma * phi) ^ (1.5)
}

Fchi <- function(Chi, phi) {
  # units: 1, 1
  Fchi <- sum(Chi * (phi %*% t(phi)))
}

free.energy <- function(phi.polymer, phi.salt, alpha, sigma, Chi, temp,
                        polymer.num, size.ratio) {
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
      k.water.size ^ 3 / (k.vol * kkB * temp) *
      (Fen(phi, polymer.num, size.ratio) + Fel(alpha, sigma, phi) + Fchi(Chi, phi))
    return(g)
}

free.energy.f <- function(x, phi.salt, alpha, sigma, Chi, temp,
                          polymer.num, size.ratio) {
  # free energy f
  free.energy(phi.polymer = x, phi.salt, alpha, sigma, Chi, temp,
              polymer.num, size.ratio)
}

free.energy.df <- function(x, phi.salt, alpha, sigma, Chi, temp,
                           polymer.num, size.ratio) {
  # df / dphi
  return( grad(free.energy.f, x, 
               phi.salt = phi.salt, 
               alpha = alpha, sigma = sigma, Chi = Chi, 
               temp = temp, polymer.num = polymer.num, size.ratio =size.ratio))
}

free.energy.ddf <- function(x, phi.salt, alpha, sigma, Chi, temp,
                            polymer.num, size.ratio) {
  return(grad(free.energy.df, x, 
              phi.salt = phi.salt, 
              alpha = alpha, sigma = sigma, Chi = Chi, 
              temp = temp, polymer.num = polymer.num, size.ratio =size.ratio))
}

free.energy.dddf <- function(x, phi.salt, alpha, sigma, Chi, temp,
                             polymer.num, size.ratio) {
  return(grad(free.energy.ddf, x, 
              phi.salt = phi.salt, 
              alpha = alpha, sigma = sigma, Chi = Chi, 
              temp = temp, polymer.num = polymer.num, size.ratio =size.ratio))
}

free.energy.funs <- function(phi, phi.salt, alpha, sigma, Chi, temp,
                             polymer.num, size.ratio) {
  # f df ddf dddf
  return(data.frame(phi = phi,
                    f = sapply(phi, free.energy.f, phi.salt, alpha, sigma, Chi, temp,
                               polymer.num, size.ratio),
                    df = sapply(phi, free.energy.df, phi.salt, alpha, sigma, Chi, temp,
                                polymer.num, size.ratio),
                    ddf = sapply(phi, free.energy.ddf, phi.salt, alpha, sigma, Chi, temp,
                                 polymer.num, size.ratio),
                    dddf = sapply(phi, free.energy.dddf, phi.salt, alpha, sigma, Chi, temp,
                                  polymer.num, size.ratio)))
}

## From concentration to Phi
phi <- function(conc, size.ratio, length.water, polym.num) {
  # units: mol/m^3, 1, m, 1
  phi <- kNa * conc * polym.num * size.ratio * length.water ^ 3
  return(phi)
}

get.phi <- function(conc, Para) {
  phis <-
    phi(conc[1:4], Para$size.ratio[1:4], Para$size, Para$polym.num[1:4])
  return(c(phis, 1 - sum(phis)))
}

get.phis <- function(concs, Para) {
  return((lapply(concs, get.phi, Para)))
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

get.free.energy_ <-
  function(temps, tot.concs, salt.concs, Chis, Paras) {
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

normed <- function(x) {
  return((x - x[1]))
}


critical.point.fun <- function(x, alpha, sigma, Chi, temp,
                               polymer.num, size.ratio ) {
  # function for critical.point()
  phi.polymer <- x[1]
  phi.salt <- x[2]
  
  output <- c(free.energy.ddf(x = x[1], phi.salt = x[2], 
                              alpha = alpha, sigma = sigma, Chi = Chi, 
                              temp = temp, polymer.num = polymer.num, size.ratio =size.ratio),
              free.energy.dddf(x = x[1], phi.salt = x[2], 
                              alpha = alpha, sigma = sigma, Chi = Chi, 
                              temp = temp, polymer.num = polymer.num, size.ratio =size.ratio))
  return(sum(output^2))
}

critical.point <- function(alpha, sigma, Chi, temp,
                           polymer.num, size.ratio) {
    # generate critical points //TODO:
    
    # guess: GOOD LUCK!!!
  guess <- c(0.02, 0.12)
    # call non-linear-equation-set solver, return c(phi.polymer, phi.salt)
  fsolve(f = critical.point.fun, x0 = guess, 
         alpha = alpha, sigma = sigma, Chi = Chi,
         temp=temp, polymer.num = polymer.num, size.ratio = size.ratio)
  }

# binodal.curve.fun <- function(s, y, interval) {
#   # function for binodal.curve()
#   # s.t. f - ax -b = 0, f' - a = 0
#   a <- s[1]
#   b <- s[2]
#   x1 <- s[3]
#   x2 <- s[4]
#   g <- function(x) spline(interval, y + a * interval + b, xout = x)$y
#   dg <- function(x) sapply(x, function(x) grad(g, x))
#   # return(sum(abs(dg(c(x1, x2)))) + sum(abs(g(c(x1, x2)))))
#   return(c(g(x1), g(x2), dg(x1), dg(x2)))
# }
binodal.curve.fun <- function(s, y, interval) {
  x1 <- s[1]
  x2 <- s[2]
  g <- function(x) spline(interval, y, xout = x)$y
  dg <- function(x) sapply(x, function(x) grad(g, x))
  a <- (g(x1) - g(x2)) / (x1 - x2)
  f <- (a - dg(x1)) ^ 2 + (a - dg(x2)) ^ 2
  return(f)
}
binodal.curve <- function(phi.polymer.seq, phi.salt.seq, ...) {
  # generate binodal curve
  
  # range of phi.salt = seq(0, critical.salt)
  
  # search binodal point return c(phi.polymer, phi.salt)
  y <- sapply(phi.polymer.seq, function(x) free.energy.f(x, phi.salt = 0.001,
                                             temp = temp,
                                             alpha = alpha,
                                             sigma = sigma,
                                             Chi = 0,
                                             polymer.num = polymer.num,
                                             size.ratio = size.ratio))
  output <- fsolve(x0 = c( 0.02, 0.15), 
        f = binodal.curve.fun, y = y, interval = phi.polymer.seq, tol = 1e-18)
  return(output)
}
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
        spline(x = func.r, y = func, xout = x)$y
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
