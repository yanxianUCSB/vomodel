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

free.energy <-
  function(phi.polymer,
           phi.salt,
           temp,
           alpha,
           sigma,
           Chi,
           polymer.num,
           size.ratio) {
    phi <-
      c(
        phi.polymer * 0.5,
        phi.polymer * 0.5,
        phi.salt * 0.5,
        phi.salt * 0.5,
        1 - phi.polymer - phi.salt
      )
    g <-
      k.water.size^3/(k.vol*kkB*temp) * 
        (Fen(phi, polymer.num, size.ratio) + Fel(alpha, sigma, phi) + Fchi(Chi, phi))
    return(g)
  }

## Concentration to Phi to Free energy
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
  protein.conc <- tot.conc * 1E3 * Para$mass.ratio[1] / Para$MW[1]  # mol/m^3
  rna.conc <- tot.conc * 1E3 * Para$mass.ratio[2] / Para$MW[2]  # mol/m^3
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

