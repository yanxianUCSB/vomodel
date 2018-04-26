#' This file contains scripts to fit experimental data with Flory Huggins Chi para.
#'       The definition of gibbs free energy is contained at vomodel.R
#'       
source('vomodel.R')
# Fn --------------------
#' fn = fn(x), here x=phi1, other parameters are known
get_chi   <- function(phi1, phi2, phi3, temp, N, lattice, sigma, p2p1) {
  if(phi1 == phi2) stop(msg = 'phi == phi2 for get_chi!')
  phi_1 <-
    c(phi1, p2p1 * phi1, phi3, phi3, 1 - phi1 - p2p1 * phi1 - 2 * phi3)
  phi_2 <-
    c(phi2, p2p1 * phi2, phi3, phi3, 1 - phi2 - p2p1 * phi2 - 2 * phi3)
  d_entropy_1 <- d_entropy(phi_1, N)
  d_entropy_2 <- d_entropy(phi_2, N)
  d_enthalpy_1 <- d_enthalpy(phi_1, temp, lattice, sigma)
  d_enthalpy_2 <- d_enthalpy(phi_2, temp, lattice, sigma)
  under_1 <- (phi_1[5]*(1 + p2p1) - phi_1[1]*(1+p2p1)^2)
  under_2 <- (phi_2[5]*(1 + p2p1) - phi_2[1]*(1+p2p1)^2)
  chi <-
    ((d_entropy_1 - d_entropy_2) + (d_enthalpy_1 - d_enthalpy_2)) / 
    (under_2 - under_1)
  return(chi)
}  # should be good
fn        <- function(phi1, phi2, phi3, temp, N, lattice, sigma, p2p1){
  chi <- get_chi(phi1, phi2, phi3, temp, N, lattice, sigma, p2p1)
  phi_1 <-
    c(phi1, p2p1 * phi1, phi3, phi3, 1 - phi1 - p2p1 * phi1 - 2 * phi3)
  phi_2 <-
    c(phi2, p2p1 * phi2, phi3, phi3, 1 - phi2 - p2p1 * phi2 - 2 * phi3)
  diff <- 
    d_gibbs(phi_2, temp, N, lattice, sigma, chi) * (phi2 - phi1) -
    gibbs(phi_2, temp, N, lattice, sigma, chi) +
    gibbs(phi_1, temp, N, lattice, sigma, chi)
}
fn.cr        <- function(phi1, phi2, phi3, temp, N, lattice, sigma, p2p1){
  sigma[1] <- sigma.cr(sigma[1], temp)
  sigma[2] <- sigma.cr(sigma[2], temp)
  phi_1    <- c(phi1, p2p1 * phi1, phi3, phi3, 1 - phi1 - p2p1 * phi1 - 2 * phi3)
  phi_2    <- c(phi2, p2p1 * phi2, phi3, phi3, 1 - phi2 - p2p1 * phi2 - 2 * phi3)
  chi      <- get_chi(phi1, phi2, phi3, temp, N, lattice, sigma, p2p1)
  diff     <- 
    d_gibbs(phi_2, temp, N, lattice, sigma, chi) * (phi2 - phi1) -
    gibbs(phi_2, temp, N, lattice, sigma, chi) +
    gibbs(phi_1, temp, N, lattice, sigma, chi)
}
# Data preparation ---------------
# this part of code transform the experimental observations into volume fraction
# and temperature in Kelvin
get_expt  <- function(path_experiment = commandArgs(trailingOnly = T)[1]) {
  library(dplyr)
  # path_experiment <- '~/Desktop/dataset.csv'
  ds <- read.csv(path_experiment) %>% 
    mutate(
      phi1 = protein * 1E-6 * 207 / k.water.conc * 1000,  # k.water.conc in mol/m^3
      phi3 = 0.5 * (nacl + 20) * 1E-3 / k.water.conc * 1000,  # 20 mM monovalent buffer salt
      temp = cloudpoint + 273.15
    ) %>% 
    select(phi1, phi3, temp) 
  ds %>% 
    saveRDS('expt.rds')
  return(ds)
}
# Simulate ------------------
sim_      <- function(x){
  x <- as.numeric(x)
  phi2 <- x[1]
  phi3 <- x[2]
  temp <- x[3]
  initial_guesses <- INITIAL_GUESS
  phi1 <- NULL
  tryCatch({
    root <- uniroot(
      f = fn,
      interval = initial_guesses,
      phi2 = phi2,
      phi3 = phi3,
      temp = temp,
      N       = POLYNUM,
      lattice = k.water.size,
      sigma   = SIGMA,
      p2p1    = p2p1
    )
    phi1 <- root$root
  }, error=function(e){
    phi1 <- NA
  }, warnings=function(w){
  }, finally = {
    phi1 <- ifelse(is.null(phi1), NA, phi1)
  })
  return(phi1)
}
simulate  <- function(output.filename.no.ext = 'simulated'){
  # get_chi(...) for dplyr::mutate
  get_chi_ <- function(phi1, phi2, phi3, temp){
    seq(1, length(phi1), 1) %>% 
      sapply(function(i){
        get_chi(phi1 = phi1[i], phi2 = phi2[i], phi3 = phi3[i], temp = temp[i], 
                N = POLYNUM, lattice = k.water.size, sigma = SIGMA, p2p1 = p2p1)
      })
  }
  expt <- get_expt()
  phi_dense <- apply(as.matrix(expt), 1, sim_)
  if(is.null(phi_dense) | 
     is.na(unique(phi_dense))){
    print('No fit found!')
  } else {
    expt %>% 
      mutate(phi2 = phi_dense) %>% 
      mutate(chi = get_chi_(phi1, phi2, phi3, temp)) -> output
    output %>% 
      write.csv(paste0(output.filename.no.ext, '.csv'), row.names = F)
    output %>% 
      saveRDS(paste0(output.filename.no.ext, '.rds'))
  }
}
sim.cr_      <- function(x){
  x <- as.numeric(x)
  phi2 <- x[1]
  phi3 <- x[2]
  temp <- x[3]
  initial_guesses <- INITIAL_GUESS
  phi1 <- NULL
  tryCatch({
    root <- uniroot(
      f = fn.cr,
      interval = initial_guesses,
      phi2 = phi2,
      phi3 = phi3,
      temp = temp,
      N       = POLYNUM,
      lattice = k.water.size,
      sigma   = SIGMA,
      p2p1    = p2p1
    )
    phi1 <- root$root
  }, error=function(e){
    phi1 <- NA
  }, warnings=function(w){
  }, finally = {
    phi1 <- ifelse(is.null(phi1), NA, phi1)
  })
  return(phi1)
}
simulate.cr  <- function(output.filename.no.ext = 'simulated'){
  # get_chi(...) for dplyr::mutate
  get_chi_ <- function(phi1, phi2, phi3, temp){
    seq(1, length(phi1), 1) %>% 
      sapply(function(i){
        get_chi(phi1 = phi1[i], phi2 = phi2[i], phi3 = phi3[i], temp = temp[i], 
                N = POLYNUM, lattice = k.water.size, sigma = SIGMA, p2p1 = p2p1)
      })
  }
  expt <- get_expt()
  phi_dense <- apply(as.matrix(expt), 1, sim.cr_)
  if(is.null(phi_dense) | 
     is.na(unique(phi_dense))){
    print('No fit found!')
  } else {
    expt %>% 
      mutate(phi2 = phi_dense) %>% 
      mutate(chi = get_chi_(phi1, phi2, phi3, temp)) -> output
    output %>% 
      write.csv(paste0(output.filename.no.ext, '.csv'), row.names = F)
    output %>% 
      saveRDS(paste0(output.filename.no.ext, '.rds'))
  }
}
# Modeling ------------
simulate.fh <- function(){
  POLYNUM       <<- c(207, 2939,   1, 1, 1)
  SIGMA         <<- c(0, 0,   0, 0, 0)  # charge per monomer
  SIZE          <<- c(1, 1, 1, 1, 1) * 0.31E-9  # m; Voorn Overbeek 1967
  MW            <<- c(22e3, 900e3)
  MOLARRATIO    <<- c(2939,  11, 0.5, 0.5, 0)  # the polycation:polyanion and cation:anion molar ratio
  p2p1          <<- (POLYNUM[2] * MOLARRATIO[2] * SIZE[2]^3 ) /
    (POLYNUM[1] * MOLARRATIO[1] * SIZE[1]^3 )  # the volume fraction ratio between specie 1 and specie 2
  INITIAL_GUESS <<- c(0.1, 0.5)
  simulate('results/FH')
}
simulate.fhvo <- function(){
  POLYNUM       <<- c(207, 2939,   1, 1, 1)
  SIGMA         <<- c(11 / 207, 1,   1, 1, 0)  # charge per monomer
  SIZE          <<- c(1, 1, 1, 1, 1) * 0.31E-9  # m; Voorn Overbeek 1967
  MW            <<- c(22e3, 900e3)
  MOLARRATIO    <<- c(2939,  11, 0.5, 0.5, 0)  # the polycation:polyanion and cation:anion molar ratio
  p2p1          <<- (POLYNUM[2] * MOLARRATIO[2] * SIZE[2]^3 ) /
    (POLYNUM[1] * MOLARRATIO[1] * SIZE[1]^3 )  # the volume fraction ratio between specie 1 and specie 2
  INITIAL_GUESS <<- c(0.1, 0.5)
  simulate('results/FHVO')
}
simulate.fhvocr <- function() {
  POLYNUM       <<- c(207, 2939,   1, 1, 1)
  SIGMA         <<- c(11 / 207, 1,   1, 1, 0)  # charge per monomer
  SIZE          <<- c(1, 1, 1, 1, 1) * 0.31E-9  # m; Voorn Overbeek 1967
  MW            <<- c(22e3, 900e3)
  MOLARRATIO    <<- c(2939,  11, 0.5, 0.5, 0)  # the polycation:polyanion and cation:anion molar ratio
  p2p1          <<- (POLYNUM[2] * MOLARRATIO[2] * SIZE[2]^3 ) /
    (POLYNUM[1] * MOLARRATIO[1] * SIZE[1]^3 )  # the volume fraction ratio between specie 1 and specie 2
  INITIAL_GUESS <<- c(0.1, 0.5)
  simulate.cr('results/FHVOCR')
}
# Launch -----
simulate.fh()
simulate.fhvo()
simulate.fhvocr()
