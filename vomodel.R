# Constant --------------
kNa                         <<- 6.02E23     # Avogadro Number
kkB                         <<- 1.38E-23    # Boltzmann Constant
ke                          <<- 1.60E-19    # Elementary charge
kEr                         <<- 4 * pi * 80 * 8.85E-12 # F/m; vaccuum permitivity = 8.85E-12
k.water.size                <<- 0.31E-9     # calculated from 18 cm^3 / mol
k.na.size                   <<- 0.235E-9    # (Marcus 1988)
k.cl.size                   <<- 0.318E-9    # (Marcus 1988)
k.water.conc                <<- 1000 / 18.01528 * 1000  # water concentration
# # Parameters
# POLYNUM       <<- c(207, 2939,   1, 1, 1)
# SIGMA         <<- c(11 / 207, 1,   1, 1, 0)  # charge per monomer
# SIZE          <<- c(1, 1, 1, 1, 1) * 0.31E-9  # m; Voorn Overbeek 1967
# MW            <<- c(22e3, 900e3)
# MOLARRATIO    <<- c(2939,  11, 0.5, 0.5, 0)  # the polycation:polyanion and cation:anion molar ratio
# p2p1          <<- (POLYNUM[2] * MOLARRATIO[2] * SIZE[2]^3 ) / 
#   (POLYNUM[1] * MOLARRATIO[1] * SIZE[1]^3 )  # the volume fraction ratio between specie 1 and specie 2

# System definition -------------
#' System is defined as phi, temp, N, lattice, sigma, chi
get_alpha  <- function(temp, lattice) {
  2 / 3 * sqrt(pi) * ((ke ^ 2 / (kEr * kkB * temp)) / lattice) ^ (3 / 2)
}  
gibbs      <- function(phi, temp, N, lattice, sigma, chi) {
  entropy <- sum(phi / N * log(phi))
  enthalpy <- -1 * get_alpha(temp, lattice) * sum(sigma * phi) ^ 1.5
  enthalpy_chi <- chi * (phi[1] + phi[2]) * phi[5]
  gibbs <- entropy + enthalpy + enthalpy_chi
}
d_entropy  <- function(phi, N) {
  p2p1 <- phi[2] / phi[1]
  (1 / N[1] * log(phi[1]) + 1 / N[1]) + 
  (1 / N[2] * log(phi[2]) + 1 / N[2]) * p2p1 +
  (1 / N[5] * log(phi[5]) + 1 / N[5]) * (-(1 + p2p1))
}  # passed
d_enthalpy <- function(phi, temp, lattice, sigma) {
  p2p1 <- phi[2]/phi[1]
  -1 * 3/2 * get_alpha(temp, lattice) * (sum(sigma * phi) ^ 0.5) * (sigma[1] + sigma[2] * p2p1)
}  # passed
d_enthalpy_chi <- function(phi, chi){
  p2p1 <- phi[2]/phi[1]
  chi * (phi[5]*(1 + p2p1) - phi[1]*(1+p2p1)^2)
}
d_gibbs   <- function(phi, temp, N, lattice, sigma, chi) {
  p2p1           <- phi[2] / phi[1]
  d_entropy      <- d_entropy(phi, N)
  d_enthalpy     <- d_enthalpy(phi, temp, lattice, sigma)
  d_enthalpy_chi <- d_enthalpy_chi(phi, chi)
  d_gibbs        <- d_entropy + d_enthalpy + d_enthalpy_chi
}  # passed
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
# Fn --------------------
#' fn = fn(x), here x=phi1, other parameters are known
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

# Data preparation ---------------
# this part of code transform the experimental observations into volume fraction
# and temperature in Kelvin
get_expt  <- function() {
  library(dplyr)
  path_experiment <- '~/Desktop/dataset.csv'
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
  initial_guesses <- c(0.1, 0.5)
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
# Modeling ------------
simulate.fh <- function(){
  POLYNUM       <<- c(207, 2939,   1, 1, 1)
  SIGMA         <<- c(0, 0,   0, 0, 0)  # charge per monomer
  SIZE          <<- c(1, 1, 1, 1, 1) * 0.31E-9  # m; Voorn Overbeek 1967
  MW            <<- c(22e3, 900e3)
  MOLARRATIO    <<- c(2939,  11, 0.5, 0.5, 0)  # the polycation:polyanion and cation:anion molar ratio
  p2p1          <<- (POLYNUM[2] * MOLARRATIO[2] * SIZE[2]^3 ) /
    (POLYNUM[1] * MOLARRATIO[1] * SIZE[1]^3 )  # the volume fraction ratio between specie 1 and specie 2
  
  simulate('FH')
}
simulate.fh()