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
sigma.cr <- function(sigma0, temp){
  sigma <- ke^2 / (kEr*kkB*temp)
  return(min(c(sigma, sigma0)))
}
