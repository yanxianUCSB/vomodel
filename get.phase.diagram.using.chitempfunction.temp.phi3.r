library(dplyr)
source('vomodel.R')
# loading ---------

load.chitemp <- function(){
  readRDS('results/chi_temp_functions.rds')
}
load.expt <- function(){
  readRDS('expt.rds')
}
# Fn --------------------
#' fn = fn(x), here x=phi1, other parameters are known
fn        <- function(x, phi1, temp, N, lattice, sigma, p2p1, chitemp){
  phi3 <- x[1]
  phi2 <- x[2]
  phi_1 <-
    c(phi1, p2p1 * phi1, phi3, phi3, 1 - phi1 - p2p1 * phi1 - 2 * phi3)
  phi_2 <-
    c(phi2, p2p1 * phi2, phi3, phi3, 1 - phi2 - p2p1 * phi2 - 2 * phi3)
  chi <- chitemp(temp)
  diff.1 <- 
    d_gibbs(phi_2, temp, N, lattice, sigma, chi) * (phi2 - phi1) -
    gibbs(phi_2, temp, N, lattice, sigma, chi) +
    gibbs(phi_1, temp, N, lattice, sigma, chi)
  diff.2 <- 
    d_gibbs(phi_2, temp, N, lattice, sigma, chi) -
    d_gibbs(phi_1, temp, N, lattice, sigma, chi) 
  return(c(diff.1, diff.2))
}
fn.cr        <- function(x, phi1, temp, N, lattice, sigma, p2p1, chitemp){
  sigma         <- c(0, 0,   0, 0, 0)  # charge per monomer
  sigma[1]      <- sigma.cr(sigma[1], temp)
  sigma[2]      <- sigma.cr(sigma[2], temp)
  fn(x, phi1, temp, N, lattice, sigma, p2p1, chitemp)
}
# fn.cr        <- function(phi1, phi2, phi3, temp, N, lattice, sigma, p2p1, chitemp){
#   sigma[1] <- sigma.cr(sigma[1], temp)
#   sigma[2] <- sigma.cr(sigma[2], temp)
#   phi_1    <- c(phi1, p2p1 * phi1, phi3, phi3, 1 - phi1 - p2p1 * phi1 - 2 * phi3)
#   phi_2    <- c(phi2, p2p1 * phi2, phi3, phi3, 1 - phi2 - p2p1 * phi2 - 2 * phi3)
#   chi <- chitemp(temp)
#   diff.1     <- 
#     d_gibbs(phi_2, temp, N, lattice, sigma, chi) * (phi2 - phi1) -
#     gibbs(phi_2, temp, N, lattice, sigma, chi) +
#     gibbs(phi_1, temp, N, lattice, sigma, chi)
#   diff.2 <- 
#     d_gibbs(phi_2, temp, N, lattice, sigma, chi) -
#     d_gibbs(phi_1, temp, N, lattice, sigma, chi) 
#   return(c(diff.1, diff.2))
# }
# Roots --------------
find.root <- function(phi1, temp, N, lattice, sigma, p2p1, fn, chitemp) {
  roots <- nleqslv::nleqslv (
    x = INITIAL_GUESS,
    fn=fn,
    phi1 = phi1, temp=temp, N=N, lattice=lattice, sigma=sigma, p2p1=p2p1, 
    chitemp=chitemp,
    method = 'Newton',
    global = 'cline',
    control = list(xtol = 1e-10, ftol = 1e-10)
  )
  roots$x
} 
find.root.here <- function(phi1, temp, sigma, fn, chitemp){
  POLYNUM       <<- c(207, 2939,   1, 1, 1)
  SIZE          <<- c(1, 1, 1, 1, 1) * 0.31E-9  # m; Voorn Overbeek 1967
  MW            <<- c(22e3, 900e3)
  MOLARRATIO    <<- c(2939,  11, 0.5, 0.5, 0)  # the polycation:polyanion and cation:anion molar ratio
  p2p1          <<- (POLYNUM[2] * MOLARRATIO[2] * SIZE[2]^3 ) /
    (POLYNUM[1] * MOLARRATIO[1] * SIZE[1]^3 )  # the volume fraction ratio between specie 1 and specie 2
  find.root(phi1, temp, N=POLYNUM, lattice = 0.31E-9, sigma = SIGMA, 
            p2p1 = p2p1, fn=fn, chitemp=chitemp)
}
find.root.here.FH <- function(phi1, temp, chitempfun) {
  SIGMA         <<- c(0, 0,   0, 0, 0)  # charge per monomer
  INITIAL_GUESS <<- INITIAL_GUESS
  find.root.here(phi1,temp,sigma=SIGMA, fn=fn, chitemp=chitempfun)
}
find.root.here.FHVO <- function(phi1, temp, chitempfun) {
  SIGMA         <<- c(11/207, 1,   1, 1, 0)  # charge per monomer
  INITIAL_GUESS <<- INITIAL_GUESS
  find.root.here(phi1,temp,sigma=SIGMA, fn=fn, chitemp=chitempfun)
}
find.root.here.FHVOCR <- function(phi1, temp, chitempfun) {
  SIGMA         <<- c(11/207, 1,   1, 1, 0)  # charge per monomer
  INITIAL_GUESS <<- INITIAL_GUESS
  find.root.here(phi1,temp,sigma=SIGMA, fn=fn.cr, chitemp=chitempfun)
}

# Simulate ===========
sim_ <- function(phi1s, temps, find.root.func, chitempfun) {
  lapply(phi1s, function(phi1){
    lapply(temps, function(temp){
      x12 <- find.root.func(phi1 = phi1, temp=temp, chitempfun)
      data.frame(phi3 = x12[1], 
                 phi2 = x12[2],
                 phi1 = phi1,
                 temp = temp,
                 chi = chitempfun(temp))
    }) %>% bind_rows() 
  }) %>% bind_rows()
}
simulate.FH <- function(phi1s, temps, chitempfun) {
  sim_(phi1s, temps, find.root.here.FH, chitempfun)
}
simulate.FHVO <- function(phi1s, temps, chitempfun) {
  sim_(phi1s, temps, find.root.here.FHVO, chitempfun)
}
simulate.FHVOCR <- function(phi1s, temps, chitempfun) {
  sim_(phi1s, temps, find.root.here.FHVOCR, chitempfun)
}
simulate.FHVO1salt <- function(phi1s, temps, chitempfun) {
  sim_(phi1s, temps, find.root.here.FHVO1salt, chitempfun)
}
simulate <- function() {
  INITIAL_GUESS <<- c(0.0008, 0.4)
  expt <- load.expt()
  chitempfunlist <- load.chitemp()
  # phi1s <- as.numeric(names(sort(table(expt$phi1),decreasing=TRUE)[1:2]))
  phi1s <- unique(expt$phi1)
  temps <- seq(min(expt$temp), max(expt$temp), (max(expt$temp)-min(expt$temp))/10)
  fh <- simulate.FH(phi1s, temps, chitempfunlist$FH) 
  write.csv(fh, 'results_temp_phi3/FH_fit.csv', row.names = F)  
  saveRDS(fh, 'results_temp_phi3/FH_fit.rds')
  fhvo <- simulate.FHVO(phi1s, temps, chitempfunlist$FHVO)
  write.csv(fhvo, 'results_temp_phi3/FHVO_fit.csv', row.names = F)  
  saveRDS(fhvo, 'results_temp_phi3/FHVO_fit.rds')
  # simulate.FHVO1salt(phi1s, temps, chitempfunlist$FHVO.1salt)
  fhvocr <- simulate.FHVOCR(phi1s, temps, chitempfunlist$FHVOCR)
  write.csv(fhvocr, 'results_temp_phi3/FHVOCR_fit.csv', row.names = F)  
  saveRDS(fhvocr, 'results_temp_phi3/FHVOCR_fit.rds')
}
simulate()