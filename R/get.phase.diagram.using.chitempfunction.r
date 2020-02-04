library(dplyr)
# loading ---------

load.chitemp <- function(){
  readRDS('results/chi_temp_functions.rds')
}
load.expt <- function(){
  readRDS('expt.rds')
}
# Fn --------------------
#' fn = fn(x), here x=phi1, other parameters are known
fn        <- function(x, phi3, temp, N, lattice, sigma, p2p1, chitemp){
  phi1 <- x[1]
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
fn.cr        <- function(x, phi3, temp, N, lattice, sigma, p2p1, chitemp){
  sigma         <- c(0, 0,   0, 0, 0)  # charge per monomer
  sigma[1]      <- sigma.cr(sigma[1], temp)
  sigma[2]      <- sigma.cr(sigma[2], temp)
  fn(x, phi3, temp, N, lattice, sigma, p2p1, chitemp)
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
find.root <- function(phi3, temp, N, lattice, sigma, p2p1, fn, chitemp) {
  roots <- nleqslv::nleqslv (
    x = INITIAL_GUESS,
    fn=fn,
    phi3 = phi3, temp=temp, N=N, lattice=lattice, sigma=sigma, p2p1=p2p1,
    chitemp=chitemp,
    method = 'Newton',
    global = 'cline',
    control = list(xtol = 1e-10, ftol = 1e-10)
  )
  roots$x
}
find.root.here <- function(phi3, temp, sigma, fn, chitemp){
  POLYNUM       <<- c(207, 2939,   1, 1, 1)
  SIZE          <<- c(1, 1, 1, 1, 1) * 0.31E-9  # m; Voorn Overbeek 1967
  MW            <<- c(22e3, 900e3)
  MOLARRATIO    <<- c(2939,  11, 0.5, 0.5, 0)  # the polycation:polyanion and cation:anion molar ratio
  p2p1          <<- (POLYNUM[2] * MOLARRATIO[2] * SIZE[2]^3 ) /
    (POLYNUM[1] * MOLARRATIO[1] * SIZE[1]^3 )  # the volume fraction ratio between specie 1 and specie 2
  find.root(phi3, temp, N=POLYNUM, lattice = 0.31E-9, sigma = SIGMA,
            p2p1 = p2p1, fn=fn, chitemp=chitemp)
}
find.root.here.FH <- function(phi3, temp, chitempfun) {
  SIGMA         <<- c(0, 0,   0, 0, 0)  # charge per monomer
  INITIAL_GUESS <<- INITIAL_GUESS
  find.root.here(phi3,temp,sigma=SIGMA, fn=fn, chitemp=chitempfun)
}
find.root.here.FHVO <- function(phi3, temp, chitempfun) {
  SIGMA         <<- c(11/207, 1,   1, 1, 0)  # charge per monomer
  INITIAL_GUESS <<- INITIAL_GUESS
  find.root.here(phi3,temp,sigma=SIGMA, fn=fn, chitemp=chitempfun)
}
find.root.here.FHVO1salt <- function(phi3, temp, chitempfun) {
  SIGMA         <<- c(11/207, 1,   1, 1, 0)  # charge per monomer
  INITIAL_GUESS <<- INITIAL_GUESS
  chitempfun2 <- function(temp){chitempfun(phi3, temp)}
  find.root.here(phi3,temp,sigma=SIGMA, fn=fn, chitemp=chitempfun2)
}
find.root.here.FHVOCR <- function(phi3, temp, chitempfun) {
  SIGMA         <<- c(11/207, 1,   1, 1, 0)  # charge per monomer
  INITIAL_GUESS <<- INITIAL_GUESS
  find.root.here(phi3,temp,sigma=SIGMA, fn=fn.cr, chitemp=chitempfun)
}

# Simulate ===========
sim_ <- function(phi3s, temps, find.root.func, chitempfun) {
  lapply(phi3s, function(phi3){
    lapply(temps, function(temp){
    x12 <- find.root.func(phi3 = phi3, temp=temp, chitempfun)
    data.frame(phi1 = x12[1],
               phi2 = x12[2],
               phi3 = phi3,
               temp = temp,
               chi = chitempfun(temp))
  }) %>% bind_rows()
  }) %>% bind_rows()
}
simulate.FH <- function(phi3s, temps, chitempfun) {
  sim_(phi3s, temps, find.root.here.FH, chitempfun)
}
simulate.FHVO <- function(phi3s, temps, chitempfun) {
  sim_(phi3s, temps, find.root.here.FHVO, chitempfun)
}
simulate.FHVOCR <- function(phi3s, temps, chitempfun) {
  sim_(phi3s, temps, find.root.here.FHVOCR, chitempfun)
}
simulate.FHVO1salt <- function(phi3s, temps, chitempfun) {
  sim_(phi3s, temps, find.root.here.FHVO1salt, chitempfun)
}
simulate <- function() {
  INITIAL_GUESS <<- c(0.001, 0.35)
  expt <- load.expt()
  chitempfunlist <- load.chitemp()
  # phi3s <- as.numeric(names(sort(table(expt$phi3),decreasing=TRUE)[1:2]))
  phi3s <- unique(expt$phi3)
  temps <- seq(min(expt$temp), max(expt$temp), (max(expt$temp)-min(expt$temp))/10)
  fh <- simulate.FH(phi3s, temps, chitempfunlist$FH)
  write.csv(fh, 'results/FH_fit.csv', row.names = F)
  saveRDS(fh, 'results/FH_fit.rds')
  fhvo <- simulate.FHVO(phi3s, temps, chitempfunlist$FHVO)
  write.csv(fhvo, 'results/FHVO_fit.csv', row.names = F)
  saveRDS(fhvo, 'results/FHVO_fit.rds')
  # simulate.FHVO1salt(phi3s, temps, chitempfunlist$FHVO.1salt)
  fhvocr <- simulate.FHVOCR(phi3s, temps, chitempfunlist$FHVOCR)
  write.csv(fhvocr, 'results/FHVOCR_fit.csv', row.names = F)
  saveRDS(fhvocr, 'results/FHVOCR_fit.rds')
}
# 1 example ======
simulate.1.eg <- function(){
  INITIAL_GUESS <<- c(0.001, 0.35)
  ds <- readRDS('results/FHVO.rds')
  eg <- ds[c(21, 26),]
  chitempfun <- function(temp){return(eg$chi[which(eg$temp == temp)])}
  phi3s <- unique(ds$phi3)
  phi3s <- seq(0.1 * min(phi3s), 5.4*max(phi3s), (max(phi3s)-min(phi3s))/20)
  temps <- eg$temp
  fhvo <- simulate.FHVO(phi3s, temps, chitempfun)
  fhvo <- bind_rows(eg, fhvo)
  write.csv(fhvo, 'results/FHVO_fit_1_example.csv', row.names = F)
  saveRDS(fhvo, 'results/FHVO_fit_1_example.rds')
}

