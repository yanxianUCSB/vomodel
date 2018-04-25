rm(list = ls())
source('vomodel.r')
library(dplyr)
# Parameters --------------
POLYNUM       <<- c(200, 3000,   1, 1, 1)
SIGMA         <<- c(11/207, .3, 1, 1, 0)  # charge per monomer
SIZE          <<- c(1, 1, 1, 1, 1) * 0.31E-9  # m; Voorn Overbeek 1967
# MW            <<- c(1, 1)
MOLARRATIO    <<- c(3000,  10, 1, 1, 0)  # the polycation:polyanion and cation:anion molar ratio
p2p1          <<- (POLYNUM[2] * MOLARRATIO[2] * SIZE[2]^3 ) / 
  (POLYNUM[1] * MOLARRATIO[1] * SIZE[1]^3 )  # the volume fraction ratio between specie 1 and specie 2

# Functions ------------
phi12phi <- function(phi1, phi3=5e-2){
  phi <- c(phi1, phi1, phi3, phi3)
  phi[5] <- 1 - sum(phi)
  return(phi)
}
test.get.phi2.when.chi.is.0 <- function(phi1) {
  temp <- 300
  N <- POLYNUM
  lattice <- 0.31E-9
  sigma <- SIGMA
  chi <- 0
  x <- seq(1e-3, 0.4, 1e-2)
  x %>% 
    sapply(function(phi2){
      get_chi(phi1, phi2, phi3 = 0.05, temp, N, lattice, sigma, p2p1 = 1) %>% 
        return()
    }) -> 
    chis
  y <- spline(x=chis, y=x, xout=0)$y
  return(y)
}  
test.get.phi2.when.chi.is.0.with.config <- function(phi1, phi3, temp = 300, 
                                          N=POLYNUM, lattice=0.31E-9, 
                                          sigma=SIGMA, chi=0) {
  x <- seq(1e-6, 0.2, 1e-4)
  x %>% 
    sapply(function(phi2){
      get_chi(phi1, phi2, phi3, temp, N, lattice, sigma, p2p1 = 1) %>% 
        return()
    }) -> 
    chis
  return(spline(x=chis, y=x, xout=chi)$y)
}  
# Get ref dataset ---------------
get_ref <- function(){
  library(dplyr)
  path <- 'ref.chimad.txt'
  read.csv(path, sep = '\t') %>% 
    mutate(phi1 = VolumeFraction / 2, phi3=Binodal / 2) %>% 
    mutate(temp = 298.15) %>% 
    arrange(phi1) %>% 
    select(phi1, phi3, temp) ->
    ds
  ds %>% saveRDS('ref.rds')
  return(ds)
}
get_ref_fh_chi_is_0 <- function(){
  library(dplyr)
  path <- 'ref.chimad.fh.chi.is.0.txt'
  read.csv(path, sep = '\t') %>% 
    mutate(phi1 = VolumeFraction / 2, phi3=Binodal / 2) %>% 
    mutate(temp = 298.15) %>% 
    arrange(phi1) %>% 
    select(phi1, phi3, temp) ->
    ds
  ds %>% saveRDS('ref.rds')
  return(ds)
}
# Test cases -------------
test.case.1 <- function() {
  # test gibbs
  sapply(seq(1e-3, 1, 1e-3), function(phi1){
    # temp <- 300
    # N <- POLYNUM
    # lattice <- 0.31E-9
    # sigma <- SIGMA
    # chi <- 0
    phi <- c(phi1/2, phi1/2, 5e-2, 5e-2)
    phi[5] <- 1 - sum(phi)
    gibbs(phi = phi, temp = 303, N = POLYNUM, lattice = 0.31E-9, sigma = SIGMA, chi = 0) %>%
      return()
  }, simplify = T) %>% 
    plot()
  print('curve should be monotonic reducing')
}  # test gibbs
# test.case.1()
test.case.3 <- function() {
  temp <- 300
  N <- POLYNUM
  lattice <- 0.31E-9
  sigma <- SIGMA
  chi <- 0
  seq(1e-3, 0.4, 1e-3) %>% 
    sapply(function(phi1){
      da1 <- numDeriv::grad(func = function(phi1){
        phi <- c(phi1, phi1, 5e-2, 5e-2)
        phi[5] <- 1 - sum(phi)
        enthalpy <- -1 * get_alpha(temp, lattice) * sum(sigma * phi) ^ 1.5
      }, x=phi1)
      phi <- c(phi1, phi1, 5e-2, 5e-2)
      phi[5] <- 1 - sum(phi)
      da2 <- -1 * 3/2 * get_alpha(temp, lattice) * (sum(sigma * phi) ^ 0.5) * (sigma[1] + sigma[2] * p2p1)
      da1 - da2
    }) %>% 
    plot
}
# test.case.3()  # d_enthalpy passes
test.case.4 <- function() {
  temp <- 300
  N <- POLYNUM
  lattice <- 0.31E-9
  sigma <- SIGMA
  chi <- 0
  seq(1e-3, 0.1, 1e-3) %>% 
    sapply(function(phi1){
      da1 <- numDeriv::grad(x=phi1, func = function(phi1){
        phi <- phi12phi(phi1)
        entropy <- sum(phi / N * log(phi))
      })
      
      phi <- phi12phi(phi1)
      
      da2 <-     
        (1 / N[1] * log(phi[1]) + 1 / N[1]) + 
        (1 / N[2] * log(phi[2]) + 1 / N[2]) * p2p1 +
        (1 / N[5] * log(phi[5]) + 1 / N[5]) * (-(1 + p2p1))
      
      da1 - da2
    }) %>% 
    plot
}
# test.case.4()  # d_entropy passes
test.case.2 <- function() {
  sapply(seq(1e-3, 0.4, 1e-3), function(phi1){
    phi <- c(phi1, phi1, 5e-2, 5e-2)
    phi[5] <- 1 - sum(phi)
    g <- d_gibbs(phi = phi, temp = 303, N = POLYNUM, lattice = 0.31E-9, 
                 sigma = SIGMA, chi = 0) 
    g2 <- numDeriv::grad(func = function(x){
      phi <- c(x, x, 5e-2, 5e-2)
      phi[5] <- 1 - sum(phi)
      gibbs(phi = phi, temp = 303, N = POLYNUM, lattice = 0.31E-9, sigma = SIGMA, chi = 0)
    }, x = phi1)
    return(g-g2)
  }, simplify = T) %>% 
    plot()
}
# test.case.2()  # d_gibbs passes
test.case.5 <- function() {
  temp <- 300
  N <- POLYNUM
  lattice <- 0.31E-9
  sigma <- SIGMA
  chi <- 0
  seq(1e-6, 0.4, 1e-3) %>% 
    sapply(function(phi1){
      
      phi2 <- .2
      phi3 <- 0
      get_chi(phi1, phi2, phi3 = phi3, temp, N, lattice, sigma, p2p1 = 1) %>% 
        return()
      
    }) %>% 
    plot(x=seq(1e-3, 0.4, 1e-3), ylab = 'chi', xlab = 'phi1')
}  # get_chi
# test.case.5()  # one case
test.case.6 <- function() {
  temp <- 300
  N <- POLYNUM
  lattice <- 0.31E-9
  sigma <- SIGMA
  chi <- 0
  phi1s <- seq(1e-6, 0.4, 1e-3)
  phi1s %>% 
    sapply(function(phi1){
      test.get.phi2.when.chi.is.0(phi1)
    }) ->
    phi2s
  plot(phi1s, phi2s)
}  # get phi2
# test.case.6()  # one binodal curve for chi = 0
test.case.7 <- function() {
  SIGMA         <<- c(0, 0,   0, 0, 0)  # No VO term
  ref <- get_ref_fh_chi_is_0()
  1:nrow(ref) %>% 
    sapply(function(row){
      phi1 <- ref$phi1[row]
      # phi3 <- ref$phi3[row]
      phi3 <- 0
      temp <- ref$temp[row]
      test.get.phi2.when.chi.is.0.with.config(phi1, phi3, temp, 
                                              N=POLYNUM, lattice=0.31E-9, 
                                              sigma=SIGMA, chi=0) %>% 
        return()
    }) ->
    phi2
  ref$phi2 <- phi2
  plot(c(ref$phi1, ref$phi2), c(ref$phi3, ref$phi3))
  return(ref)
}  # get_chi
# ref <- test.case.7()  # one binodal curve for chi = 0
# ref
# Test fn -------------------
# phi3 = 0.01
# temp = 300
# sigma = 0 (FH, no VO)
# phi1, phi2 =  0.00155, 0.050
#               0.00455, 0.035
# phi3 = 0.0004504
# phi1, phi2 =  7.458e-05, 0.32
# 
test.fn <- function() {
  phi1 <- 7.458e-05
  phi3 <- 0.0004504
  temp <- 288
  N <- POLYNUM
  lattice <- 0.31E-9
  sigma <- SIGMA
  p2p1 <- p2p1
  phi2s <- seq(1e-5, 0.3, 1e-3)
  phi2s %>% 
    sapply(function(phi2){
      fn(phi1, phi2, phi3, temp, N, lattice, sigma, p2p1)
    }, simplify = T) ->
    fn
  plot(phi2s, fn)
}
test.fn()

