# This file contains scripts to fit experimental data with Flory Huggins Chi para.
#       The definition of gibbs free energy is contained at vomodel.R
#

# Fn --------------------
#' Get chi value from given constrains
#'
#' A, B, C, D, E are five species defined in README
#'
#' @param phi1 A number. volume fraction of A in dilute phase
#' @param phi2 A number. volume fractino of A in dense phase
#' @param phi3 A number. volume fraction of C monovalent cation
#' @param p2p1 A number. the volume fraction ratio between specie 1 and specie 2
#' @inheritParams gibbs
#'
#' @return If inputs are valid, the output will be a length-one numeric vector.
#' @export
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
#' Phase separation function used for root finding
#'
#' \eqn{fn(x) = \frac{\partial G}{\partial x} \times (x - x0) - (G(x0) - G(x))}
#'
#' @inheritParams get_chi
#'
#' @return If inputs are valid, the output will be a length-one numeric vector
#' @export
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
#' Phase separation function (Counter-ion Release) used for root finding
#'
#' @inheritParams get_chi
#'
#' @return
#' @export
#'
#' @examples
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

# Simulate ------------------
#' Get Dense Phase \phi from Dilute Phase
#'
#' Find \phi1 given \phi2 and \phi3
#'
#' @param x A numberic vector of length 3. \phi2, \phi3 and T
#'
#' @return If root is found, the output will be a number \phi1.
#'         Otherwise it will be NA
#'
sim_      <- function(x, root_fn = fn){
  x <- as.numeric(x)
  phi2 <- x[1]
  phi3 <- x[2]
  temp <- x[3]
  initial_guesses <- INITIAL_GUESS
  phi1 <- NULL
  tryCatch({
    root <- uniroot(
      f = root_fn,
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
#' Simulate from Experiment
#'
#' @param output.filename.no.ext A string.
#' @param root_fn A callable function. Either fn or fn.cr.
#' @param expt A data.frame. Default is get
#' @param bWrite A boolean. Whether to write results to files or not.
#'
#' @return
#'
simulate  <- function(output.filename.no.ext = 'simulated',
                      root_fn = fn,
                      expt = get_expt(),
                      bWrite = TRUE){
  # get_chi(...) for dplyr::mutate
  get_chi_ <- function(phi1, phi2, phi3, temp){
    seq(1, length(phi1), 1) %>%
      sapply(function(i){
        get_chi(phi1 = phi1[i], phi2 = phi2[i], phi3 = phi3[i], temp = temp[i],
                N = POLYNUM, lattice = k.water.size, sigma = SIGMA, p2p1 = p2p1)
      })
  }
  expt <- get_expt()
  phi_dense <- apply(as.matrix(expt), 1, function(x) sim_(x, root_fn))
  if(is.null(phi_dense) |
     is.na(unique(phi_dense))){
    print('No fit found!')
  } else {
    output <- expt %>%
      mutate(phi2 = phi_dense) %>%
      mutate(chi = get_chi_(phi1, phi2, phi3, temp))

    if (bWrite) {
        output %>%
            write.csv(paste0(output.filename.no.ext, '.csv'), row.names = F)
        output %>%
            saveRDS(paste0(output.filename.no.ext, '.rds'))
    }

    return(output)
  }
}
# Modeling ------------
#' Flory Huggin Modeling
#'
#' @inheritParams simulate
#'
#' @return A data.frame.
#' @export
#'
simulate.fh <- function(expt = get_expt()){
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
#' Voorn-Overbeek Flory Huggins Modeling
#'
#' @inheritParams simulate
#'
#' @return A data.frame.
#' @export
#'
simulate.fhvo <- function(expt = get_expt()){
  POLYNUM       <<- c(207, 2939,   1, 1, 1)
  SIGMA         <<- c(11 / 207, 1,   1, 1, 0)  # charge per monomer
  SIZE          <<- c(1, 1, 1, 1, 1) * 0.31E-9  # m; Voorn Overbeek 1967
  MW            <<- c(22e3, 900e3)
  MOLARRATIO    <<- c(2939,  11, 0.5, 0.5, 0)  # the polycation:polyanion and cation:anion molar ratio
  p2p1          <<- (POLYNUM[2] * MOLARRATIO[2] * SIZE[2]^3 ) /
    (POLYNUM[1] * MOLARRATIO[1] * SIZE[1]^3 )  # the volume fraction ratio between specie 1 and specie 2
  INITIAL_GUESS <<- c(0.1, 0.5)
  simulate('results/FHVO', root_fn = fn, expt = get_expt())
}
#' Voorn-Overbeek Flory Huggins Counter-ion Release Modeling
#'
#' @inheritParams simulate
#'
#' @return A data.frame.
#' @export
#'
simulate.fhvocr <- function(expt = get_expt()) {
  POLYNUM       <<- c(207, 2939,   1, 1, 1)
  SIGMA         <<- c(11 / 207, 1,   1, 1, 0)  # charge per monomer
  SIZE          <<- c(1, 1, 1, 1, 1) * 0.31E-9  # m; Voorn Overbeek 1967
  MW            <<- c(22e3, 900e3)
  MOLARRATIO    <<- c(2939,  11, 0.5, 0.5, 0)  # the polycation:polyanion and cation:anion molar ratio
  p2p1          <<- (POLYNUM[2] * MOLARRATIO[2] * SIZE[2]^3 ) /
    (POLYNUM[1] * MOLARRATIO[1] * SIZE[1]^3 )  # the volume fraction ratio between specie 1 and specie 2
  INITIAL_GUESS <<- c(0.1, 0.5)
  simulate('results/FHVOCR', root_fn = fn.cr, expt = get_expt())
}

