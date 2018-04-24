source('vomodel.r')

# Parameters
POLYNUM       <<- c(1000, 1000,   1, 1, 1)
SIGMA         <<- c(0.15, 0.15,   1, 1, 0)  # charge per monomer
SIZE          <<- c(1, 1, 1, 1, 1) * 0.31E-9  # m; Voorn Overbeek 1967
MW            <<- c(1, 1)
MOLARRATIO    <<- c(1,  1, 0.5, 0.5, 0)  # the polycation:polyanion and cation:anion molar ratio
p2p1          <<- (POLYNUM[2] * MOLARRATIO[2] * SIZE[2]^3 ) / 
  (POLYNUM[1] * MOLARRATIO[1] * SIZE[1]^3 )  # the volume fraction ratio between specie 1 and specie 2

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

# get_chi
get_chi   <- function(phi1, phi2, phi3, temp, N, lattice, sigma, p2p1) {
  phi_1 <-
    c(phi1, p2p1 * phi1, phi3, phi3, 1 - phi1 - p2p1 * phi1 - 2 * phi3)
  d_entropy_1 <-
    1 / N[1] * log(phi_1[1]) + 1 / N[1] + (1 / N[2] * log(phi_1[2]) + 1 / N[2]) * p2p1
  d_enthalpy_1 <- 
    3/2 * get_alpha(temp, lattice) * (sum(sigma * phi_1) ^ 0.5) * (sigma[1] + sigma[2] * p2p1)
  under_1 <- 
    phi_1[5] * (p2p1 + 1) + (phi_1[1] + phi_1[2]) * p2p1
  
  phi_2 <-
    c(phi2, p2p1 * phi2, phi3, phi3, 1 - phi2 - p2p1 * phi2 - 2 * phi3)
  d_entropy_2 <-
    1 / N[1] * log(phi_2[1]) + 1 / N[1] + (1 / N[2] * log(phi_2[2]) + 1 / N[2]) * p2p1
  d_enthalpy_2 <- 
    3/2 * get_alpha(temp, lattice) * (sum(sigma * phi_2) ^ 0.5) * (sigma[1] + sigma[2] * p2p1)
  under_2 <- 
    phi_2[5] * (p2p1 + 1) + (phi_2[1] + phi_2[2]) * p2p1
  
  chi <-
    ((d_entropy_1 - d_entropy_2) + 
       (d_enthalpy_1 - d_enthalpy_2)) / (under_2 - under_1)
}
test.get_chi <- function(){
  get_ref() %>%
    mutate(group = (phi1 - phi1[which.max(phi3)])<=0) ->
    ds
  dsref <- data.frame(
    phi1 = ds$phi1[ds$group],
    phi2 = rev(ds$phi1[!ds$group]),
    phi3 = ds$phi3[!ds$group],
    temp = ds$temp[ds$group]
  )
  apply(dsref %>% as.matrix(), 1, function(x){
    get_chi(
      phi1 = x[1],
      phi2 = x[2],
      phi3 = x[3],
      temp = x[4],
      N       = POLYNUM,
      lattice = k.water.size,
      sigma   = SIGMA,
      p2p1    = p2p1
    )
  }) -> chi
  return(
    dsref %>% mutate(chi = chi)
  )
}

fn        <- function(phi1, phi2, phi3, temp, N, lattice, sigma, p2p1){
  chi <- 0
  phi_1 <-
    c(phi1, p2p1 * phi1, phi3, phi3, 1 - phi1 - p2p1 * phi1 - 2 * phi3)
  phi_2 <-
    c(phi2, p2p1 * phi2, phi3, phi3, 1 - phi2 - p2p1 * phi2 - 2 * phi3)
  diff <- 
    d_gibbs(phi_2, temp, N, lattice, sigma, chi, p2p1) * (phi2 - phi1) -
    gibbs(phi_2, temp, N, lattice, sigma, chi) +
    gibbs(phi_1, temp, N, lattice, sigma, chi)
}
sim_ <- function(x){
  x <- as.numeric(x)
  phi2 <- x[1]
  phi3 <- x[2]
  temp <- x[3]
  initial_guesses <- c(0.01, 0.05)
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
# Thu Apr 19 10:15:02 2018  Test ---------------------- 
test <- function(){
  library(ggplot2)
  # example <- as.matrix(get_expt())[2, 1:3]
  example <- as.numeric(as.matrix(get_ref())[2, 1:3])
  print(example)
  phi2 <- example[1]
  phi3 <- example[2]
  temp <- example[3]
  initial_guesses <- c(0.005, 0.01)
  x <- seq(initial_guesses[1], initial_guesses[2], 1e-3)
  y <- sapply(x, function(x){
    fn(
      phi1 = x,
      phi2 = phi2,
      phi3 = phi3,
      temp = temp,
      N       = POLYNUM,
      lattice = k.water.size,
      sigma   = SIGMA,
      p2p1    = p2p1
    )
  }, simplify = T)
  
  ggplot(data.frame(x,y)) +
    geom_line(aes(x,y))+
    geom_hline(yintercept=0)+
    scale_x_log10()
}
test2 <- function(){
  ref <- get_ref()
  phi_dense <- apply(ref, 1, function(x) sim_(x))
  if(F){
    print('No fit found!')
  } else {
    ref %>% 
      mutate(phi2 = phi_dense) ->
      ds
    ds %>% 
      write.csv('simulated.csv', row.names = F)
    return(ds)
  }
}
# simulate()
testds <- test2() %>% select(-temp)
print(testds)
plot(testds)

# plot(get_ref())
