source('vomodel.r')

# Parameters
POLYNUM       <<- c(100, 100,   1, 1, 1)
SIGMA         <<- c(1 / 1, 1,   1, 1, 0)  # charge per monomer
SIZE          <<- c(1, 1, 1, 1, 1) * 0.31E-9  # m; Voorn Overbeek 1967
MW            <<- c(1, 1)
MOLARRATIO    <<- c(1,  1, 0.5, 0.5, 0)  # the polycation:polyanion and cation:anion molar ratio
p2p1          <<- (POLYNUM[2] * MOLARRATIO[2] * SIZE[2]^3 ) / 
  (POLYNUM[1] * MOLARRATIO[1] * SIZE[1]^3 )  # the volume fraction ratio between specie 1 and specie 2

get_ref <- function(){
  library(dplyr)
  path <- 'ref.chimad.txt'
  read.csv(path, sep = '\t') %>% 
    select(phi1 = VolumeFraction, phi3=Binodal) %>% 
    mutate(temp = 298.15) ->
    ds
  ds %>% saveRDS('ref.rds')
  return(ds)
}
# Thu Apr 19 10:15:02 2018  Test ---------------------- 
test <- function(){
  library(ggplot2)
  # example <- as.matrix(get_expt())[2, 1:3]
  example <- as.matrix(get_ref())[2, 1:3]
  print(example)
  phi2 <- example[1]
  phi3 <- example[2]
  temp <- example[3]
  initial_guesses <- phi2 + (1-phi2) * c(0, 0.5)
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
  phi_dense <- apply(as.matrix(ref), 1, sim_)
  if(is.null(phi_dense) | 
     is.na(unique(phi_dense))){
    print('No fit found!')
  } else {
    expt %>% 
      mutate(phi2 = phi_dense) %>% 
      write.csv('simulated.csv', row.names = F)
  }
}
# simulate()
test()