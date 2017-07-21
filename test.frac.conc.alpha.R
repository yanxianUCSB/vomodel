rm(list = ls())
source('vomodel.R')
source('para.spruijt.R')
test.frac.conc.alpha <- function(totl.conc, nacl.conc, Para) {
  print(Para)
  
  phis <-
    do.call(rbind, get.phis(lapply(totl.conc, get.conc, nacl.conc, Para), Para))
  print(paste('tot.conc =', totl.conc, 'mg/mL'))
  print(paste('nacl.conc =', nacl.conc, 'mM'))
  print(paste('>>protein.frac =', phis[, 1]))
  print(paste('>>rna.frac =', phis[, 2]))
  print(paste('>>nacl.frac =', phis[, 3]))
  print(paste('>>nacl.frac =', phis[, 4]))
  print(paste('>>water.frac =', phis[, 5]))
  
  
  print(paste('conc =', get.conc(totl.conc, nacl.conc, Para), 'mol/m^3'))
  print(paste(
    'alpha =',
    alpha(temp, Para$size),
    'J',
    '[REF = 3.655 in Voorn M, Complex Coacervation, PhD Thesis]'
  ))
  
}

totl.conc <- 2
nacl.conc <- 0
temp <- 300
test.frac.conc.alpha(totl.conc, nacl.conc, Para)
