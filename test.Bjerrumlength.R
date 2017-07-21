source('vomodel.R')
test.Bjerrumlength <- function(temp = 300){
  # Bjerrum length at 300K ~0.7nm
  lB <- ke ^ 2 / (kEr * kkB * temp)
  print(paste('Bjerrum length', temp, 'K:', lB / 1E-9, 'nm'))
}