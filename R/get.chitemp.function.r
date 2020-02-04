most.frequent <- function(variable){
  names(sort(table(variable),decreasing=TRUE)[1])
}
get.chitemp <- function(file) {
  readRDS(paste0('results/', file, '.rds')) -> ds
  chi <- ds$chi
  temp <- ds$temp
  return(data.frame(chi=chi, temp=temp))
}
get.chitemp.1salt <- function(file) {
  readRDS(paste0('results/', file, '.rds')) -> ds
  ds <- ds[which(ds$phi3 == most.frequent(ds$phi3)),]
  chi <- ds$chi
  temp <- ds$temp
  return(data.frame(chi=chi, temp=temp))
}
mod.chitemp <- function(temp, chi){
  mod <- lm(y ~ x, data=data.frame(
    x = 1/temp,
    y = chi
  ))
}
fit.chitemp <- function(mod){
  return(function(temp){
    predict.lm(mod, newdata = data.frame(x = 1/temp))
  })
}

get.chitemp.function <- function(file) {
  chitemp <- get.chitemp(file)
  mod.chitemp <- mod.chitemp(chitemp$temp, chitemp$chi)
  return(fit.chitemp(mod.chitemp))
}
get.chitemp.function.1salt <- function(file) {
  chitemp <- get.chitemp.1salt(file)
  mod.chitemp <- mod.chitemp(chitemp$temp, chitemp$chi)
  return(fit.chitemp(mod.chitemp))
}

FH <- get.chitemp.function('FH')
FHVO <- get.chitemp.function('FHVO')
FHVO.1salt <- get.chitemp.function.1salt('FHVO')
FHVOCR <- get.chitemp.function('FHVOCR')

saveRDS(list(
  FH=FH,
  FHVO=FHVO,
  FHVO.1salt=FHVO.1salt,
  FHVOCR=FHVOCR
), 'results/chi_temp_functions.rds')

