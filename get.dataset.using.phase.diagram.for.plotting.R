library(dplyr)
# data =========
most.freq <- function(variable, n){
  names(sort(table(variable),decreasing=TRUE)[1:n])
}
get.dataset <- function(path = 'results/'){
  lapply(c('FH', 'FHVO', 'FHVOCR'), function(mod){
    list(readRDS(paste0(path, mod, '_fit.rds')) %>% 
           mutate(mod = mod, group = 'Fit'),
         readRDS(paste0('results/', mod, '.rds')) %>% 
           mutate(mod = mod, group = 'Expt')) %>% 
      bind_rows()
  }) %>% bind_rows() 
}
get.dataset.FHVO <- function(path = 'results/'){
  get.dataset(path) %>% filter(mod == 'FHVO')
}
saveds_ <- function(dataset, filename){
  write.csv(dataset, paste0(filename, '.csv'), row.names = F)
  saveRDS(dataset, paste0(filename, '.rds'))
}
save.dataset.FHVO <- function(path = 'results/'){
  get.dataset.FHVO(path) %>% 
    select(phi1, phi3, temp, group, mod) %>% 
    filter(phi3 %in% most.freq(phi3, 2)) %>% 
    saveds_(paste0(path, 'FHVO_fit_plot'))
}
save.dataset <- function(path = 'results/'){
  get.dataset(path) %>% 
    select(phi1, phi3, temp, group, mod) %>% 
    filter(phi3 %in% most.freq(phi3, 2)) %>% 
    saveds_(paste0(path, '_fit_plot'))
}
save.dataset.FHVO('results/')
save.dataset.FHVO('results_temp_phi3/')
save.dataset('results/')
save.dataset('results_temp_phi3/')
