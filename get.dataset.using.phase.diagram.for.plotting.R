library(dplyr)
# data =========
most.freq <- function(variable, n){
  names(sort(table(variable),decreasing=TRUE)[1:n])
}
get.dataset <- function(path = 'results/'){
  lapply(c('FH', 'FHVO', 'FHVOCR'), function(mod){
    list(readRDS(paste0(path, mod, '_fit.rds')) %>% 
           mutate(mod = mod, group = 'Fit'),
         readRDS(paste0(path, mod, '.rds')) %>% 
           mutate(mod = mod, group = 'Expt')) %>% 
      bind_rows()
  }) %>% bind_rows() %>% 
    filter(
      phi3 %in% most.freq(phi3, 2)
    ) %>% 
    mutate(
      phi3 = as.character(phi3)
      # phi3 = paste("\u03d53 =", as.character(round(phi3 * 1000, 2)))
    ) 
}
get.dataset.FHVO <- function(path = 'results/'){
  get.dataset(path) %>% 
    filter(mod == 'FHVO')
}
save.dataset <- function(dataset, filename){
  write.csv(dataset, paste0(filename, '.csv'), row.names = F)
  saveRDS(dataset, paste0(filename, '.rds'))
}
save.dataset.FHVO <- function(path = 'results/'){
  get.dataset.FHVO(path) %>% 
    select(x = phi1,
           y = temp,
           group = group,
           group2 = phi3) %>% 
    save.dataset(paste0(path, 'FHVO_fit_plot'))
}
save.dataset.FHVO()