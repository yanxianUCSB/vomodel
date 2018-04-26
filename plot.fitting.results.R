library(ggplot2)
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
      phi3 = paste("\u03d53 =", as.character(round(phi3 * 1000, 2)))
    ) 
}
# Plotting -------------
plot <- function(){
  dataset <- get.dataset()
  ggplot() +
    geom_point(aes(x = phi1, y = temp, color = 'Expt', pch = phi3), 
               data = dataset %>% filter(group == 'Expt')) +
    geom_point(aes(x = phi2, y = temp, color = 'Expt', pch = phi3), 
               data = dataset %>% filter(group == 'Expt')) + 
    geom_line(aes(x = phi1, y = temp, color = 'Fit'),
              data = dataset %>% filter(group == 'Fit')) +
    geom_line(aes(x = phi2, y = temp, color = 'Fit'),
              data = dataset %>% filter(group == 'Fit')) +
    facet_wrap(~mod) +
    scale_x_log10()+
    labs(
      x = expression(paste(phi["1"])),
      y = expression(paste("T (K)")),
      pch = NULL,
      color = NULL
    ) +
    theme_bw()+
    ggsci::scale_color_npg()
}
plot()
ggsave('results/fitting.png', dpi = 300, width = 5, height = 3)