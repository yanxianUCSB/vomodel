library(ggplot2)

# Plotting -------------
plot <- function(dataset){
  ggplot() +
    geom_point(aes(x = phi1, y = temp, color = phi3), 
               data = dataset %>% filter(group == 'Expt')) +
    # geom_point(aes(x = phi2, y = temp, color = phi3), 
    #            data = dataset %>% filter(group == 'Expt')) + 
    geom_line(aes(x = phi1, y = temp, color = phi3),
              data = dataset %>% filter(group == 'Fit')) +
    # geom_line(aes(x = phi2, y = temp, color = phi3),
    #           data = dataset %>% filter(group == 'Fit')) +
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
plot.mod <- function(){
  plot(get.dataset()) +
    facet_wrap(~mod) 
  ggsave('results/fitting.png', dpi = 1000, width = 5, height = 3)
  ggsave('results/fitting_wide.png', dpi = 1000, width = 8, height = 3)
}
plot.FHVO <- function(){
  g <- plot(get.dataset() %>% 
         filter(mod == 'FHVO'))
  ggsave('results/fitting_fhvo.png', dpi = 1000, width = 5, height = 3)
  saveRDS(g, 'results/fitting_fhvo.rds')
}
plot.mod()
plot.FHVO()