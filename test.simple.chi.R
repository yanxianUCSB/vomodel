rm(list = ls())
source('vomodel.R')

g.simple <- function(x, y, chi) {
    x * 0.002 * log( 0.5 * x) + 2 * y * log(0.5*y) + (1-x-y)*log(1-x-y) -
        3.622*(0.075*x + y)^1.5 + 0.5 * x^2 * chi
}
dg.simple <- function(x, y, chi) {
    0.002 * log(0.5*x) + 0.002 + 2*y*log(0.5*y) + (1-x-y)*log(1-x-y)-
        3.622*1.5*0.075*(0.075*x+y)^0.5 + chi*x
}

f.simple <- function(x2y, x1, chi) {
    return(
        c(
            dg.simple(x2y[1], x2y[2], chi) - dg.simple(x1, x2y[2], chi),
            (x2y[1] - x1) * dg.simple(x1, x2y[2], chi) - 
                (g.simple(x2y[1], x2y[2], chi) - g.simple(x1, x2y[2], chi))
        )
    )
}

f.complex <- function(x2y, x1, chi) {
    Chi <- matrix(rep(0, 25), 5, 5)
    Chi[1,2] <- chi
    Chi[2,1] <- chi
    out <- binodal.curve.fun_(x2y, x1, Chi = Chi,
                       alpha = 3.622,
                       sigma = c(0.15, 0.15, 1, 1, 0),
                       polymer.num = c(1000, 1000, 1, 1, 1),
                       size.ratio = rep(1, 5),
                       molar.ratio = rep(1, 5)
                       )
}

ds.complex <- do.call(rbind, lapply(seq(-0.5, 0.5, 0.1), function(chi){
    
    ds <- do.call(rbind, lapply(seq(1e-7, 4e-2, 1e-4), function(x1) {
        
        x2y <- nleqslv(
                x = c(0.1, 0.1),
                fn = f.complex,
                global = 'pwldog',
                control = list(xtol = 1e-10),
                x1 = x1,
                chi = chi
            )$x
        
        data.frame(
            x1 = x1, 
            x2 = x2y[1],
            y = x2y[2],
            chi = chi
        ) %>% 
            filter(abs(x1 - x2) > 1e-4)
        
    }))
}))
ds.simple <- do.call(rbind, lapply(seq(-0.5, 0.5, 0.1), function(chi){
    
    ds <- do.call(rbind, lapply(seq(1e-7, 4e-2, 1e-4), function(x1) {
        
        x2y <- nleqslv(
                x = c(0.1, 0.1),
                fn = f.simple,
                global = 'pwldog',
                control = list(xtol = 1e-10),
                x1 = x1,
                chi = chi
            )$x
        
        data.frame(
            x1 = x1, 
            x2 = x2y[1],
            y = x2y[2],
            chi = chi
        ) %>% 
            filter(abs(x1 - x2) > 1e-4)
        
    }))
}))

ds <- rbind(
    ds.simple %>% mutate(type = 'simple'),
    ds.complex %>% mutate(type = 'complex')
    )

g <- ggplot(ds, aes(y = y, group = chi)) +
    geom_line(aes(x = x1, col = chi), lwd = 1.5) +
    geom_line(aes(x = x2, col = chi), lwd = 1.5) +
    scale_color_continuous(guide = 'legend', breaks = unique(ds$chi)) +
    labs(x = 'Polymer Frac.', y = 'Salt Frac.', title = 'Sigma = 0.15',
         col = 'Chi 1-2') +
    facet_wrap(~type)
g <- theme.title.text.2(g)

print(g)

ggsave('test.simple.chi.png', width = 5, height = 5)
