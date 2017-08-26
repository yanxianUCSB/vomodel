rm(list = ls())
source('vomodel.testing.R')
library(yxplot)
library(ggplot2)

g.simple.en <- function(x, y, chi) {
    x * 0.001 * log( 0.5 * x) +  y * log(0.5*y) + (1-x-y)*log(1-x-y) 
}
g.simple.el <- function(x, y, chi) {
    -3.622*(0.15*x + y)^1.5 
}
g.simple.chi <- function(x, y, chi) {
    (0.5*x) * (0.5*x) * chi * 2
}
g.simple <- function(x, y, chi) {
    # 2 * 0.5*x * 0.001 * log(0.5*x) + 2 * 0.5*y * log(0.5*y) + (1-x-y)*log(1-x-y) -
        # 3.622*(1*x + y)^1.5 + (0.5*x)^2 * chi
    g.simple.en(x, y, chi) + g.simple.el(x, y, chi) + g.simple.chi(x, y, chi)
}
dg.simple <- function(x, y, chi) {
    0.001 * log(0.5*x) + 0.001 -log(1-x-y) -1 -
        3.622*1.5*0.15*(0.15*x+y)^0.5 + 2*0.5*chi*(0.5*x) * 2
}

f.simple <- function(x2y, x1, chi) {
    return(
        c(
            dg.simple(x2y[1], x2y[2], chi) - dg.simple(x1, x2y[2], chi),
            (x2y[1] - x1) * dg.simple(x2y[1], x2y[2], chi) - 
                (g.simple(x2y[1], x2y[2], chi) - g.simple(x1, x2y[2], chi))
        )
    )
}


# chi <- -0.1
# Chi <- matrix(rep(0, 25), 5, 5)
# Chi[1,2] <- chi
# Chi[2,1] <- chi
# a <- gibbs.d(0.1, 0.1, Chi = Chi,
#            alpha = 3.622,
#            sigma = c(0.15, 0.15, 1, 1, 0),
#            polymer.num = c(1000, 1000, 1, 1, 1),
#            size.ratio = rep(1, 5),
#            molar.ratio = rep(1, 5))
# b <- dg.simple(0.1, 0.1, -0.1)
# print(c(a, b, a-b))
# stop()
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

ds.simple <- do.call(rbind, lapply(c(-0.2, 0, 0.2), function(chi){
    
    ds <- do.call(rbind, lapply(seq(1e-5, 0.01, 1e-4), function(x1) {
        
        sln <- nleqslv(
                x = c(0.05, 0.002),
                fn = f.simple,
                method = 'Newton',
                global = 'cline',
                # control = list(xtol = 1e-8),
                x1 = x1,
                chi = chi
            )
        
        x2y <- sln$x
        
        if (sln$termcd == 3) return()
        
        data.frame(
            x1 = x1, 
            x2 = x2y[1],
            y = x2y[2],
            chi = chi
        ) %>% 
            filter(abs(x1 - x2) > 1e-3)
        
    }))
}))
ds.complex <- do.call(rbind, lapply(c(-0.2, 0, 0.2), function(chi){
    
    ds <- do.call(rbind, lapply(seq(1e-5, 0.01, 1e-4), function(x1) {
        
        sln <- nleqslv(
                x = c(0.05, 0.002),
                fn = f.complex,
                method = 'Newton',
                global = 'cline',
                # control = list(xtol = 1e-8),
                x1 = x1,
                chi = chi
            )
        
        x2y <- sln$x
        
        if (sln$termcd == 3) return()
        
        data.frame(
            x1 = x1, 
            x2 = x2y[1],
            y = x2y[2],
            chi = chi
        ) %>% 
            filter(abs(x1 - x2) > 1e-3)
        
    }))
}))

ds <- rbind(
    ds.simple %>% mutate(type = 'simple'),
    ds.complex %>% mutate(type = 'complex')
    )

g <- ggplot(ds, aes(y = y, group = chi)) +
    geom_point(aes(x = x1, col = chi), lwd = 1.5) +
    geom_point(aes(x = x2, col = chi), lwd = 1.5) +
    scale_color_continuous(guide = 'legend', breaks = unique(ds$chi)) +
    labs(x = 'Polymer Frac.', y = 'Salt Frac.', title = 'Sigma = 0.15',
         col = 'Chi 1-2') +
    facet_wrap(~type)
g <- theme.title.text.2(g)

print(g)

ggsave('test.simple.chi.png', width = 5, height = 5)
