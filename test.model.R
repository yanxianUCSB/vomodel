# test.model
rm(list = ls())
library(dplyr)

jet.colors <-
    colorRampPalette(
        c(
            "#00007F",
            "blue",
            "#007FFF",
            "cyan",
            "#7FFF7F",
            "yellow",
            "#FF7F00",
            "red",
            "#7F0000"
        )
    )

x  <- seq(0, 1, 0.01)
k <- seq(-10, -5, 0.1)
# k <- -5
ds <- expand.grid(x = x, k = k)
ds$dy <- log(ds$x / (1 - ds$x)) + ds$k * sqrt(ds$x) + 5
ds$y <- with(ds, x * log(x) + (1-x) * log(1-x) + k * sqrt(x) + 5*x )
# ds$dy.num <- numericDeriv(y ~ x * log(x) + (1-x) * log(1-x) + k * sqrt(x) + 5*x, theta = c('x'),)

library(ggplot2)

g <- ggplot(ds, aes(x = x, y = k, fill = y)) + geom_tile() +
    scale_fill_gradientn(colors = jet.colors(7)) 
    

# plot(ds)
