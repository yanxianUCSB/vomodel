# test spruijt
library(dplyr)
library(yxplot)
source('./vomodel.R')
source('para.spruijt.R')
# free energy landscape at temperature and concentration
totl.conc <- seq(0.1, 2400, 1)
# nacl.conc <- seq(0, 1000, 100)
nacl.conc <- 100
temps <- c(100)
tot.concs <- totl.conc
salt.concs <- nacl.conc
ds <-
  get.free.energy_(
    temps = temps,
    tot.concs = tot.concs,
    salt.concs = salt.concs,
    list(Chi),
    list(Para)
  ) %>%
  group_by(sConc, temp) %>%
  mutate(free.energy.change = free.energy - free.energy[1]) %>%
  mutate(entropy.change = entropy - entropy[1],
         enthalpy.el.change = enthalpy.el - enthalpy.el[1]) %>%
  mutate(free.energy.dif = c(diff(free.energy)/diff(tot.concs),  NA),
         free.energy.difdif = c(diff(free.energy.dif) / diff(tot.concs)[-1], NA)) %>%
  ungroup()

library(ggplot2)
peek <- function(peek.col) {
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
  g <-
    ggplot(ds, aes(
      x = tConc,
      y = temp,
      fill = (eval(parse(text = peek.col)))
    )) +
    geom_tile() +
    scale_fill_gradientn(colors = jet.colors(7), guide = guide_colorbar(title = peek.col)) +
    facet_wrap(~ sConc)
  g <- theme.title.text.1(g)
  return(g)

}
peek.curve <- function(peek.col) {
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
  g <-
    ggplot(ds, aes(
      x = tConc,
      y = (eval(parse(text = peek.col))),
      group = temp
    )) +
    geom_line(aes(col = temp)) +
    labs(y = '') +
    scale_color_gradientn(colors = jet.colors(7), guide = guide_legend(title = peek.col)) +
    facet_wrap(~ sConc)
  g <- theme.title.text.1(g)
  return(g)

}

g <- list()

g[[length(g)+1]] <- peek.curve('entropy')
g[[length(g)+1]] <- peek.curve('enthalpy.el')
g[[length(g)+1]] <- peek.curve('free.energy.change')
g[[length(g)+1]] <- peek.curve('free.energy.dif')
g[[length(g)+1]] <- peek.curve('free.energy.difdif')
g[[length(g)+1]] <- peek('entropy')
g[[length(g)+1]] <- peek('enthalpy.el')
g[[length(g)+1]] <- peek('free.energy.change')
g[[length(g)+1]] <- peek('free.energy.dif')
g[[length(g)+1]] <- peek('free.energy.difdif')



# multiplot(g[[1]], g[[2]], g[[3]], cols = 1)
multiplot(g[[4]], g[[5]])
# ggsave(filename = 'test.spruijt.curve.png', width = 9, height = 5)
# ggsave(filename = 'test.spruijt.tile.png', width = 8, height = 8)

