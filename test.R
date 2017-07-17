# test
library(dplyr)
library(yxplot)
source('./vomodel.R')

normed <- function(x) {
    return((x - x[1]))
}
# Parameters
# c(+, -, w, sp, sn)
Para <- list()
Para$size <- k.water.size
Para$charge.den <- c(11 / 207, 1, 1, 1, 0)
Para$polym.num <-
    c(207, 900E3 / (324 - 18), 1, 1, 1)  # polym.num. REF[Yanxian's Notebook]
Para$size.ratio <- c(1, 1, 1, 1, 1)
Chi <- matrix(rep(0, 25), 5, 5)  # test Chi = 0 first

# # Bjerrum length at 300K ~0.7nm
# temp <- 300
# lB <- ke ^ 2 / (kEr * kkB * temp)
# cat('Bjerrum length', temp, 'K:', lB)

# # # Conc vs Phis
# totl.conc <- seq(10, 5000, 1)
# nacl.conc <- 50
# phis <- do.call(rbind, lapply(totl.conc, function(t.conc){
#     conc <- get.conc(t.conc, nacl.conc)
#     phis <- get.phi(conc, Para)
#     return(phis)
# }))
# plot(totl.conc, phis[, 1])

# # Entropy vs Conc
# entropy <- sapply(totl.conc, function(t.conc){
#     conc <- get.conc(t.conc, nacl.conc)
#     phis <- get.phi(conc, Para)
#
#     return(Fen(phis, Para$polym.num, Para$size.ratio))
# })
# plot(totl.conc, entropy)

# # Enthalpy vs Conc
# enthapy <- sapply(totl.conc, function(t.conc){
#     conc <- get.conc(t.conc, nacl.conc)
#     phis <- get.phi(conc, Para)
#     return(Fel(alpha = alpha(temp = temp, size = Para$size), sigma = Para$charge.den, phi = phis))
# })
# plot(totl.conc, enthapy)

# # Entropy + Enthalpy
# plot(totl.conc, enthapy + entropy)

# # free energy
# free.energy <- sapply(totl.conc, function(t.conc){
#     return(get.free.energy2(tempC = temp - 273, tot.conc = t.conc, salt.conc = nacl.conc,
#                             Chi = Chi, Para = Para))
# })
# plot(totl.conc, free.energy)

# free energy landscape at temperature and concentration
temps <- seq(200, 400, 1)
tot.concs <- seq(100, 6000, 100)
salt.concs <- c(1, 10000)
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
    ungroup()

library(ggplot2)
# ## Curve
# g <- ggplot(ds, aes(x = (tConc), y = free.energy, group = temp))
# # g <- g + geom_point(aes(col = sConc, pch = as.factor(temp)))
# g <- g + geom_line(aes(col = (temp)), lwd = 1.25)
# g <- g + facet_wrap(~sConc)
# # g <- g + facet_wrap(~temp)
# g <- g + scale_color_gradient(low = "#56B1F7", high = '#132B43')
# g <- theme.title.text.2(g)
# ggsave('./free.energy.change.png', height = 4, width = 8)
#
# g <- ggplot(ds, aes(x = log(tConc), y = entropy, group = temp))
# # g <- g + geom_point(aes(col = sConc, pch = as.factor(temp)))
# g <- g + geom_line(aes(col = (temp)), lwd = 1.25)
# g <- g + facet_wrap(~sConc)
# g <- g + scale_color_gradient(low = "#56B1F7", high = '#132B43')
# g <- theme.title.text.2(g)
# ggsave('./entropy.change.png', height = 4, width = 8)
#
# g <- ggplot(ds, aes(x = log(tConc), y = enthalpy.el, group = temp))
# # g <- g + geom_point(aes(col = sConc, pch = as.factor(temp)))
# g <- g + geom_line(aes(col = (temp)), lwd = 1.25)
# g <- g + facet_wrap(~sConc)
# g <- g + scale_color_gradient(low = "#56B1F7", high = '#132B43')
# g <- theme.title.text.2(g)
# ggsave('./enthalpy.el.change.png', height = 4, width = 8)

## Tile
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
        facet_wrap( ~ sConc) 
    g <- theme.title.text.1(g)
    ggsave(paste0(peek.col, '.tile.png'),
           height = 4,
           width = 8)
}

peek('entropy')
peek('enthalpy.el')
peek('free.energy')
peek('free.energy.change')



#
# total.conc <- 0.25
# nacl.conc <- 50
# conc <- get.conc(total.conc, nacl.conc)  # 10uM protein, 30ug/mL RNA, 55555mol/m^3 water, 50mM monovalent salt
# # volume fraction
# phis <- phi(conc, Para$size.ratio, Para$size, Para$polym.num)
# sum(phi(conc, Para$size.ratio, Para$size, Para$polym.num))
# # free.energy sign
# get.free.energy(temp, conc, Chi, Para)
#
# # LCST
# temps <- seq(100, 500, 1)
# energy <- get.free.energy(temps, conc, Chi, Para)
# plot(temps, energy)
#
#
# Fens <- Fen(phis, Para$polym.num, Para$size.ratio)