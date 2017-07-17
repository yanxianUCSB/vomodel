# test
source('./github/vomodel/vomodel.R')

# Parameters
# c(+, -, w, sp, sn)
Para <- list()
Para$size <- k.water.size 
Para$charge.den <- c(11/207, 1, 0, 1, 1)
Para$polym.num <- c(207, 900E3 / (324 - 18), 1, 1, 1)  # polym.num. REF[Yanxian's Notebook]
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
temps <- c(700, 900, 1100, 1500, 3000)
tot.concs <- seq(100, 16000, 100)
salt.concs <- c(1, 10000)
ds <- get.free.energy_(temps = temps, tot.concs = tot.concs, salt.concs = salt.concs, list(Chi), list(Para)) %>%
    group_by(sConc, temp) %>%
    mutate(free.energy.change = free.energy - free.energy[1]) %>%
    ungroup()


library(ggplot2)
g <- ggplot(ds, aes(x = log(tConc), y = free.energy.change, group = temp))
# g <- g + geom_point(aes(col = sConc, pch = as.factor(temp)))
g <- g + geom_line(aes(col = (temp)), lwd = 1.25)
g <- g + facet_wrap(~sConc)
# g <- g + facet_wrap(~temp)
g <- g + scale_color_gradient(low = "#56B1F7", high = '#132B43')
g <- theme.title.text.2(g)
ggsave('./github/vomodel/free.energy.change.png', height = 4, width = 8)

g <- ggplot(ds, aes(x = log(tConc), y = entropy, group = temp))
# g <- g + geom_point(aes(col = sConc, pch = as.factor(temp)))
g <- g + geom_line(aes(col = (temp)), lwd = 1.25)
g <- g + facet_wrap(~sConc)
g <- g + scale_color_gradient(low = "#56B1F7", high = '#132B43')
g <- theme.title.text.2(g)
ggsave('./github/vomodel/entropy.change.png', height = 4, width = 8)

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