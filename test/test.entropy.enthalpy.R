rm(list = ls())
library(ggplot2)
library(reshape2)
source('vomodel.R')
# Parameters
# c(+, -, w, sp, sn)
Para <- list()
Para$MW <- c(1, 1, 0, 0, 0)  # PDAEMA, PAA, Na, Cl, water, g/mol
Para$mass.ratio <- c(0.5, 0.5, 0, 0, 0)
Para$water.size <- k.water.size  # meter
Para$charge.den <-
    c(1, 1, 1, 1, 0)  # protein, rna, na, cl, h2o
Para$polym.num <-
    c(1, 1, 1, 1, 1)  # polym.num. REF[Yanxian's Notebook]
Para$size.ratio <- c(1, 1, 1, 1, 1)
chi.ab <- 0
Chi <- matrix(c(0, chi.ab, 0, 0, 0, chi.ab, 0, 0, 0, 0, rep(0, 10)), 5, 5)  # test Chi = 0 first


test.entropy.enthalpy <- function(totl.conc, nacl.conc, Para) {
    # Entropy vs Conc
    entropy <- sapply(totl.conc, function(t.conc) {
        conc <- get.conc(t.conc, nacl.conc, Para)
        phis <- get.phi(conc, Para)
        
        return(Fen(phis, Para$polym.num, Para$size.ratio))
    })
    
    
    # Enthalpy vs Conc
    enthalpy <- sapply(totl.conc, function(t.conc) {
        conc <- get.conc(t.conc, nacl.conc, Para)
        phis <- get.phi(conc, Para)
        return(Fel(
            alpha = get.alpha(temp = temp, size = Para$water.size),
            sigma = Para$charge.den,
            phi = phis
        ))
    })
    
    # Entropy + Enthalpy
    ds <-
        data.frame(
            conc = totl.conc,
            entropy = entropy,
            enthalpy = enthalpy,
            energy = entropy + enthalpy
        )
    
    return(ds)
}


totl.conc <- seq(0, 50, 1e-3)
nacl.conc <- 0.0001
temp <- 300
ds <- test.entropy.enthalpy(totl.conc, nacl.conc, Para)
ds <- melt(data = ds, id.vars = c("conc"))
g <- ggplot(ds, aes(x = conc, y = value, col = variable)) +
    geom_line() +
    labs(x = 'Mol/L for polymer at nacl 100 mM, temp 300 K')
print(g)
ggsave('test.entropy.enthalpy.png', width = 5, height = 5)
