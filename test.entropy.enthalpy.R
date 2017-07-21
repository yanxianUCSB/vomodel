library(ggplot2)
library(reshape2)
source('vomodel.R')
source('para.spruijt.R')

test.entropy.enthalpy <- function(totl.conc, nacl.conc, Para){

# Entropy vs Conc
entropy <- sapply(totl.conc, function(t.conc){
  conc <- get.conc(t.conc, nacl.conc, Para)
  phis <- get.phi(conc, Para)
  
  return(Fen(phis, Para$polym.num, Para$size.ratio))
})


# Enthalpy vs Conc
enthalpy <- sapply(totl.conc, function(t.conc){
  conc <- get.conc(t.conc, nacl.conc, Para)
  phis <- get.phi(conc, Para)
  return(Fel(alpha = alpha(temp = temp, size = Para$size), sigma = Para$charge.den, phi = phis))
})

# Entropy + Enthalpy
ds <- data.frame(conc = totl.conc, entropy = entropy, enthalpy = enthalpy, energy = entropy + enthalpy)

return(ds)
}


totl.conc <- seq(0, 2500, 5)
nacl.conc <- 0.0001
temp <- 300
ds <- test.entropy.enthalpy(totl.conc, nacl.conc, Para)
ds <- melt(data = ds, id.vars = c("conc"))
g <- ggplot(ds, aes(x = conc, y = value, col = variable)) +
  geom_line()
print(g)
# ggsave('test.entropy.enthalpy.png', width = 5, height = 5)
