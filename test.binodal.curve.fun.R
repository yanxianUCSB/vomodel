# test binodal curve function
rm(list = ls())
source('vomodel.R')
library(yxplot)
library(ggplot2)
library(rootSolve)
library(pracma)
source('test.peek.para.R')
DEBUG <<- T

p <- lapply(phi.polymer.seq[which(phi.polymer.seq < 0.003)], function(phi.polymer1){
  q <- lapply(phi.polymer.seq[which(phi.polymer.seq > 0.003)], function(phi.polymer2){
  bnf.out <- binodal.curve.fun(c(phi.polymer1, phi.polymer2),
                       phi.salt = phi.salt,
                   temp = temp,
                   alpha = alpha,
                   sigma = sigma,
                   Chi = 0,
                   polymer.num = polymer.num,
                   size.ratio = size.ratio,
                   guess.critical.point = c(phi.polymer=0.01, phi.salt=0.15)
                   )
    out <- list(
      phi.polymer1 = phi.polymer1,
      phi.polymer2 = phi.polymer2,
      f1 = bnf.out[1],
      f2 = bnf.out[2]
                )
    return(out)
  })
  return(do.call(rbind, q))
})

p <- (do.call(rbind, p))
p <- as.data.frame.matrix(as.matrix(p))
p <- as.data.frame.matrix(sapply(p, as.numeric))
head(p)
str(p)
g <- ggplot(p, aes(x = phi.polymer1, y = phi.polymer2, fill = (f2))) +
  geom_tile() +
  scale_fill_gradient2(low = 'blue', high = 'red')
print(g)