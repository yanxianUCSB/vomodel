# test binodal curve
rm(list = ls())
source('vomodel.R')
library(yxplot)
library(ggplot2)
library(rootSolve)
library(pracma)
source('test.peek.para.R')
DEBUG <<- T
phi.salt <- 0.1
a <- mean(sapply(phi.polymer.seq, gibbs.d, phi.salt, 
                 temp = temp,
                 alpha = alpha,
                 sigma = sigma,
                 Chi = 0,
                 polymer.num = polymer.num,
                 size.ratio = size.ratio
                 ))
b <- mean(sapply(phi.polymer.seq, gibbs, phi.salt, temp = temp,
                 alpha = alpha,
                 sigma = sigma,
                 Chi = 0,
                 polymer.num = polymer.num,
                 size.ratio = size.ratio)[1:10])

# binodal.curve.fun(c(a, b), 
#                   phi.salt = phi.salt,
#                   temp = temp,
#                   alpha = alpha,
#                   sigma = sigma,
#                   Chi = 0,
#                   polymer.num = polymer.num,
#                   size.ratio = size.ratio)

# 
# p <- binodal.curve.fun(c(0.01,0.1),
#                        phi.salt = 0.15,
#                    temp = temp,
#                    alpha = alpha,
#                    sigma = sigma,
#                    Chi = 0,
#                    polymer.num = polymer.num,
#                    size.ratio = size.ratio,
#                    guess.critical.point = c(phi.polymer=0.01, phi.salt=0.15)
#                    )
p <- binodal.curve(phi.polymer.seq,
                   temp = temp,
                   alpha = alpha,
                   sigma = sigma,
                   Chi = 0,
                   polymer.num = polymer.num,
                   size.ratio = size.ratio,
                   guess.critical.point = c(phi.polymer=0.01, phi.salt=0.15)
                   )
print(p)
par(mfrow = (c(1,2)))
plot(p[,1])
plot(p[,2])

plot(y$phi, y$f)
abline(a = p$par[1], b = p$par[2], col = 'red')
