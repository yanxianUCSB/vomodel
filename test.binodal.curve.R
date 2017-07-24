# test binodal curve
rm(list = ls())
source('vomodel.R')
library(yxplot)
library(ggplot2)
library(rootSolve)
library(pracma)
source('test.peek.para.R')
DEBUG <<- F
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
# p <- sapply(seq(0.0001, 0.017, 1e-3), function(phi.salt){
#   out <- binodal.curve.fun(c(0.01,0.1),
#                        phi.salt = phi.salt,
#                    temp = temp,
#                    alpha = alpha,
#                    sigma = sigma,
#                    Chi = 0,
#                    polymer.num = polymer.num,
#                    size.ratio = size.ratio,
#                    guess.critical.point = c(phi.polymer=0.01, phi.salt=0.15)
#                    )
#   return(c(phi.salt = phi.salt, f1 = out[1], f2 = out[2]))
# })

p <- binodal.curve(phi.polymer.seq,
                   temp = temp,
                   alpha = alpha,
                   sigma = sigma,
                   Chi = 0,
                   polymer.num = polymer.num,
                   size.ratio = size.ratio,
                   guess.critical.point = c(phi.polymer=0.01, phi.salt=0.15),
                   epsilon = 1E-8
                   )
print(p)
head(p)
str(p)

p <- as.data.frame.matrix(p)
g <- ggplot(p, aes(y = phi.salt)) +
  geom_point(aes(x = phi.polymer1), col = 'red') +
  geom_point(aes(x = phi.polymer2), col = 'blue')
print(g)