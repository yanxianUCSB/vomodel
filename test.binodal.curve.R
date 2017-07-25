# test binodal curve
rm(list = ls())
source('vomodel.R')
library(yxplot)
library(ggplot2)
library(rootSolve)
library(pracma)
library(minpack.lm)
source('test.peek.para.R')
DEBUG <<- T
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

p <- binodal.curve_(
                   temp = temp,
                   alpha = alpha,
                   sigma = sigma,
                   Chi = 0,
                   polymer.num = polymer.num,
                   size.ratio = size.ratio,
                   epsilon = 1E-8,
                   guess.critical.point = c(phi.polymer=0.01, phi.salt=0.15),
                   binodal.guess = c(0.1, 0.1)
                   )

# 
p <- as.data.frame.matrix(p)
g <- ggplot(p, aes(y = phi.salt)) +
  geom_line(aes(x = phi.polymer.1), col = 'red') +
  geom_line(aes(x = phi.polymer.2), col = 'blue') +
  labs(x = 'Polymer [phi]',
       y = 'Salt [phi]')
g <- theme.title.text.1(g)
print(g)
ggsave('test.binodal.curve.png', width = 5, height = 5)
