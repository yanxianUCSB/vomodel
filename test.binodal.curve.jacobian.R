# test.binodal.curve.jacobian.R
# test binodal curve function
rm(list = ls())
source('vomodel.R')
library(yxplot)
library(ggplot2)
library(rootSolve)
library(pracma)
source('test.peek.para.R')
DEBUG <<- T
phi.salt <- 0.13
guess <- c(0.001, 0.02)

j1 <- binodal.curve.jacobian(guess, phi.salt = phi.salt,
                            temp = temp,
                            alpha = alpha,
                            sigma = sigma,
                            Chi = 0,
                            polymer.num = polymer.num,
                            size.ratio = size.ratio,
                            guess.critical.point = c(phi.polymer=0.01, phi.salt=0.15))
j1inv <- inv(j1)

j2 <- jacobian(binodal.curve.fun, guess,                        
               phi.salt = phi.salt,
               temp = temp,
               alpha = alpha,
               sigma = sigma,
               Chi = 0,
               polymer.num = polymer.num,
               size.ratio = size.ratio,
               guess.critical.point = c(phi.polymer=0.01, phi.salt=0.15))
j2inv <- inv(j2)

print(c('gibbs.d', gibbs.d(guess[1], guess[2], temp = temp,
                         alpha = alpha,
                         sigma = sigma,
                         Chi = 0,
                         polymer.num = polymer.num,
                         size.ratio = size.ratio,
                         guess.critical.point = c(phi.polymer=0.01, phi.salt=0.15))))
print(c('gibbs.d.num', grad(gibbs, guess[1], 
                                phi.salt = guess[2], 
                                temp = temp,
                         alpha = alpha,
                         sigma = sigma,
                         Chi = 0,
                         polymer.num = polymer.num,
                         size.ratio = size.ratio,
                         guess.critical.point = c(phi.polymer=0.01, phi.salt=0.15))))

print(j1)
print(j2)
print(j1inv)
print(j2inv)

print(c('det j1', det(j1)))
print(c('det j1inv', det(j1inv)))
print(c('det j2', det(j2)))
print(c('det j2inv', det(j2inv)))