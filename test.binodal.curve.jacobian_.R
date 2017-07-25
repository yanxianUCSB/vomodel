
# test jacobian
rm(list = ls())
source('vomodel.R')
library(yxplot)
library(ggplot2)
library(rootSolve)
library(pracma)
library(minpack.lm)
source('test.peek.para.R')
DEBUG <<- T

print(binodal.curve.jacobian_(c(0.1, 0.1), 0.05, 
                              temp = temp,
                              alpha = alpha,
                              sigma = sigma,
                              Chi = 0,
                              polymer.num = polymer.num,
                              size.ratio = size.ratio,
                              epsilon = 1E-8,
                              guess.critical.point = c(phi.polymer=0.01, phi.salt=0.15),
                              binodal.guess = c(0.1, 0.1)))

# [,1]        [,2]
# [1,] -0.13853311 0.348242483
# [2,] -0.01150555 0.008953414