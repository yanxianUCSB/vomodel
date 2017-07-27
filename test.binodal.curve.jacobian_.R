
# test jacobian
rm(list = ls())
source('vomodel.R')
library(yxplot)
library(ggplot2)
library(rootSolve)
library(pracma)

phi.polymer.seq <- c(seq(1E-5, 0.14, 0.001))
binodal.guess <- c(0.1, 0.1)  # phi.polymer.2, phi.salt
phi.salt <- 0.150
temp <- 300
polymer.num <- c(1000, 1000, 1, 1, 1)
alpha <- 3.655
sigma <- c(0.44, 0.44, 1, 1, 0)
size.ratio <- c(1, 1, 1, 1, 1)
#  (k.vol * kkB * arg$temp) / k.water.size ^ 3  = 1
k.vol = k.water.size^3 / kkB / temp

DEBUG <<- T

print(binodal.curve.jacobian_(c(0.1, 0.1), 0.05, 
                              temp = temp,
                              alpha = alpha,
                              sigma = sigma,
                              Chi = 0,
                              polymer.num = polymer.num,
                              size.ratio = size.ratio,
                              molar.ratio = c(1, 1, 1, 1, 0),
                              epsilon = 1E-8,
                              guess.critical.point = c(phi.polymer=0.01, phi.salt=0.15),
                              binodal.guess = c(0.1, 0.1)))
cat('if your gibbs and gibbs derivative functions does not change then you will get\n')
cat(
"            [,1]        [,2]
[1,] -0.13853311 0.348242483
[2,] -0.01150555 0.008953414")

