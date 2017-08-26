rm(list = ls())
# options(digits = 22)
source('vomodel.R')
source('para.voorn.chi.R')
DEBUG <- TRUE
# TEST 
phi.polymer <- seq(0.001, 0.4, 0.001)
phi.salt <- 0.000155
temp <- 300
chipw <- -0.9
Chi <- matrix(c(.0,.0,0,0,chipw,
                .0,.0,0,0,0,
                0,0,0,0,0,
                0,0,0,0,0,
                chipw,0,0,0,0), 5,5)


sysprop <- system.properties

# TEST gibbs and gibbs dirivative
ds <- gibbs.funs(phi.polymer, 
                 phi.salt = phi.salt,
                 temp = temp,
                 Chi = Chi,
                 alpha = get.alpha(temp, k.water.size),
                 sigma = sysprop$sigma,
                 polymer.num = sysprop$polymer.num,
                 size.ratio = sysprop$size.ratio,
                 molar.ratio = sysprop$molar.ratio)
print(ds)
plot(ds[-2])
head(ds)