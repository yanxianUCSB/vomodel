# test.freeEnergy.grad free energy
# using pracama::grad() to generate derivative

rm(list = ls())
# options(digits = 22)
source('vomodel.R')
source('para.voorn.chi.R')
DEBUG <- TRUE

# TEST CONSTANTS
source('test.Bjerrumlength.R')
# TEST Parameters
source('test.ABSkpq.R')
test.ABSkpq(sysprop = system.properties)
# TEST jacobian
source('test.binodal.curve.jacobian_.R')
test.binodal.curve.jacobian_()




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
stop()
ds2 <- free.energy.funs(
    phi.polymer[1], 
    phi.salt = phi.salt,
    sysprop = system.properties,
    fitting.para = fitting.para,
    temp = temp,
    Chi = Chi,
    alpha = get.alpha(temp, k.water.size),
    sigma = sysprop$sigma,
    polymer.num = sysprop$polymer.num,
    size.ratio = sysprop$size.ratio,
    molar.ratio = sysprop$molar.ratio
)
plot(ds2)
par(mfrow = c(2,2))
plot(ds2$f - ds$f)
plot(ds2$df - ds$df)
plot(ds2$ddf - ds$ddf)
plot(ds2$dddf - ds$dddf)
