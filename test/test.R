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
