system.properties <- list(
    polymer.num   = c(1000, 1000,   1, 1, 1),
    sigma         = c(1, 1,   1, 1, 0),
    size.ratio    = rep(1, 5),
    MW            = c(1, 1),
    charge.ratio  = c(1, 1),
    molar.ratio   = c(1, 1, 1, 1, 0),
    # the polycation:polyanion and cation:anion molar ratio
    water.size    = k.water.size,
    Chi = matrix(c(.0,.0,0,0,0,
                   .0,.0,0,0,0,
                   0,0,0,0,0,
                   0,0,0,0,0,
                   0,0,0,0,0), 5,5)
)
fitting.para <- list(
    epsilon = 1E-8 , 
    sampling.start = 1e-5,
    sampling.gap = 1e-5 ,
    sampling.end = 5e-0,
    critical.point.guess = c(phi.polymer = 0.01, phi.salt = 0.001) ,
    c.point.temp.fun = c.point.temp.fun(c.point.temp(system.properties, fitting.para)) ,
    binodal.guess = 0.1,  # phi.polymer.2
    condensation = F,
    counterion.release = F
)