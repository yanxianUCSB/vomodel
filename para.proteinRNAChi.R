test.chi.pw <- 0.0

system.properties <- list(
    polymer.num   = c(207, 900e3 / 306.2,   1, 1, 1),
    sigma         = c(11 / 207, 0.5,   1, 1, 0),
    size.ratio    = c(k.amino.acid.length, k.dna.contour.unit.length, k.na.size, k.cl.size, k.water.size) / k.water.size,
    MW            = c(22e3, 900e3),
    charge.ratio  = c(1, 1),
    molar.ratio   = c(900e3 / 306.2,  11, 0.5, 0.5, 0),
    # the polycation:polyanion and cation:anion molar ratio
    water.size    = k.water.size,
    Chi = matrix(c(.0,.0,0,0,test.chi.pw,
                   .0,.0,0,0,0,
                   0,0,0,0,0,
                   0,0,0,0,0,
                   test.chi.pw,0,0,0,0), 5,5),
    lattice.spacing = k.water.size
)
fitting.para <- list(
    epsilon = 1E-8 , 
    sampling.start = 1e-5,
    sampling.gap = 1e-7 ,
    critical.point.guess = c(phi.polymer = 0.01, phi.salt = 0.0001) ,
    c.point.temp.fun = c.point.temp.fun(c.point.temp(system.properties, fitting.para)) ,
    binodal.guess = 0.05  # phi.polymer.2
)