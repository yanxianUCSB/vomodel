k.vol                       <<- 1           # sample volume, doesn't affect the binodal curve
k.dna.contour.unit.length   <<- 0.33E-9     # 0.33 nm  (Blackburn 2006)
k.lysozyme.density          <<- 0.70E-6     # m^3/g (Perkins 1986)
k.amino.acid.length         <<- (k.lysozyme.density * 14307 / 129 / kNa)^(1/3)

system.properties <- list(
    polymer.num   = c(207, 900e3 / 306.2,   1, 1, 1),
    sigma         = c(11 / 207, 1,   1, 1, 0),
    size.ratio    = c(k.amino.acid.length, k.dna.contour.unit.length, k.na.size, k.cl.size, k.water.size) / k.water.size,
    MW            = c(22e3, 900e3),
    # charge.ratio  = c(1, 1),
    molar.ratio   = c(900e3 / 306.2,  11, 0.5, 0.5, 0),
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
    sampling.end = 5e-0,
    sampling.gap = 1e-5 ,
    critical.point.guess = c(phi.polymer = 0.001, phi.salt = 0.003) ,
    c.point.temp.fun = c.point.temp.fun(c.point.temp(system.properties, fitting.para)) ,
    binodal.guess = 0.1,  # phi.polymer.2
    condensation = T,
    counterion.release = T,
    default.critical.point = list(phi.polymer=0.1, phi.salt=0.01)
)

