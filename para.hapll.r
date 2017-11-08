k.vol                       <<- 1           # sample volume, doesn't affect the binodal curve

system.properties <- list(
    polymer.num   = c(1000, 1000,   1, 1, 1),
    sigma         = c(1, 1, 1, 1, 0),
    size.ratio    = c(1, 1, 1, 1, 1),     # Estimated from 0.50 nm, 1.0 nm, ion size, source: wikipedia
    MW            = c(146-18, 776-18),   # MW of lysine and monomer hyaluronic acid
    molar.ratio   = c(1, 1, 1, 1, 0),    # the polycation:polyanion and cation:anion molar ratio
    water.size    = k.water.size,
    Chi = matrix(rep(0, 25), 5, 5)
)
fitting.para <- list(
    epsilon = 1E-8 , 
    sampling.start = 1e-10,
    sampling.end = 0.05,
    sampling.gap = 5e-6 ,
    critical.point.guess = c(phi.polymer = 0.001, phi.salt = 0.003) ,
    binodal.guess = c(1e-3, 1e-3),  
    condensation  = F,
    counterion.release = T,
    default.critical.point = list(phi.polymer=0.1, phi.salt=0.1)
)

