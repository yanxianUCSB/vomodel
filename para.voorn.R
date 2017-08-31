system.properties <- list(
    polymer.num   = c(1000, 1000,   1, 1, 1),  # Degree of polymerization
    sigma         = c(0.15, 0.15,   1, 1, 0),  # charge density from 0-1
    size.ratio    = rep(1, 5),  # length of the cubic molecule relative to water
    MW            = c(1, 1),  # molecular weight, only used in transforming molar concentration to mass conc
    molar.ratio   = c(1, 1, 1, 1, 0),  # relative molar ratio to match charge neutrality
    water.size    = k.water.size,  # constant, 0.31 nm
    Chi = matrix(rep(0, 25), 5, 5)  # Chi matrix
)
fitting.para <- list(
    epsilon = 1E-8,  # only used in numerical solving function, xtol = epsilon
    sampling.start = 1e-5,  # starting valume fraction of binodal curve 
    sampling.gap = 1e-5 ,  # sampling interval
    sampling.end = 1e-2,  # end point, usually 0.01 is enough
    binodal.guess = c(1e-4, 1e-4),  # initial binodal guess, the program will attempt to guess on multiple initial guess
    condensation = F,  # control for condensation mode
    counterion.release = F,  # control for counterion release mode
    default.critical.point = list(phi.polymer=0.05, phi.salt=0.02)  # the end point for fitting
)