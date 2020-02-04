# Data preparation ---------------
# this part of code transform the experimental observations into volume fraction
# and temperature in Kelvin
get_expt  <- function(path_experiment = commandArgs(trailingOnly = T)[1]) {
    library(dplyr)
    # path_experiment <- '~/Desktop/dataset.csv'
    ds <- read.csv(path_experiment) %>%
        mutate(
            phi1 = protein * 1E-6 * 207 / k.water.conc * 1000,  # k.water.conc in mol/m^3
            phi3 = 0.5 * (nacl + 20) * 1E-3 / k.water.conc * 1000,  # 20 mM monovalent buffer salt
            temp = cloudpoint + 273.15
        ) %>%
        select(phi1, phi3, temp)
    ds %>%
        saveRDS('expt.rds')
    return(ds)
}


# Launch get.chi.from.experiment.R-----
source('R/get.chi.from.experiment.R')
simulate.fh()
simulate.fhvo()
simulate.fhvocr()

# Launch get.chitemp.function.r-----
source('R/get.chitemp.function.r')
FH <- get.chitemp.function('FH')
FHVO <- get.chitemp.function('FHVO')
FHVO.1salt <- get.chitemp.function.1salt('FHVO')
FHVOCR <- get.chitemp.function('FHVOCR')
saveRDS(list(
    FH=FH,
    FHVO=FHVO,
    FHVO.1salt=FHVO.1salt,
    FHVOCR=FHVOCR
), 'results/chi_temp_functions.rds')


source('R/get.phase.diagram.using.chitempfunction.r')
simulate()
simulate.1.eg()

source('R/get.phase.diagram.using.chitempfunction.temp.phi3.r')
simulate()

# Launch get.dataset.using.phase.diagram.for.plotting.R-----
source('R/get.dataset.using.phase.diagram.for.plotting.R')
save.dataset.FHVO('results/')
save.dataset.FHVO('results_temp_phi3/')
save.dataset('results/')
save.dataset('results_temp_phi3/')


source('exec/plot.fitting.results.R')
