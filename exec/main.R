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
