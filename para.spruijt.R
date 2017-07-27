# Parameters
# c(+, -, w, sp, sn)
Para <- list()
Para$MW <- c(1, 1, 0, 0, 0)  # PDAEMA, PAA, Na, Cl, water, g/mol
Para$mass.ratio <- c(0.5, 0.5, 0, 0, 0)
Para$water.size <- k.water.size  # meter
Para$charge.den <-
  c(1, 1, 1, 1, 0)  # protein, rna, na, cl, h2o
Para$polym.num <-
  c(1, 1, 1, 1, 1)  # polym.num. REF[Yanxian's Notebook]
Para$size.ratio <- c(1, 1, 1, 1, 1)
chi.ab <- 0
Chi <- matrix(c(0, chi.ab, 0, 0, 0, chi.ab, 0, 0, 0, 0, rep(0, 10)), 5, 5)  # test Chi = 0 first
