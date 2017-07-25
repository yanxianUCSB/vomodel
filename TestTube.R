library(R6)
TestTube <- R6Class(
  'TestTube',
  private = list(
    ..chemicals = NULL,
    # names
    ..concs = c(0),
    # mol / L
    ..temp = 300,
    # K
    ..pH = 7.0,
    ..ionic.strength = 0,
    ..volume = 1  # m^3
  ),
  public = list(
    initialize = function() {
      
    },
    add.chemical = function(chemical, conc) {
      assert_that(is.string(chemical))
      if (chemical %in% private$..chemicals) {
        private$..concs[which(private$..chemicals == chemical)] <-
          private$..concs[which(private$..chemicals == chemical)] + conc
      } else {
        private$..chemicals <- append(private$..chemicals, chemical)
        assert_that(is.number(conc), conc >= 0)
        private$..concs <- append(private$..cons, conc)
      }
    },
    remove.chemical = function(chemical) {
      assert_that(
        is.string(chemical),
        chemical %in% private$..chemicals,
        msg = paste(chemical, 'is not in the TestTube')
      )
      index <- which(chemical == private$..chemicals)
      private$..chemicals[-index] -> private$..chemicals
      private$..concs[-index] -> private$..concs
    },
    add.solvent = function(vol) {
      assert_that(is.numeric(vol))
    }
  ),
  active = list(
    chemicals = function() {
      list(chemicals =
             private$..chemicals,
           concs = private$..concs)
    },
    concs = function() {
      private$..concs
    }
  )
)

# test
testtube <- TestTube$new()
testtube$add.chemical('NaCl', 1)
testtube$remove.chemical('NaCl')
testtube$chemicals
