# TEST parameters
test.ABSkpq <- function(sysprop) {
    ds <- pp(sysprop = sysprop)
    print(c('kpq = ', ds$kpq))
    print(c('S = ', ds$S))
    print(c('A = ', ds$A))
    print(c('B = ', ds$B))
}