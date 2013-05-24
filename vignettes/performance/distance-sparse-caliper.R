library(profr)
library(optmatch)

load("setup.rda")

gc()
benchmark.sparse.caliper <- profr(result.sparse.caliper <- match_on(x =
predicted, z = DATA$Z, caliper = 1))

save(file = "distance-sparse-caliper.rda",
     benchmark.sparse.caliper,
     result.sparse.caliper)
