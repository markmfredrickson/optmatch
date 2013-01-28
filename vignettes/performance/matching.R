library(profr)
library(optmatch)

load("setup.rda")
load("distance.rda")

gc()
benchmark.matching.dense <- profr(result.matching.dense <-
  fullmatch(result.dense))


gc()
benchmark.matching.sparse <- profr(result.matching.sparse <-
  fullmatch(result.sparse.caliper))


gc()
benchmark.matching.sparse.strat <- profr(result.matching.sparse.strat <-
  fullmatch(result.sparse.within))


gc()
benchmark.matching.pairmatch <- profr(result.matching.pairmatch <-
  fullmatch(result.dense))


save(file = "matching.rda",
  benchmark.matching.dense,
  result.matching.dense,
  benchmark.matching.sparse,
  result.matching.sparse,
  benchmark.matching.sparse.strat,
  result.matching.sparse.strat,
  benchmark.matching.pairmatch,
  result.matching.pairmatch)


