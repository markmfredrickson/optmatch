library(profr)
library(optmatch)

load("setup.rda")

gc()
benchmark.dense <- profr(result.dense <- match_on(x = predicted, z = DATA$Z))

gc()
benchmark.sparse.caliper <- profr(result.sparse.caliper <- match_on(x =
predicted, z = DATA$Z, caliper = 1))

gc()
benchmark.sparse.within <- profr(result.sparse.within <- match_on(x =
predicted, z = DATA$Z, within = exactMatch(Z ~ I(X3 == "a" | X3 == "b"), data = DATA)))

gc()
benchmark.formula <- profr(result.formula <- match_on(x =
  Z ~ X1 + X2 + X3, data = DATA))

gc()
benchmark.glm <- profr(result.glm <- match_on(x =
  glm(Z ~ X1 + X2 + X3, data = DATA)))


save(file = "distance.rda",
  result.dense,
  benchmark.dense,
  result.sparse.caliper,
  benchmark.sparse.caliper,
  result.sparse.within,
  benchmark.sparse.within,
  result.formula,
  benchmark.formula,
  result.glm,
  benchmark.glm)
