library(profr)
library(optmatch)

load("setup.rda")

gc()

benchmark.dense <- profr(result.dense <- match_on(x = predicted, z = DATA$Z))

save(file = "distance-dense.rda",
     benchmark.dense,
     result.dense)
