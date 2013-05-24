load("setup.rda")
library(profr)
library(optmatch)

gc()
# we have to use a smaller interval to actually catch the timings on this fast
# operation.
benchmark.mdist <- profr(result.mdist <- mdist(model, 
 structure.fmla = ~ X3), interval = 0.001)

save(file = "mdist.rda",
     benchmark.mdist,
     result.mdist)
