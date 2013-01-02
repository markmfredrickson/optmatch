# the testthat library executes test within the optmatch namespace, 
# so it can't detect if we forget to export methods
# R CMD check tests, on the other hand, use the package externally

library(optmatch)

data(nuclearplants)

### as.matrix

# sparse
tmp <- match_on(pr ~ date + cost, data = nuclearplants, within = exactMatch(pr ~ pt, data = nuclearplants))

tmp.m <- as.matrix(tmp)

stopifnot(dim(tmp.m) == c(10,22))
stopifnot(class(tmp.m) == "matrix")

# dense
tmp <- match_on(pr ~ date + cost, data = nuclearplants)

tmp.m <- as.matrix(tmp)

stopifnot(dim(tmp.m) == c(10,22))
stopifnot(class(tmp.m) == "matrix")


