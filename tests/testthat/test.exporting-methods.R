# the testthat library executes test within the optmatch namespace,
# so it can't detect if we forget to export methods
# R CMD check tests, on the other hand, use the package externally

context('exporting methods')

test_that("sparse", {
  data(nuclearplants)
  tmp <- match_on(pr ~ date + cost, data = nuclearplants, within = exactMatch(pr ~ pt, data = nuclearplants))

  tmp.m <- as.matrix(tmp)

  expect_equal(dim(tmp.m), c(10,22))
  expect_equal(class(tmp.m), "matrix")
})

test_that("dense", {
  data(nuclearplants)
  tmp <- match_on(pr ~ date + cost, data = nuclearplants)

  tmp.m <- as.matrix(tmp)

  expect_equal(dim(tmp.m), c(10,22))
  expect_equal(class(tmp.m), "matrix")
})
