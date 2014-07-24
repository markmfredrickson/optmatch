################################################################################
# Reversion tests for specific issues.
################################################################################

context("Reversion Tests")

test_that("Issue 82", {
  data(nuclearplants)

  # this would fail in the old optmatch
  x <- match_on(cap, z = pr, data = nuclearplants)
  expect_is(x, "DenseMatrix")

  # now test if we have an actual variable called cap
  
  n <- dim(nuclearplants)[1]
  cap <- rep(1, n)
  names(cap) <- 1:n
  
  y <- match_on(cap, z = pr, data = nuclearplants)
  
  expect_true(all(y == 0))
})
