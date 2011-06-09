################################################################################
### Tests for the InfinitySparseMatrix class
###############################################################################
library(testthat)

context("InfinitySparseMatrix tests")

test_that("ISM Basics", {
  A <- new("InfinitySparseMatrix", c(1,2,3), cols = c(1,2, 2), rows = c(1,3,4))  
  expect_is(A, "InfinitySparseMatrix")
  expect_equal(dim(A), c(2,2))

  # converting to the equivalent matrix
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2)
  expect_equal(as.matrix(A), m)

})

