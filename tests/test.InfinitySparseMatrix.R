################################################################################
### Tests for the InfinitySparseMatrix class
###############################################################################
library(testthat)

context("InfinitySparseMatrix tests")

test_that("ISM Basics", {
  A <- new("InfinitySparseMatrix", c(1,2,3), cols = c(1,2, 2), rows = c(1,2,4))  
  expect_is(A, "InfinitySparseMatrix")
  expect_equal(dim(A), c(2,2))

})

