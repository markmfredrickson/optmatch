################################################################################
### Tests for the InfinitySparseMatrix class
###############################################################################
library(testthat)

context("InfinitySparseMatrix tests")

test_that("ISM Basics", {
  A <- makeInfinitySparseMatrix(c(1,2,3), cols = c(1,2, 2), rows = c(1,1,2))
  expect_is(A, "InfinitySparseMatrix")
  expect_equal(dim(A), c(2,2))
  
  # converting to the equivalent matrix
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2)
  expect_equivalent(as.matrix(A), m)

  # converting from a matrix to a ISM
  expect_equivalent(as.InfinitySparseMatrix(m), A)
  # and back again
  expect_equivalent(as.matrix(as.InfinitySparseMatrix(m)), m)
  # expect_equal(as(m, "InfinitySparseMatrix"), A)

  # a more complicated examples, missing an entire row/col
  w <- matrix(c(1,Inf,2, 3, Inf, 4), nrow = 3)
  B <- as.InfinitySparseMatrix(w)
  expect_equivalent(as.matrix(B), w)

  y <- matrix(c(1,2,3,Inf, Inf, Inf), nrow = 3)
  D <- as.InfinitySparseMatrix(y)
  expect_equivalent(as.matrix(D), y)

  # the as() technique should be equivalent
  expect_equivalent(as(D, "matrix"), y)
  expect_equivalent(A, as(m, "InfinitySparseMatrix"))

})

test_that("ISM Handles Names", {
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2)
  rownames(m) <- c("A", "B")
  colnames(m) <- c("C", "D")

  expect_equal(as.matrix(as(m, "InfinitySparseMatrix")), m)
  
})
