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

test_that("Math Ops", {
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2)
  A <- as.InfinitySparseMatrix(m)

  # scalar math
  expect_equivalent(as.matrix(A + 1), m + 1)
  expect_equivalent(as.matrix(A - 1), m - 1)
  expect_equivalent(as.matrix(A * 2), m * 2)
  expect_equivalent(as.matrix(A / 2), m / 2)
  
  # matrix element wise math
  expect_equivalent(as.matrix(A + A), m + m)

  # Inf - Inf or Inf / Inf gives NA
  mm <- m - m
  mm[is.na(mm)] <- Inf

  md <- m / m
  md[is.na(md)] <- Inf

  expect_equivalent(as.matrix(A - A), mm)
  expect_equivalent(as.matrix(A * A), m * m)
  expect_equivalent(as.matrix(A / A), md)

  # The harder case is when the matrix has non-identical row/col ids

  q <- matrix(c(1, 2, Inf, 4), nrow = 2, ncol = 2)
  B <- as.InfinitySparseMatrix(q)

  expect_equivalent(as.matrix(A + B), m + q)
  expect_equivalent(as.matrix(A * B), m * q)

  # TODO, make up temp matrices for sub and div

  # dense + sparse => sparse
  Aq = A + q
  expect_is(Aq, "InfinitySparseMatrix")
  expect_equivalent(as.matrix(Aq), m + q)

  # make sure it works the other direction (and with mult)
  qA = q * A
  expect_is(qA, "InfinitySparseMatrix")
  expect_equivalent(as.matrix(qA), q * m)

  # names should be sticky across arithmatic
  # TODO, math should reorder by names in case that changes things
  colnames(A) <- paste("C", 1:2, sep = "")
  rownames(A) <- paste("T", 1:2, sep = "")
  colnames(q) <- paste("C", 1:2, sep = "")
  rownames(q) <- paste("T", 1:2, sep = "")

  Aq = A + q
  expect_equal(colnames(Aq), c("C1", "C2"))
  expect_equal(rownames(Aq), c("T1", "T2"))

})
