################################################################################
### Tests for objects implementing the distance specification protocol
################################################################################
library(testthat)

context("Distance Specification Protocol")

test_that("Basic types are DistanceSpeficiations", {
  expect_true(isClassUnion("DistanceSpecification"))

  m <- matrix(1, nrow = 2, ncol = 2)
  expect_true(is(m, "DistanceSpecification")) # expect_is fails for S3 classes?
  expect_true(is(as.InfinitySparseMatrix(m), "DistanceSpecification"))
  B <- rep(c(0,1), each = 5)
  names(B) <- letters[1:10] # either B or Z must have names
  expect_true(is(exactMatch(B, rep(c(0,1), 5)),
    "DistanceSpecification"))
  
})

test_that("Matrix => nodes and arcs", {
  # tests:
  expect_true(isGeneric("prepareMatching"))
  
  # matrix and matrix.csr implement prepareMatching
  # each returns the proper number of nodes and arcs
  
  # test data: 4 arcs (2 pairs unmatchable)
  m <- matrix(c(1, Inf, 1, 2, 2, Inf), nrow = 2, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E")

  m.result <- prepareMatching(m)
  expect_equal(dim(m.result), c(4, 3))
  expect_equal(unique(m.result$treatment), as.factor(c("D", "E")))
  expect_equal(unique(m.result$control), as.factor(c("A", "B", "C")))


})

test_that("ISM => nodes and arcs", {
  m <- matrix(c(1, Inf, 1, 2, 2, Inf), nrow = 2, ncol = 3)
  A.nonames <- as.InfinitySparseMatrix(m)

  expect_error(prepareMatching(A.nonames))

  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E")
  A <- as.InfinitySparseMatrix(m)

  res.ISM <- prepareMatching(A)

  expect_equal(dim(res.ISM), c(4, 3))
  expect_equal(unique(res.ISM$treatment), as.factor(c("D", "E")))
  expect_equal(unique(res.ISM$control), as.factor(c("A", "B", "C")))

})

test_that("Subproblems", {
  m <- matrix(1, nrow = 2, ncol = 2)
  expect_false(subproblems(m))
  expect_false(subproblems(as.InfinitySparseMatrix(m)))
  
  B <- rep(c(0,1), each = 5)
  names(B) <- letters[1:10]
  em <- exactMatch(B, rep(c(0,1), 5))
  res.em <- subproblems(em)
  expect_equal(length(res.em), 2)
  
  expect_true(all(sapply(res.em, function(i) { is(i,
    "DistanceSpecification")})))

  m1 <- matrix(0, nrow = 2, ncol = 3,
    dimnames = list(c("b", "d"), c("a", "c", "e")))

  m2 <- matrix(0, nrow = 3, ncol = 2,
    dimnames = list(c("f", "h", "j"), c("g", "i")))

  expect_equal(as.matrix(res.em[[1]]), m1)
  expect_equal(as.matrix(res.em[[2]]), m2)
})
