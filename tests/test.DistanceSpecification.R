################################################################################
### Tests for objects implementing the distance specification protocol
################################################################################
library(testthat)

context("Distance Specification Protocol")

test_that("Basic types are DistanceSpeficiations", {
  expect_true(isClassUnion("DistanceSpecification"))

  m <- matrix(1, nrow = 2, ncol = 2)
  expect_true(is(m, "DistanceSpecification")) # expect_is fails for S3 classes?

  
})

test_that("DistSpec => nodes and arcs", {
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
