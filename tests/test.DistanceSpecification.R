################################################################################
### Tests for objects implementing the distance specification protocol
################################################################################
library(testthat)
library(SparseM)

context("Distance Specification Protocol")

test_that("Basic types are DistanceSpeficiations", {
  expect_true(isClassUnion("DistanceSpecification"))

  m <- matrix(1, nrow = 2, ncol = 2)
  expect_true(is(m, "DistanceSpecification")) # expect_is fails for S3 classes?

  m.csr <- as.matrix.csr(m)
  expect_is(m.csr, "DistanceSpecification")

  # make up a optmatch.dlist
  # odl <- mdist(pr ~ t1 + t2, data = nuclearplants)
  # expect_is(odl, "DistanceSpecification")
  
})

test_that("DistSpec => nodes and arcs", {
  # tests:
  expect_true(isGeneric("prepareMatching"))
  
  # matrix and matrix.csr implement prepareMatching
  # each returns the proper number of nodes and arcs
  
  # test data: 2t + 3c = 5 nodes, 4 arcs (2 pairs unmatchable)
  m <- matrix(c(1, Inf, 1, 2, 2, Inf), nrow = 2, ncol = 3)
  # for the sparse representation, we use 0 to notate unmatchable
  m.csr <- as.matrix.csr(matrix(c(1, 0, 1, 2, 2, 0), nrow = 2, ncol = 3))

})
