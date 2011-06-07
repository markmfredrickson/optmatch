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
