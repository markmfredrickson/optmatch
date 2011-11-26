################################################################################
# Caliper Tests
################################################################################

library(testthat)
context("Caliper")

test_that("Caliper return values", {
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2,
              dimnames = list(treated = c("A", "B"),
                              control = c("C", "D")))
  A <- as.InfinitySparseMatrix(m)

  # use the Mahalanobis distance mdist method
  result <- caliper(2, A)
  expect_is(result, "DistanceSpecification")
  
  expect_equal(result@.Data, c(0,0))

  # make sure that matrix input does same thing
  expect_equal(caliper(2, A), caliper(2, m))
})

test_that("Caliper exclusion", {
  m <- matrix(c(3,Inf, 1, 3), nrow = 2, ncol = 2,
              dimnames = list(treated = c("A", "B"),
                              control = c("C", "D")))
  A <- as.InfinitySparseMatrix(m)

  # use the Mahalanobis distance mdist method
  result <- caliper(2, A, exclude = c("B"))


  m2 <- matrix(c(Inf,Inf, 0, 0), nrow = 2, ncol = 2,
              dimnames = list(treated = c("A", "B"),
                              control = c("C", "D")))

  expect_equal(as.matrix(result), m2)
  
})


