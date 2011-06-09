################################################################################
### R/Fortran Interface Tests
################################################################################

library(testthat)
context("R/Fortran Interface")

test_that("fmatch accepts DistanceSpecifications", {
  # the goal of this matrix is that there is a clear match to make
  # A:D, B:E, C:F
  m <- matrix(c(1, Inf, 2,
                2, 1, Inf,
                3, 2, 1), nrow = 3, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E", "F")

  res <- fmatch(m, 3, 3)

  # expect_equal(length(res), 3)
  
})
