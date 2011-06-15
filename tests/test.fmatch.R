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

  expect_equal(length(res), 6)
  expect_equal(length(levels(res)), 3)
  expect_equal(res["A"], res["D"])
  expect_false(res["A"] == res["B"])

  ### for use at the R prompt
  test.data <- prepareMatching(m)
  
})

test_that("Fortran output => factor", {
  m <- matrix(c(1, Inf, 2,
                2, 1, Inf,
                3, 2, 1), nrow = 3, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E", "F")
  test.data <- prepareMatching(m)
  matches <- c(1,0,1,0,0,1,1)
  # AB - D, C - EF

  res <- matches2factor(test.data, matches)

  expect_equal(length(res), 6)
  expect_equal(length(levels(res)), 2)

  # make sure the correct items are paired
  # this is probably overkill, but better safe than sorry
  expect_equivalent(res["A"], res["B"])
  expect_equivalent(res["A"], res["D"])
  expect_equivalent(res["B"], res["D"])
  expect_equivalent(res["C"], res["E"])
  expect_equivalent(res["C"], res["F"])
  expect_equivalent(res["E"], res["F"])

  expect_false(res["A"] == res["C"])

  expect_equal(c("A", "B", "C", "D", "E", "F"), names(res))

})
