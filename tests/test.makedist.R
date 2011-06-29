################################################################################
# Tests for makedist function
################################################################################

library(testthat)
context("Makedist tests")

test_that("Checking input", {
  # Z should have exactly two levels
  expect_error(makedist(rep(1,10), rep(1,10), identity))
  expect_error(makedist(rep(0,10), rep(1,10), identity))
  expect_error(makedist(rep(c(1,2,3), 3), rep(1,9), identity))

  # Z and data should be same length
  expect_error(makedist(rep(c(1,0), 5), c(1,2,3), identity))
  expect_error(makedist(rep(c(1,0), 5), data.frame(c(1,2,3), c(4,5,6))))

  # no NA's in Z
  expect_error(makedist(c(NA, 1, 0, 1, 0), c(1,2,3,4,5), identity))

  # Z and/or data should have rownames
  expect_error(makedist(c(rep(1, 5), rep(0, 5)), 1:10, `-`))
})

test_that("No mask => dense matrix", {
  data <- c(1:5, 2:9)
  names(data) <- letters[1:13]
  z <- c(rep(0, 5), rep(1, 8))

  # this is what we should get. the equivalent of outer
  m <- t(outer(X = data[z == 1], Y = data[z == 0], FUN = `-`))

  res <- makedist(z, data, `-`)

  expect_equal(dim(res), c(5, 8))
  expect_is(res, "matrix")
  expect_equivalent(res, m)

  # same basic test, with a data frame
  data.df <- data.frame(a = data, xyz = 1:13)
  aminus <- function(treat, control) { treat$a - control$a }

  res.df <- makedist(z, data.df, aminus)
  expect_equal(res.df, res)
  expect_equivalent(res.df, m)
})
