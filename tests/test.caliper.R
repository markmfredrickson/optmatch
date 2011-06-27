################################################################################
# Caliper Tests
################################################################################

library(testthat)
context("Caliper")

test_that("Caliper return values", {
  n <- 16
  Z <- c(rep(0, n/2), rep(1, n/2))
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- rep(c(0,1), n/2)
  test.data <- data.frame(Z, X1, X2, B)

  # use the Mahalanobis distance mdist method
  result <- caliper(.2, Z ~ X1 + X2, data = test.data)
  expect_is(result, "optmatch.dlist")
  
})

test_that("Caliper exclusion", {
  n <- 16
  Z <- c(rep(0, n/2), rep(1, n/2))
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- rep(c(0,1), n/2)
  test.data <- data.frame(Z, X1, X2, B)
  rownames(test.data) <- letters[1:n]
  # test exclusion functionality
  result <- caliper(.2, Z ~ X1 + X2, data = test.data, exclude = c("a", "b"))
  m <- result$m
  expect_true(all(m[, "a"]) == 0)
  expect_true(all(m[, "b"]) == 0)
  
})

