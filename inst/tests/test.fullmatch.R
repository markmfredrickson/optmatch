################################################################################
# Fullmatch tests
################################################################################

library(testthat)
library(optmatch)

context("fullmatch function")

test_that("No cross strata matches", {
  # test data
  Z <- rep(c(0,1), 4)
  B <- rep(c(0,1), each = 4)
  distances <- 1 + exactMatch(Z ~ B)

  res <- fullmatch(distances)
  expect_false(any(res[1:4] %in% res[5:8]))
})
