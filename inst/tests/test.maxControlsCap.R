################################################################################
# maxControlsCap: finding the best input for fullmatch
################################################################################

library(testthat)

context("maxControlsCap function")

test_that("basics", {
  # needs to take both simple and stratified problems.
  n <- 16
  Z <- rep(c(0,1), n/2)
  my.names <- paste(rep(c("C", "T"), n/2), 1:16, sep = "")
  names(Z) <- my.names
  B <- c(rep(0, n/2), rep(1, n/2))

  em <- exactMatch(B, treatment = Z) # factor, factor implementation

  res <- maxControlsCap(em)
  
  expect_equal(length(res$strictest.feasible.max.controls), 2) # two level problem
  expect_true(!is.na(res$strictest.feasible.max.controls))
})
