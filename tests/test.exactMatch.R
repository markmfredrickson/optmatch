################################################################################
# Tests for exactMatch function: a function to create InfinitySpareMatrices
################################################################################

library(testthat)

context("exactMatch function")

test_that("Exact Match on Factors", {
  n <- 16
  Z <- rep(c(0,1), n/2)
  my.names <- paste(rep(c("C", "T"), n/2), 1:16, sep = "")
  names(Z) <- my.names

  W <- rnorm(16)
  B <- c(rep(0, n/2), rep(1, n/2))
  test.data <- data.frame(Z, W, B)

  res <- exactMatch(B, Z) # factor, factor implementation

  # the resulting matrix should be block diagonal
  m0 <- matrix(0, nrow = n/4, ncol = n/4)
  mInf <- matrix(Inf, nrow = n/4, ncol = n/4)

  tmp1 <- cbind(m0, mInf)
  tmp2 <- cbind(mInf, m0)
  m <- rbind(tmp1, tmp2)

  expect_equivalent(as.matrix(res), m)
  expect_equal(dim(res), c(8,8))

  expect_error(exactMatch(B, rep(1:(n/4), 4)))
  expect_error(exactMatch(B, c(Z, 0)))
  expect_error(exactMatch(c(B, 1), Z))

  # row and column names
  expect_equal(rownames(res), my.names[Z == 1])
  expect_equal(colnames(res), my.names[Z == 0])
  

})
