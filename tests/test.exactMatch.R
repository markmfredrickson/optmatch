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

  res <- exactMatch(B, treatment = Z) # factor, factor implementation

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

test_that("Exact match on formula", {
  n <- 16
  Z <- rep(c(0,1), n/2)
  my.names <- paste(rep(c("C", "T"), n/2), 1:16, sep = "")
  names(Z) <- my.names

  W <- rnorm(16)
  B <- c(rep(0, n/2), rep(1, n/2))
  test.data <- data.frame(Z, W, B)

  res <- exactMatch(Z ~ B)

  # the resulting matrix should be block diagonal
  m0 <- matrix(0, nrow = n/4, ncol = n/4)
  mInf <- matrix(Inf, nrow = n/4, ncol = n/4)

  tmp1 <- cbind(m0, mInf)
  tmp2 <- cbind(mInf, m0)
  m <- rbind(tmp1, tmp2)

  expect_equivalent(as.matrix(res), m)
  expect_equal(dim(res), c(8,8))

  res.data <- exactMatch(Z ~ B, data = test.data)
  expect_equal(res.data, res)

  # combine mulitiple factors into a single factor
  B2 <- rep(c(0,1), 4, each = 2)

  # combine them by hand into a single factor
  BB <- B + 2 * B2
  res.bb <- exactMatch(BB, Z)

  res.multi <- exactMatch(Z ~ B + B2)

  expect_equal(res.bb, res.multi)
  
})

test_that("Use proper environment or data.frame", {
  n <- 16
  Z <- rep(c(0,1), n/2)
  W <- rnorm(16)
  B <- c(rep(0, n/2), rep(1, n/2))
  test.data <- data.frame(a = Z, x = W, c = B)

  names(Z) <- letters[1:n]
  rownames(test.data) <- letters[1:n]

  res.envir <- exactMatch(Z ~ B)
  res.df <- exactMatch(a ~ c, data = test.data)

  expect_identical(res.envir, res.df)
  
})

test_that("Makes correct mask", {
  # this data gave me problems with a makedist() test.
  # it should produces a matrix with a 2x3 0 matrix in
  # the upper left and a 3x2 0 m matrix in the lower right
  # it was producing a 3x3 and a 2x2 for some reason.

  set.seed(20110629)
  data <- data.frame(z = rep(c(1,0), 5),
                     y = rnorm(10),
                     b = rep(c(1,0), each = 5))
  rownames(data) <- letters[1:10]
  Y <- data$z
  A <- data$b
  names(Y) <- rownames(data)
  names(A) <- rownames(data)

  reference <- matrix(c(rep(c(0,0,0,Inf,Inf), 2), 
                        rep(c(Inf, Inf, Inf, 0, 0), 3)), 
                      nrow = 5, ncol = 5)

  mask.df <- exactMatch(z ~ b, data = data)
  expect_equal(length(mask.df), 3*2 + 2*3) # sizes of the 0 blocks

  mask.fac <- exactMatch(A, Y)
  expect_equal(length(mask.fac), 12)

  expect_identical(mask.df, mask.fac)
  
})
