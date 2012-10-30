################################################################################
# Tests for makedist function
################################################################################

library(testthat)
context("Makedist tests")

# this is inefficient (in that it duplicates data entries many times)
# but good enough for small tests
absdiff <- function(index, data, z) {
  abs(data[index[,1]] - data[index[,2]])
}

test_that("Checking input", {
  # Z should have exactly two levels
  expect_error(makedist(rep(1,10), rep(1,10), identity), "Treatment")
  expect_error(makedist(rep(0,10), rep(1,10), identity), "Treatment")
  expect_error(makedist(rep(c(1,2,3), 3), rep(1,9), identity), "Treatment")

  # Z and data should be same length
  expect_error(makedist(rep(c(1,0), 5), c(1,2,3), identity), "length")
  expect_error(makedist(rep(c(1,0), 5), data.frame(c(1,2,3), c(4,5,6))),
               "length")

  # no NA's in Z
  expect_error(makedist(c(NA, 1, 0, 1, 0), c(1,2,3,4,5), identity), "NA")

  # Z and/or data should have rownames
  expect_error(makedist(c(rep(1, 5), rep(0, 5)), 1:10, absdiff), "names")
})

test_that("No within => dense matrix", {
  data <- c(1:5, 2:9)
  names(data) <- letters[1:13]
  z <- c(rep(0, 5), rep(1, 8))

  # this is what we should get. the equivalent of outer
  m <- abs(outer(X = data[z == 1], Y = data[z == 0], FUN = `-`))

  res <- makedist(z, data, absdiff)

  expect_equal(dim(res), c(8, 5))
  expect_is(res, "matrix")
  expect_equivalent(res@.Data, m)

  # same basic test, with a data frame
  data.df <- data.frame(a = data, xyz = 1:13)
  aminus <- function(index, data, z) { abs(data[index[,1], "a"] - data[index[,2], "a"]) }

  res.df <- makedist(z, data.df, aminus)
  expect_equal(res.df, res)
  expect_equivalent(as.matrix(res.df), m)
})

test_that("Mask => ISM result", {
  set.seed(20110629)
  data <- data.frame(z = rep(c(1,0), 5),
                     y = rnorm(10),
                     b = rep(c(1,0), each = 5))
  rownames(data) <- letters[1:10]
 
  yminus <- function(index, data, z) { data[index[,1], 'y'] - data[index[,2], 'y'] }

  upper.left <- makedist(data$z[data$b == 1], data[data$b == 1,], yminus)
  lower.right <- makedist(data$z[data$b == 0], data[data$b == 0,], yminus)
  upper <- cbind(upper.left, matrix(Inf, nrow = 3, ncol = 3))
  lower <- cbind(matrix(Inf, nrow = 2, ncol = 2), lower.right)
  m <- rbind(upper, lower)

  test.within <- exactMatch(z ~ b, data = data)

  res <- makedist(data$z, data, yminus, within = test.within)

  expect_equal(length(res), length(test.within))

  expect_equivalent(as.matrix(res), m)

  # withins should match the data on treatment and control names
  data2 <- data
  rownames(data2) <- letters[11:20]
  test.within.bad <- exactMatch(z ~ b, data = data2)
  
  expect_error(makedist(data$z, data, yminus, within = test.within.bad),
               "names")

  # repeat previous test with bad row and column names respectively
  data3 <- data
  rownames(data3) <- c("foo", rownames(data[-1,]))
  test.within.bad.treat <- exactMatch(z ~ b, data = data3)
  
  expect_error(makedist(data$z, data, yminus, within = test.within.bad.treat),
               "names")

  data4 <- data
  rownames(data3) <- c(rownames(data)[1:9], "bar")
  test.within.bad.cntrl <- exactMatch(z ~ b, data = data3)
  
  expect_error(makedist(data$z, data, yminus, within = test.within.bad.cntrl),
               "names")

})

test_that("makedist works on single column data.frames", {
  set.seed(20110707)
  data <- data.frame(z = rep(c(1,0), 5),
                     y = rnorm(10),
                     b = rep(c(1,0), each = 5))
  rownames(data) <- letters[1:10]
  f <- function(index, data, z) { abs(as.vector(data[index[,1], 1] - data[index[,2], 1])) }

  res <- makedist(data$z, subset(data, T, select = 2), f)
  expect_true(all(res != 0)) # makes sure res <- ... worked
})

test_that("Z can be a numeric, logical, or two level factor", {
  n <- 16
  Z <- numeric(n)
  Z[sample.int(n, n/2)] <- 1
  X1 <- rnorm(n, mean = 5)
  
  names(X1) <- letters[1:n]

  res.one <- makedist(Z, X1, absdiff) 
  res.logical <- makedist(as.logical(Z), X1, absdiff) 
  expect_identical(res.one, res.logical)
  
  res.factor <- makedist(as.factor(Z), X1, absdiff) 
  expect_identical(res.one, res.factor)
  
  Y <- rep(1:4, n/4)
  expect_error(makedist(as.factor(Y), X1, absdiff), "Treatment")
})

test_that("distancefn specification", {
  Z <- rep(c(1,0), 10)
  X <- rep(c(5,10), 10)
  names(Z) <- names(X) <- letters[1:20]

  expect_true(all(makedist(Z, X, absdiff) == 5))

  # should get all args
  absdiffz <- function(index, data, z) {
    if(missing(z)) { stop("Missing z")}
    abs(data[index[,1]] - data[index[,2]])
  }

  expect_true(all(makedist(Z, X, absdiffz) == 5))
  
})

