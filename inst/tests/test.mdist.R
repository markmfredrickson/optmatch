################################################################################
# Mdist: methods to create distance specifications, possibliy sparse
################################################################################

library(testthat)

context("mdist function")

test_that("Distances from glms", {
  # note: mdist.glm uses mdist.numeric, so if these tests fail, check that as well.

  n <- 16
  Z <- numeric(n)
  Z[sample.int(n, n/2)] <- 1
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- rep(c(0,1), n/2)

  test.glm <- glm(Z ~ X1 + X2 + B, family = binomial()) # the coefs should be zero or so
   
  result.glm <- mdist(test.glm)
  
  expect_true(is(result.glm, "optmatch.dlist"))
  # can't combine s3 classes in a class union: expect_true(is(result.glm, "DistanceSpecification"))
  expect_equal(length(result.glm), 1) # not stratified
  
})

test_that("Distances from formulas", {
  n <- 16
  Z <- rep(c(0,1), n/2)
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- as.factor(rep(c(0,1), each = n/2))

  test.data <- data.frame(Z, X1, X2, B)

  result.fmla <- mdist(Z ~ X1 + X2 + B, data = test.data)
  expect_true(is(result.fmla, "optmatch.dlist"))

  # test pulling from the environment, like lm does
  result.envir <- mdist(Z ~ X1 + X2 + B)
  expect_equivalent(result.fmla, result.envir)

  expect_error(mdist(~ X1 + X2, data = test.data))
  expect_error(mdist(Z ~ 1, data = test.data))

  # NB: these were written for the S4 match_on function. They fail for the S3 mdist function.
  # checking diferent classes of responses
#  res.one <- mdist(Z ~ X1) 
#  res.logical <- mdist(as.logical(Z) ~ X1)
#  expect_identical(res.one, res.logical)

#  res.factor <- mdist(as.factor(Z) ~ X1)
#  expect_identical(res.one, res.factor)

  # stratifying
  res.strat <- mdist(Z ~ X1 + X2 | B)
  expect_is(res.strat, "optmatch.dlist")
  expect_equal(length(res.strat), 2 )

})

test_that("Distances from functions", {
  n <- 16
  Z <- c(rep(0, n/2), rep(1, n/2))
  X1 <- rep(c(1,2,3,4), each = n/4)
  B <- rep(c(0,1), n/2)
  test.data <- data.frame(Z, X1, B)

  # NB: match_on takes a different kind of function. In that version, treateds and controls
  # are of equal length, one for each treated and control pair in the final matrix (basically, 
  # outer is called before)
  sdiffs <- function(treatments, controls) {
      abs(outer(treatments$X1, controls$X1, `-`))
  }
  
  result.function <- mdist(sdiffs, structure.fmla = Z ~ B, data = test.data)
  expect_equivalent(lapply(result.function, dim), list(c(4,4), c(4,4)))

  # no treatment indicator
  expect_error(mdist(sdiffs, data = test.data))

  # no data
  expect_error(mdist(sdiffs, structure.fmla = Z ~ 1))
  
})
 
test_that("Errors for numeric vectors", {
  expect_error(mdist(1:10))
})

test_that("Bigglm distances", {
  if (require('biglm')) {
    n <- 16
    Z <- c(rep(0, n/2), rep(1, n/2))
    X1 <- rnorm(n, mean = 5)
    X2 <- rnorm(n, mean = -2, sd = 2)
    B <- rep(c(0,1), n/2)
    test.data <- data.frame(Z, X1, X2, B)

    bgps <- bigglm(Z ~ X1 + X2, data = test.data, family = binomial())
    res.bg <- mdist(bgps, structure.fmla = Z ~ 1, data = test.data)

    # compare to glm
    res.glm <- mdist(glm(Z ~ X1 + X2, data = test.data, family = binomial()))
    expect_equivalent(res.bg, res.glm) 

    # structure.fmla arg required
    expect_error(mdist(bgps, data = test.data))
  }
})

test_that("Jake found a bug 2010-06-14", {
  ### Issue appears to be a missing row.names/class

  jb.sdiffs <- function(treatments, controls) {
    abs(outer(treatments$X1, controls$X2, `-`))
  }
  
  n <- 16
  Z <- c(rep(0, n/2), rep(1, n/2))
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- rep(c(0,1), n/2)
  test.data <- data.frame(Z, X1, X2, B)

  absdist1 <- mdist(jb.sdiffs, structure.fmla = Z ~ 1, data = test.data)
  # failing because fmatch is in transition, commentb back in later
  expect_true(length(pairmatch(absdist1)) > 0)
 
})

test_that("General optmatch.dlist tests", {
  # test data: 8 arcs (2 pairs unmatchable in each subgroup)
  m1 <- m2 <- matrix(c(1, Inf, 1, 2, 2, Inf), nrow = 2, ncol = 3)

  colnames(m1) <- c("A", "B", "C")
  rownames(m1) <- c("D", "E")

  colnames(m2) <- c("f", "g", "h")
  rownames(m2) <- c('i', "j")

  odl <- list(m1 = m1, m2 = m2)
  class(odl) <- c("optmatch.dlist", "list")

  expect_equal(dim(odl), c(4,6))

  expect_equal(dimnames(odl), list(treated = c("D", "E", "i", 'j'), control = c("A", "B", "C", "f", "g", "h")))

  matrix.expected <- matrix(c(1,Inf,Inf,Inf, 1,2,Inf,Inf, 2,Inf,Inf,Inf,
                            Inf,Inf,1,Inf, Inf,Inf,1,2, Inf,Inf,2,Inf),
                            nrow = 4, ncol = 6,
                            dimnames = list(treated = c("D", "E", "i", "j"),
                                         control = c("A", "B", "C", "e", "f", "g")))

  expect_equivalent(as.matrix(odl), matrix.expected)

  subset.expected <- matrix(c(1, Inf, 1, Inf, Inf, 2), nrow = 2, ncol = 3,
    dimnames = list(treated = c("D", "i"), control = c("A", "B", "h")))

  expect_equal(subset(odl, subset = c(T,F,T,F), select = c(T,T,F,F,F,T)), subset.expected)

})

test_that("update() optmatch.dlist objects", {
  test.data <- data.frame(Z = rep(c(0,1), 10),
                          X1 = 1:20,
                          B = rep(c(0,1), each = 10))

  basic <- mdist(Z ~ X1, data = test.data, structure.fmla = ~ B)
  expect_equal(length(basic), 2)

  test.data$B <- rep(1:5, each = 4)
  updated <- update(basic)
  expect_equal(length(updated), 5)
})
