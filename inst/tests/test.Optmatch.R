################################################################################
# Tests for the optmatch object and basic methods
################################################################################

library(testthat)
context("Optmatch object")

test_that("Object creation", {
  dist <- diag(5)
  dimnames(dist) <- list(letters[1:5], letters[6:10])
  ms <- list(c(a = 1, f = 1, b = 1 , g = 2), c(c = 1, h = 1, d = 2, i = 2, e = 3, j = 3))  
  res.opt <- makeOptmatch(dist, ms, NULL)

  expect_equal(length(res.opt), 10)
  expect_is(res.opt, "factor")
  expect_is(res.opt, "optmatch")
  
  # two levels of matches shouldn't be 1.NA, 2.NA, should just be NA
  ms2 <- list(c(a = 1, b = 2, c = 1 , d = NA), c(e = 1, f = 2, g = 3, h = 1, i = 2, j = NA))
  res.opt2 <- makeOptmatch(dist, ms2, NULL)
  expect_true(all(is.na(res.opt2[c("d", "j")])))

})

test_that("Object subsetting", {
  dist <- diag(5)
  dimnames(dist) <- list(letters[1:5], letters[6:10])
  ms <- list(c(a = 1, f = 1, b = 1 , g = 2), c(c = 1, h = 1, d = 2, i = 2, e = 3, j = 3))  
  res.opt <- makeOptmatch(dist, ms, NULL)

  expect_equal(names(res.opt[1:4]), c("a", "f", "b", "g"))
  expect_equal(length(res.opt[c("a", "b")]), 2)
  
})

test_that("Matched distances", {
  # see R/matched.distances.R for the function
  # it is only called by makeOptmatch internally, so putting the tests here
  # start with an easy case: 
  dist <- matrix(Inf, nrow = 5, ncol = 5)
  diag(dist) <- 1:5
  dimnames(dist) <- list(letters[1:5], letters[6:10])
  dist.match <- as.factor(c(1.1,1.1,1.2,1.2,2.1,2.1,2.2,2.2,2.3,2.3))
  names(dist.match) <- c("a","f","b","g","c","h","d","i","e","j")
  class(dist.match) <- c("optmatch", "factor")

  res.md <- matched.distances(dist.match, dist)
  expect_equivalent(as.vector(res.md), 1:5)

  # now an ISM version
  dist.i <- as.InfinitySparseMatrix(dist)
  res.mdi <- matched.distances(dist.match, dist.i)
  expect_equivalent(as.vector(res.mdi), 1:5)
  
})
