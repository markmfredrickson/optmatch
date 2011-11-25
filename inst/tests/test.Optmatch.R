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
  res.opt2 <- makeOptmatch(ms2, c("a", "c", "e", "f", "g"), NULL)
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
