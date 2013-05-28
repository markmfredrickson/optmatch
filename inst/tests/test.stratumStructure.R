################################################################################
# Tests for utility functions
################################################################################

if(require("testthat")) {
context("StratumStructure tests")

test_that("Basics", {
  m <- matrix(c(1,Inf,Inf,Inf, 1,Inf,Inf,Inf, Inf,1,Inf,Inf, Inf,Inf,1,Inf, Inf,Inf,Inf,1), nrow = 4,
              dimnames = list(LETTERS[1:4], letters[22:26]))

  fm <- fullmatch(m)
  res.ss <- stratumStructure(fm)

  # current implementation has res.ss as an array with attributes
  # casting to numeric is just easy to to test.
  expect_equal(as.numeric(res.ss), c(3,1))

})
  
}

test_that("Helper functions to compute small bits about matches", {
  Z <- rep(c(0,1), 8)
  B <- rep(c("foo", "bar"), each = 8)
  names(Z) <- names(B) <- letters[1:16]

  # effective sample size should be equal 8
  res.pairs <- pairmatch(exactMatch(Z ~ B), data = Z)  
  expect_equal(effectiveSampleSize(res.pairs), 8)

  expect_equal(effectiveSampleSize.default(res.pairs, Z), 8)
  # effective sample size should be 2/(1/4)
  res.mixed <- fullmatch(exactMatch(Z~B), max.controls = 2, omit.fraction = 1/4, data = Z)
  expect_equal(sum(!is.na(res.mixed)), 14) 
  
  # should throw away 2 controls and create 2 groups with 2 treated, 6 groups total
  # the sum of the harmonic means is: 4 * 2/(1/1 + 1/1) + 2 * 2 / (1/2 + 1/1) = 20/3
  expect_equal(effectiveSampleSize(res.mixed), 20/3)

  #todo: pass a Z argument, say after reordering a match, an error if contrast.group is not found
  # ok, I actually found it really hard to break up the optmatch data, but ust in case it happens
  tmp <- as.numeric(res.pairs)
  names(tmp) <- names(res.pairs)
  expect_error(effectiveSampleSize.factor(tmp), "contrast.group")
  expect_error(effectiveSampleSize(tmp))
  
  expect_equal(effectiveSampleSize(tmp, Z), 8)
})


test_that("Correct output for full failures", {
  expect_equal(effectiveSampleSize(rep(NA, 10), rep(c(T,F), 5)), 0)

})
