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
