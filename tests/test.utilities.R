################################################################################
# Tests for utility functions
################################################################################

if(require("testthat")) {
context("Utility Functions")

test_that("toZ", {
  Z <- rep(c(T,F), 5) # correct representation  
  names(Z) <- letters[1:10]

  expect_identical(Z, toZ(Z))

  nZ <- as.numeric(Z) ; names(nZ) <- names(Z)  
  expect_identical(Z, toZ(nZ))

  fZ <- as.factor(Z) ; names(fZ) <- names(Z)
  expect_identical(Z, toZ(fZ))

  cZ <- as.character(Z) ; names(cZ) <- names(Z)
  expect_identical(Z, toZ(cZ))

})
  
}
