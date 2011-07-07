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

  # expected failures
  expect_error(toZ(c(1,2,3,4)))
  expect_error(toZ(as.factor(c(1,2,3,4))))
  expect_error(toZ(as.character(c(1,2,3,4))))

  # single column data.frames and matrices should be valid
  dZ <- data.frame(Z)
  expect_identical(Z, toZ(dZ))

  mZ <- matrix(Z, ncol = 1) ; rownames(mZ) <- names(Z)
  expect_identical(Z, toZ(mZ))
})
  
}
