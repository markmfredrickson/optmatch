################################################################################
# Tests for utility functions
################################################################################

context("Utility Functions")

test_that("toZ", {

  # Per discussion on issue #124, toZ will only accept numeric 0/1/NA or
  # logical TRUE/FALSE/NA.

  Z <- rep(c(TRUE, FALSE), 5) # correct representation
  names(Z) <- letters[1:10]

  expect_identical(Z, toZ(Z))

  nZ <- as.numeric(Z)
  names(nZ) <- names(Z)
  expect_identical(Z, toZ(nZ))

  Zmiss <- Z
  Zmiss[1] <- NA
  expect_identical(Zmiss, toZ(Zmiss))

  # expected failures
  fZ <- as.factor(Z)
  names(fZ) <- names(Z)
  expect_error(toZ(fZ))

  cZ <- as.character(Z)
  names(cZ) <- names(Z)
  expect_error(toZ(cZ))

  expect_error(toZ(c(1,2,3,4)))
  expect_error(toZ(as.factor(c(1,2,3,4))))
  expect_error(toZ(as.character(c(1,2,3,4))))

  expect_error(toZ(c(0,0,1,1,2)))
  expect_error(toZ(c(0,0,2,2,2)))
  expect_error(toZ(c(0,0)))
  expect_error(toZ(c(1,1)))

  # single column data.frames and matrices should be valid
  dZ <- data.frame(Z)
  expect_identical(Z, toZ(dZ))

  mZ <- matrix(Z, ncol = 1) ; rownames(mZ) <- names(Z)
  expect_identical(Z, toZ(mZ))

})

test_that("#159 - toZ for labelled", {
  if (requireNamespace("haven", quietly = TRUE)) {
    Z <- haven::labelled(c(0,1,0), c(Treatment = 1, Control = 0))
    expect_is(toZ(Z), "logical")
    expect_length(Z, 3)
    Znum <- as.numeric(Z)
    expect_identical(toZ(Z), toZ(Znum))
  } else {
    expect_true(TRUE) # avoiding empty test warning
  }
})
