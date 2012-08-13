################################################################################
# Tests for interaction with RItools
################################################################################

if(require("testthat") && require("RItools")) {
context("RItools interaction")

test_that("Summary function adds RItools info", {

  ### Testing ways to create glm objects and send them to RItools for balance tests

  test.data <- data.frame(Z = rep(c(0,1),10), X = c(1,2,NA,3,4,NA,5,NA,9,10))

  ## Method 1: passing a matrix/data.frame or similar

  # fill.NAs is a layer of indirection between the glm object and RItools
  # putting th test here because this could happen for other ways to create glms, not just fill.NAs
  test.df <- glm(fill.NAs(Z ~ X, data = test.data), family = binomial)

  test.pm <- pairmatch(mdist(test.df))

  # may cause an error if the function can't get at the proper data that was used to create the glm object
  res <- summary(test.pm, test.df)

  expect_true(!is.null(res$balance))

  ## Method 2: using a formula and a data argument
  
  test.fmla <- glm(Z ~ X, data = test.data, family = binomial)
  test.pm <- pairmatch(mdist(test.fmla))

  res <- summary(test.pm, test.fmla)
  expect_true(!is.null(res$balance))

  ## Method 3: Using enviornment
  test.env <- with(test.data, glm(Z ~ X, family = binomial))
  test.pm <- pairmatch(mdist(test.env))

  res <- summary(test.pm, test.env)
  expect_true(!is.null(res$balance))
  
})

}
