################################################################################
# Tests for interaction with RItools
################################################################################

if(require("testthat") && require("RItools")) {
context("RItools interaction")

test_that("Summary function adds RItools info", {

  ### Testing ways to create glm objects and send them to RItools for balance tests

  test.data <- data.frame(Z = rep(c(0,1),10), X = 1:10)
  test.data.na <- test.data ; test.data.na$X[c(3,6,8)] <- NA

  ## Method 1: passing a matrix/data.frame or similar

  # fill.NAs is a layer of indirection between the glm object and RItools
  # putting th test here because this could happen for other ways to create glms, not just fill.NAs
  test.df <- glm(fill.NAs(Z ~ sin(X), data = test.data.na), family = binomial)

  match.df <- pairmatch(match_on(test.df), data = test.data.na)

  # may cause an error if the function can't get at the proper data that was used to create the glm object
  res <- summary(match.df, test.df)

  expect_true(!is.null(res$balance))

  ## Method 2: using a formula and a data argument
  
  test.fmla <- glm(Z ~ sin(X), data = test.data, family = binomial)
  match.fmla <- pairmatch(match_on(test.fmla), data = test.data)

  res <- summary(match.fmla, test.fmla)
  expect_true(!is.null(res$balance))

  ## Method 3: Using enviornment
  test.env <- with(test.data, glm(Z ~ sin(X), family = binomial))
  match.env <- pairmatch(match_on(test.env), data = test.data)

  res <- summary(match.env, test.env)
  expect_true(!is.null(res$balance))

  ## Now trying Method 2 with mising values

  test.fna <- glm(Z ~ sin(X), data = test.data.na, family = binomial)
  match.fna <- pairmatch(match_on(test.fna), data = test.data.na)

  res <- summary(match.fna, test.fna)
  expect_true(!is.null(res$balance))

  ## Finally, we should issue a useful error message when users do not pass the data argument.
  match.nodata <- pairmatch(match_on(test.fna)) # data not passed

  expect_error(summary(match.nodata, test.fna), "'data' argument")
  
})

}
