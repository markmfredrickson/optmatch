################################################################################
# fill.NAs tests
################################################################################

library(testthat)
context("fill.NAs")

test_that("Basic Tests", {
  # Takes and returns a data frame
  expect_is(fill.NAs(data.frame(1)), "data.frame")
  
  # A formula alone is not allowed
  expect_error(fill.NAs(y ~ x))
 
  sample.df <- data.frame(a = 1:100, b = 100:1, c = rep(c(1,2, NA, 3, 4), 20))

  # takes a formula and a data.frame, returns a data frame
  result <- fill.NAs(a ~ b, sample.df)
  expect_is(result, "data.frame") # no missingness
  
  # simple calls should be equivalent to model.frame
  expect_equal(length(result), 2)
  
  # Adds additional columns for missing data indicators
  result <- fill.NAs(sample.df)
  expect_equal(length(colnames(result)), 4)
  
  # the last column should be TRUE every 3 unit
  expect_identical(result[[4]], rep(c(F, F, T, F, F), 20))
  
  # column name should be c.NA
  expect_identical(colnames(result)[4], "c.NA")
})

test_that("Function expansion", {
  library(splines)
  # for variables encapsulated in functions, only the variable should be expanded into a NA column
  sample.df <- data.frame(a = 1:100, b = 100:1, c = rep(c(1,2, NA, 3, 4), 20))

  result <- fill.NAs(a ~ ns(c, df = 3), sample.df)
  expect_equal(length(result), 5)
  expect_equal(colnames(result)[1], "a")

  ## right number of columns if 2 of the same variable used
  imputed.fmla <- fill.NAs(a ~ log(c) + sqrt(c), data = sample.df)
  expect_equal(dim(imputed.fmla)[2],  4)

})

test_that("Matrices are valid", {
  sample.df <- as.matrix(data.frame(a = 1:100, b = 100:1, c = rep(c(1,2, NA,
  3, 4), 20)))

  result <- fill.NAs(a ~ ns(c, df = 3), sample.df)
  expect_equal(length(result), 5)
  expect_equal(colnames(result)[1], "a")
})

test_that("Results pass to lm()", {
  sample.df <- data.frame(a = 1:100, c = rep(c(1,2, NA, 3, 4), 20))
  
  imputed.fmla <- fill.NAs(a ~ log(c), data = sample.df)
  imputed.frame <- fill.NAs(sample.df)

  m1 <- lm(imputed.fmla)
  m2 <- lm(a ~ log(c) + c.NA, data = imputed.frame)
 
  # for some reason log(c) appears as `log(c)`. I strip these
  # out and treat the results as equal otherwise
  expect_identical(gsub("`", "", names(m1$coef)), names(m2$coef))

})

test_that("Response not imputed by default", {
    
  #### Do not impute response, only covariates
  naresponse.df <- data.frame(Y = c(1, 2, 3, NA, 5), X = c(10, 20, NA, 40, 50))
  imputed.response <- fill.NAs(Y ~ X, naresponse.df)
  expect_true(any(is.na(imputed.response$Y)))
  expect_true(!any(is.na(imputed.response$X)))

  #### Impute when all.covs = T
  
  # formula style
  imputed.all <- fill.NAs(Y ~ X, naresponse.df, all.covs = T)
  expect_true(!any(is.na(imputed.all)))

  # model frame style
  imputed.all <- fill.NAs(naresponse.df, all.covs = T)
  expect_true(!any(is.na(imputed.all)))

})

test_that("Transform, then impute", {
  
  #### Transform, then impute ####
  #### turning off tests for now. the strategy is to use model.matrix before
  #### imputing
  transform.df <- data.frame(Y = c(1,2,3,4,5), X1 = c(2,2,4, NA, 4), X2 = c(NA, 10, 20, 30, NA))
  imputed.transform <- fill.NAs(Y ~ X1 * X2, data = transform.df)
  # should have 6 columns Y, X1, X2, X2:X3, X1.NA, and X2.NA
  expect_equal(dim(imputed.transform)[2], 6) 
  expect_identical(imputed.transform$X1 , c(2,2,4,3,4))
  expect_identical(imputed.transform$X2 , c(20, 10, 20, 30, 20))
  expect_equal(imputed.transform[["X1:X2"]], c(50, 20, 80, 50, 50))
  
  i2.transform <- fill.NAs(Y ~ X1, data = transform.df)
  expect_equal(length(i2.transform), 3)
})


