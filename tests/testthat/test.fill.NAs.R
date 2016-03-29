################################################################################
# fill.NAs tests
################################################################################

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
  expect_equal(dim(fill.NAs(sample.df))[2], 4)
  expect_equal(dim(fill.NAs(sample.df, all.covs = T))[2], 4)

  result <- fill.NAs(sample.df)
  # the last column should be TRUE every 3 unit
  expect_identical(result[[4]], rep(c(F, F, T, F, F), 20))

  # column name should be c.NA
  expect_identical(colnames(result)[4], "c.NA")
})

test_that("Function expansion", {
  if (require(splines)) {
    # for variables encapsulated in functions, only the variable should be expanded into a NA column
    sample.df <- data.frame(a = 1:100, b = 100:1, c = rep(c(1,2, NA, 3, 4), 20))

    result <- fill.NAs(a ~ ns(c, df = 3), sample.df)
    expect_equal(length(result), 5)
    expect_equal(colnames(result)[1], "a")

    ## right number of columns if 2 of the same variable used
    imputed.fmla <- fill.NAs(a ~ log(c) + sqrt(c), data = sample.df)
    expect_equal(dim(imputed.fmla)[2],  4)
  }

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

test_that("response variables with complex names", {
  data(nuclearplants)
  nuclearplants$cost[1] <- NA
  nuclearplants$cap[2] <- NA
  m <- lm(cost + t1 ~ cap + pr, data=nuclearplants)
  d <- model.frame(m, na.action=na.pass)
  # Name of response in this model is now `cost + t1`
  # Renaming the column to ensure special characters aren't
  # causing problems.
  d1 <- d
  names(d1)[1] <- "costplust1"

  expect_true(all(fill.NAs(d, all.covs=TRUE) == fill.NAs(d1, all.covs=TRUE)))


  # Addressing issue #100
  m2 <- lm(cbind(cost, t1) ~ cap + pr, data=nuclearplants)
  d2 <- model.frame(m2, na.action=na.pass)

  ## Disabling for now. See issue 104 for details on probable solution
  ## d3 <- d2
  ## names(d3)[1] <- "cbind"
  ## expect_true(all(fill.NAs(d2, all.covs=TRUE), fill.NAS(d3, all.covs=TRUE)))

})

test_that("strata() function handling", {

  set.seed(20150624)

  data.full <- data.frame(z = c(rep(1, 10), rep(0, 10)),
                          x = rnorm(20),
                          s = sample(c("A", "B", "C"), size = 20, replace = TRUE),
                          t = sample(c("UP", "DOWN"), size = 20, replace = TRUE))
  data.full$x[c(1, 2, 11)] <- NA

  # basic strata handling without NAs
  res1 <- fill.NAs(z ~ x + strata(s), data = data.full)
  expect_equal(dim(res1), c(20, 4)) # do not expand strata variable
  expect_false(any(is.na(res1)))

  res2 <- fill.NAs(z ~ x + strata(s) + strata(t), data = data.full)
  expect_equal(dim(res2), c(20, 5))
  expect_false(any(is.na(res2)))

  # now, let's knock out some strata levels
  data.NAs <- data.full
  data.NAs$s[sample(1:20, size = 3)] <- NA

  res3 <- fill.NAs(z ~ x + strata(s), data = data.NAs)
  expect_equal(sum(is.na(res3$s)), 3)
  # The following line should not error per #103
  glm(z ~ x + x.NA + strata(s), family=binomial, data=res3)

  res4 <- fill.NAs(z ~ x + strata(s, na.group = TRUE), data = data.NAs)
  expect_false(any(is.na(res4$s)))
  # Again, this should not error per #103
  glm(z ~ x + x.NA + strata(s, na.group=TRUE), family=binomial, data=res4)

  ## checking for terms attribute on the returned data.frame
  tt <- terms(res1)
  expect_false(is.null(attr(tt, "specials")$strata)) # the strata term is marked as such

  # if we spell things out, we should get a model on the imputed values
  xx <- glm(z ~ x + x.NA + strata(s), data = res1, family = binomial)
  expect_true(all(names(coef(xx))[1:3] %in% c("(Intercept)", "x", "x.NATRUE")))
  ## 4th coeff sometimes turns up as "strata(s)B", other times as "strata(s)s=s=B"
  ## the latter is less desirable, but no time to explore circumventing it.
  ## (Seen with: survival_2.37-7 ; R 3.1.2; x86_64-apple-darwin12.6.0, 64-bit.)
  expect_true(grep(glob2rx("strata*B"), names(coef(xx)))==4)
  expect_true(grep(glob2rx("strata*C"), names(coef(xx)))==5)

  # does not work yet:
  # yy <- glm(res1)
  # expect_equivalent(xx, yy)

  ## imputation should be per stratum
  expect_false(all(res1$x == (fill.NAs(z ~ x, data = data.full))$x))
})

test_that("Checking for fix to factors in fill.nas, mentioned in #103", {

  data(nuclearplants)
  nuclearplants$t1[1] <- NA
  f <- fill.NAs(pr ~ cap + factor(t1), data=nuclearplants)
  glm(f, family=binomial)

})
