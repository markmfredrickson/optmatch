#############################################################################
# Tests that should not be pushed to CRAN
#############################################################################

context("Not for CRAN")

test_that("scores with bigglm", {
  # Currently, scores() for bigglm is not supported
  if (require(biglm)) {
    n <- 16
    test.data <- data.frame(Z = rep(0:1, each = n/2),
                            X1 = rnorm(n, mean = 5),
                            X2 = rnorm(n, mean = -2, sd = 2))

    bgps <- bigglm(Z ~ X1 + X2, data = test.data, family = binomial())
    m1 <- lm(Z ~ scores(bgps), data=test.data)
    m2 <- lm(Z ~ predict(bgps, newdata=test.data), data=test.data)

    expect_equal(fitted(m1), fitted(m2))

    test.data2 <- test.data
    test.data2$X1[1] <- NA

    # predict's na.action=na.pass is ignored in bigglm, so we'll need to impute
    # beforehand to get the same results
    bgps2 <- bigglm(Z ~ X1 + X2, data = test.data2, family = binomial())
    expect_warning(m3 <- lm(Z ~ scores(bgps2), data=test.data2),
                   "Imputation and refitting of bigglm objects")
    m4 <- lm(Z ~ predict(bgps2, newdata=fill.NAs(test.data2)), data=test.data2)

    expect_equal(fitted(m3), fitted(m4))

    suppressWarnings(bgps3 <- bigglm(Z ~ X1*X2, data = test.data2, family = binomial()))
    fill.test.data <- fill.NAs(Z ~ X1*X2, data=test.data2)
    expect_warning(m5 <- lm(Z ~ scores(bgps3), data=fill.test.data),
                   "Imputation and refitting of bigglm objects")
    m6 <- lm(Z ~ predict(bgps3, newdata=fill.test.data), data=test.data2)

    expect_equal(fitted(m5), fitted(m6))
  }

})

test_that("match_on with bigglm distances", {
  if (require(biglm)) {
    n <- 16
    test.data <- data.frame(Z = rep(0:1, each = n/2),
                            X1 = rnorm(n, mean = 5),
                            X2 = rnorm(n, mean = -2, sd = 2),
                            B = rep(c(0,1), times = n/2))


    bgps <- bigglm(Z ~ X1 + X2, data = test.data, family = binomial())
    res.bg <- match_on(bgps, data = test.data)

    # compare to glm
    res.glm <- match_on(glm(Z ~ X1 + X2, data = test.data, family = binomial()))
    expect_equivalent(res.bg, res.glm)
  }
})

### next tests disabled due to an odd scoping-related error
### occurring only within the test, not in interactive use.
### (The problem is that `test.design` can't be found when
### the model is re-fit inside of `scores()`.)
###test_that("scores with svyglm (survey package) objects",
if (FALSE) {
  if (require(survey)) {
    n <- 16
    test.data <- data.frame(Z = rep(0:1, each = n/2),
                            X1 = rnorm(n, mean = 5),
                            X2 = rnorm(n, mean = -2, sd = 2),
                            B = rep(c(0,1), times = n/2))
    test.design <- svydesign(~1, probs=1, data=test.data)
    sglm <- svyglm(Z ~ X1 + X2, test.design, family = quasibinomial())
    expect_silent(scores(sglm, newdata=sglm$data))
}
}
###)
###test_that("match_on with svyglm (survey package) objects",
if (FALSE) {
  if (require(survey)) {
    n <- 16
    test.data <- data.frame(Z = rep(0:1, each = n/2),
                            X1 = rnorm(n, mean = 5),
                            X2 = rnorm(n, mean = -2, sd = 2),
                            B = rep(c(0,1), times = n/2))
    test.design <- svydesign(~1, probs=1, data=test.data)
    sglm <- svyglm(Z ~ X1 + X2, test.design, family = binomial())
    expect_silent(res.svy0 <- match_on(sglm, data=test.data, standardization.scale=1))
    expect_silent(res.svy1 <- match_on(sglm, data=test.data, standardization.scale=svy_sd))
    expect_silent(res.svy2 <- match_on(sglm, data=test.data, standardization.scale=svy_mad))
    ## Comparisons to glm: currently failing, disabled pending investigation.
    ## First step: figure out whether sglm/aglm are returning the same
    ## linear predictors.
    if (FALSE)
    {
        aglm <- glm(Z ~ X1 + X2, test.data, family = binomial())
    res.glm0 <- match_on(aglm, data=test.data, standardization.scale=1)
        res.glm1 <- match_on(aglm, data=test.data, standardization.scale=stats::sd)
        mad_lo  <- function(x) {
            ctr  <- quantile(x, probs=0.5, type=1)
            stats::mad(x, center=ctr, low=TRUE)
                       }
    res.glm2 <- match_on(aglm, data=test.data, standardization.scale=mad_lo)

    expect_equivalent(res.svy0, res.glm0)
    expect_equivalent(res.svy1, res.glm1)
    expect_equivalent(res.svy2, res.glm2)
    }
}
}
###)

test_that("cox model testing", {
  # this is the test case from:
  # https://github.com/markmfredrickson/optmatch/issues/44

  if (require(survival)) {
    coxps <- predict(coxph(Surv(start, stop, event) ~ age + year + transplant + cluster(id), data=heart))
    names(coxps) <- row.names(heart)
    coxmoA <- match_on(coxps, z = heart$event, caliper = 1)
    expect_true(max(coxmoA) <= 1)

    coxmoC <- match_on(coxps, within = exactMatch(event ~ transplant, data = heart), z = heart$event, caliper = 1)
    expect_true(max(coxmoC) <= 1)
  }
})
