################################################################################
# Match: methods to create distance specifications, possibliy sparse
################################################################################

context("match_on function")

test_that("Distances from glms", {

  n <- 16
  Z <- numeric(n)
  Z[sample.int(n, n/2)] <- 1
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- rep(c(0,1), n/2)

  test.glm <- glm(Z ~ X1 + X2 + B, family = binomial()) # the coefs should be zero or so

  result.glm <- match_on(test.glm)

  expect_true(validDistanceSpecification(result.glm))
  expect_equal(length(result.glm), (n/2)^2)

  # what about classes that inherit from glm?
  class(test.glm) <- c("foo", class(test.glm))

  result.foo <- match_on(test.glm)

  expect_true(validDistanceSpecification(result.foo))
  expect_equal(length(result.foo), (n/2)^2)

  # test that we can use GLMs with an attached data.frame
  df <- data.frame(Z. = Z, X1. = X1, X2. = X2, B. = B)
  expect_error(glm(Z. ~ X1. + X2. + B., family = binomial()))

  model.attach <- with(df, glm(Z. ~ X1. + X2. + B., family = binomial()))
  res.attach <- match_on(model.attach)

  expect_true(isTRUE(all.equal(res.attach, result.foo,
                               check.attributes = FALSE)))

})

test_that("Missingness in treatment", {
  set.seed(548243)
  Z <- sample(0:1, 10, TRUE)
  X <- rnorm(10)

  g <- glm(Z ~ X, family=binomial)

  expect_equal(sum(dim(match_on(g))), 10)

  Z.na <- Z; X.na <- X

  Z.na[1] <- NA
  X.na[2] <- NA

  g2 <- glm(Z ~ X.na, family=binomial)
  # Should have full recovery.
  expect_equal(sum(dim(match_on(g2))), 10)

  g3 <- glm(Z.na ~ X, family=binomial)
  # With missing treatment, drop observation.
  expect_equal(sum(dim(match_on(g3))), 9)

  g4 <- glm(Z.na ~ X.na, family=binomial)
  expect_equal(sum(dim(match_on(g4))), 9)

  d1 <- data.frame(Z,X.na)
  d2 <- data.frame(Z.na, X)
  d3 <- data.frame(Z.na, X.na)
  rm("Z","X","Z.na","X.na")

  # Passing data should have same recovery structure as before
  h1 <- glm(Z ~ X.na, family=binomial, data=d1)
  expect_equal(sum(dim(match_on(h1, data=d1))), 10)

  h2 <- glm(Z.na ~ X, family=binomial, data=d2)
  expect_equal(sum(dim(match_on(h2, data=d2))), 9)

  h3 <- glm(Z.na ~ X.na, family=binomial, data=d3)
  expect_equal(sum(dim(match_on(h3, data=d3))), 9)
})

test_that("Distances from formulas", {
  n <- 16
  Z <- rep(c(0,1), n/2)
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- as.factor(rep(c(0,1), each = n/2))

  test.data <- data.frame(Z, X1, X2, B)

  result.fmla <- match_on(Z ~ X1 + X2 + B, data = test.data)
  expect_true(validDistanceSpecification(result.fmla))

  # test pulling from the environment, like lm does
  result.envir <- match_on(Z ~ X1 + X2 + B)
  expect_equivalent(result.fmla, result.envir)

  expect_error(match_on(~ X1 + X2, data = test.data))
  expect_error(match_on(Z ~ 1, data = test.data))

  # checking diferent classes of responses
  res.one <- match_on(Z ~ X1)
  res.logical <- match_on(as.logical(Z) ~ X1)
  expect_equivalent(res.one, res.logical)

  tol <- 10^(9L - getOption("digits")) * sqrt(.Machine$double.eps)

  # euclidean distances
  # first, compute what the distances should be for the data.
  euclid <- as.matrix(dist(test.data[,-1], method = "euclidean", upper = T))
  z <- as.logical(Z)
  euclid <- euclid[z, !z]
  expect_true(all(abs(match_on(Z ~ X1 + X2 + B, method = "euclidean") - euclid) <
                      tol)) # there is some rounding error, but it is small

  # factor-related
  f0 <- as.factor(rep(1:4, each=n/4))
  f1 <- as.factor(rep(rep(1:2, each=2),n/4))
  f2 <- as.factor(rep(3:4, each=n/2))

  # Euclidean distances on a single factor should be 1 or 0
  tmp <- match_on(Z~f1, method="euclidean")
  expect_true(all(abs(tmp) < tol | abs(tmp - 1) < tol))

  euclid2 <- match_on(Z~f0, method="euclidean")
  expect_true(all(abs(euclid2) < tol | abs(euclid2 - 1) < tol))

  # with 2 orthogonal factors, distances should be 0, 1 or 2
  tmp <- match_on(Z ~ f1 + f2, method="euclidean")
  expect_true(all(abs(tmp) < tol | abs(tmp - 1) < tol | abs(tmp - sqrt(2)) < tol))

  # with a single 2 level factor, numeric version to give same distances
  expect_equal(as.matrix(match_on(Z~as.numeric(f1), method="euclidean")), as.matrix(match_on(Z~f1, method="euclidean")))

  # passing a function name for method
  #  expect_true(all(abs(match_on(Z ~ X1 + X2 + B, method = optmatch:::compute_euclidean) - euclid) < tol)) # there is some rounding error, but it is small
  # removed - no longer support user-functions in method

  # Mahalanobis distances involving factors

  expect_equal(as.matrix(match_on(Z~as.numeric(f1), method="mahalanobis")), as.matrix(match_on(Z~f1, method="mahalanobis")))

  # excluding matches combined with a formula
  stratify <- exactMatch(Z ~ B)
  res.strat <- match_on(Z ~ X1 + X2, within = stratify)
  expect_is(res.strat, "InfinitySparseMatrix")
  expect_equal(length(res.strat), 2 * (n/4)^2)

})

test_that("Issue 87: NA's in data => unmatchable, but retained, units in distances", {
  d <- data.frame(z  = c(1,1,1,0,0),
                  x1 = c(7, 9, NA, -1, 4),
                  x2 = c(1, 0, 0, 0, 1))

  rownames(d) <- c("A", "B", "C", "y", "z")
  g <- function(x) {
    tmp <- rep("OTHER", length(g))
    tmp[is.finite(x)]  <- "FINITE"
    tmp[!is.finite(x)] <- "INF"
    tmp[is.na(x)] <- "NA"
    return(tmp)
  }

  f <- function(method) {
    v <- as.matrix(match_on(z ~ x1 + x2, data = d, method = method))
    g(v)
  }

  expectedM <- c("FINITE", "FINITE", "INF", "FINITE", "FINITE", "INF")

  expect_equivalent(f("mahalanobis"), expectedM)
  expect_equivalent(f("euclid"), expectedM)
  expect_equivalent(f("rank_mahal"), expectedM)

  cal1 <- caliper(match_on(z~x1, data=d), width=1e3)
  expect_equivalent(g(as.matrix(match_on(z ~ x1 + x2, data = d,
                                         within=cal1, method = "mahalanobis"))),
                    expectedM)

  ## now with numeric method:
  zz <- d$z
  x1 <- d$x1
  names(zz) <- names(x1) <- rownames(d)

  v <- as.matrix(match_on(x1, z = zz))
  expect_equivalent(g(v), expectedM)

  ## glm should have the opposite behavior: automatically imputing
  expect_warning(v <- as.matrix(match_on(glm(z ~ x1 + x2, data = d, family = binomial))), "fitted")
  expect_equivalent(g(v), rep("FINITE", 6))
})

# while the formula method often handles mahalanobis distances, separating the tests for clarity
test_that("Mahalanobis distance calculations", {
  badData <- data.frame(Z = rep(c(0,1), 10),
                        all1 = as.factor(rep(1,20)),
                        badf1 = c(rep(1,3), rep(0,7)),
                        badf2 = c(rep(0,3), rep(1,7)))

  expect_error(match_on(Z ~ all1, data = badData), "contrasts can be applied only to factors with 2 or more levels")

  # even though the supplied data is a bad idea, it should work using the svd() decomposition
  res <- match_on(Z ~ badf1 + badf2, data = badData)
})

test_that("Distances from functions", {
  n <- 16
  Z <- c(rep(0, n/2), rep(1, n/2))
  X1 <- rep(c(1,2,3,4), each = n/4)
  B <- rep(c(0,1), n/2)
  test.data <- data.frame(Z, X1, B)

  sdiffs <- function(index, data, z) {
    abs(data[index[,1], "X1"] - data[index[,2], "X1"])
  }

  result.function <- match_on(sdiffs, z = Z, data = test.data)
  expect_equal(dim(result.function), c(8,8))

  # the result is a blocked matrix:
  # | 2 1 |
  # | 3 2 |

  expect_equal(mean(result.function), 2)

  # no treatment indicator
  expect_error(match_on(sdiffs, data = test.data))

  # no data
  expect_error(match_on(sdiffs, z = Z))

})

###### Using mad() instead of sd() for GLM distances
###
###result <- match_on(glm(pr ~ t1 + t2 + cost, data = nuclearplants, family = binomial()))
###
#### this is an odd test, but a simple way to make sure mad is running, not SD().
#### I would like a better test of the actual values, but it works
###test(mean(result$m) > 2)
###
###

test_that("Errors for numeric vectors", {
  expect_error(match_on(1:10))
})

###### Stratifying by a pipe (|) character in formulas
###
###main.fmla <- pr ~ t1 + t2
###strat.fmla <- ~ pt
###combined.fmla <- pr ~ t1 + t2 | pt
###
###result.main <- match_on(main.fmla, structure.fmla = strat.fmla, data = nuclearplants)
###result.combined <- match_on(combined.fmla, data = nuclearplants)
###
###test(identical(result.main, result.combined))
###
###### Informatively insist that one of formulas specify the treatment group
###shouldError(match_on(~t1+t2, structure.fmla=~pt, data=nuclearplants))
###test(identical(match_on(pr~t1+t2, structure.fmla=~pt, data=nuclearplants),
###               match_on(~t1+t2, structure.fmla=pr~pt, data=nuclearplants))
###     )
###### Finding "data" when it isn't given as an argument
###### Caveats:
###### * data's row.names get lost when you don't pass data as explicit argument;
###### thus testing with 'all.equal(unlist(<...>),unlist(<...>))' rather than 'identical(<...>,<...>)'.
###### * with(nuclearplants, match_on(fmla)) bombs for obscure scoping-related reasons,
###### namely that the environment of fmla is the globalenv rather than that created by 'with'.
###### This despite the facts that identical(fmla,pr ~ t1 + t2 + pt) is TRUE and that
###### with(nuclearplants, match_on(pr ~ t1 + t2 + pt)) runs fine.
###### But then with(nuclearplants, lm(fmla)) bombs too, for same reason, so don't worry be happy.
###attach(nuclearplants)
###test(all.equal(unlist(result.fmla),unlist(match_on(fmla))))
###test(all.equal(unlist(result.main),unlist(match_on(main.fmla, structure.fmla=strat.fmla))))
###test(all.equal(unlist(result.combined),unlist(match_on(combined.fmla))) )
###detach("nuclearplants")
###test(identical(fmla,pr ~ t1 + t2 + pt))
###test(all.equal(unlist(result.fmla),unlist(with(nuclearplants, match_on(pr ~ t1 + t2 + pt)))))
###test(identical(combined.fmla, pr ~ t1 + t2 | pt))
###test(all.equal(unlist(result.combined), unlist(with(nuclearplants, match_on(pr ~ t1 + t2 | pt)))))
###test(all.equal(unlist(result.fmla), unlist(with(nuclearplants[-which(names(nuclearplants)=="pt")],
###                                           match_on(update(pr ~ t1 + t2 + pt,.~.-pt + nuclearplants$pt))
###                                           )
###                                      )
###          )
###     )
###test(all.equal(unlist(result.combined), unlist(with(nuclearplants, match_on(pr ~ t1 + t2, structure.fmla=strat.fmla)))))

# test_that("Bigglm distances", {
# Moved to test.notforCRAN.R

test_that("Numeric: simple differences of scores", {
  # note: the propensity score method depends on this method as well, so if
  # those tests start failing, check here.

  scores <- rep(7, 12)
  z <- rep(c(0,1), 6)
  names(z) <- names(scores) <- letters[1:12]

  expect_true(all(match_on(scores, z = z) == 0))

  expect_true(all(match_on(z * 2, z = z) == 2))

  expect_true(all(match_on(z * -2, z = z) == 2))

  # proper errors
  expect_error(match_on(scores), "treatment")
  expect_error(match_on(scores, z = c(1,2)), "length")
  expect_error(match_on(c(1,2,3,4), z = c(0,1,0,1)), "names")

  # pass a caliper width, limits the comparisons that are going to be made.
  # the scores are going to be computed using abs diff, we can use the tools
  # developed in feasible.R

  scores <- rep(1:3, each = 4)
  names(scores) <- letters[1:12]

  # first, test the helper function scoreCaliper that generates the list of allowed comparisons.
  scres <- scoreCaliper(scores, z, 1)
  # match_on(scores, z) shows that there are 8 comparisons of size 2
  expect_equal(length(scres), 28)
  # every entry should be zero, the match_on function will compute the actual differences
  expect_equal(sum(scres), 0)
  # repeating above using non-integer values
  expect_equal(length(scoreCaliper(scores/3, z, 1/3)), 28)

  # try it as a matrix
  expect_equivalent(as.matrix(scres), matrix(c(0,0,0,0,Inf,Inf,
                                               0,0,0,0,Inf,Inf,
                                               rep(0, 6),
                                               rep(0, 6),
                                               Inf, Inf, 0,0,0,0,
                                               Inf,Inf,0,0,0,0), nrow = 6, ncol = 6))

  # repeat with match_on
  expect_equal(length(match_on(scores, z = z, caliper = 1)), 28) # 6 * 6 - 8
  expect_equal(length(match_on(scores, z = z, caliper = 1.5)), 28)

  # combine the caliper width with an within argument
  b <- rep(1:3, 4)
  ez <- exactMatch(z ~ b)

  res <- match_on(scores, z = z, caliper = 1, within = ez)
  expect_equal(length(res), 9)

  # next test includes treated and control units that are excluded entirely
  # with caliper = 1
  scores2 <- c(scores, -50, 100, 200, -100)
  z2 <- c(z, 1,0,1,0)
  names(scores2) <- names(z2) <- letters[1:(length(scores2))]
  res <- match_on(scores2, z = z2, caliper = 1)
  expect_equal(length(res), 28) # effectively same result as without the new units

  # caliper must be of length 1
  expect_error(match_on(scores2, z = z2, caliper = c(1,2)), "scalar")
})

test_that("Numeric, issues with vector names", {
  # #189

  # If z is not named, it should copy over x's names
  z <- c(0, 0, 1, 1)
  x <- c("a" = NA, "b" = 1, "c" = 2, "d" = 3)
  m1 <- match_on(x, z = z)
  expect_equal(dim(m1), c(2, 2))
  expect_equal(m1, match_on(z ~ x, method = "euclidean",
                            data = data.frame(x,z)),
               check.attributes = FALSE)

  # If z is named, and the same as x, everything should work.
  z <- c("a" = 0,  "b" = 0, "c" = 1, "d" = 1)
  x <- c("a" = NA, "b" = 1, "c" = 2, "d" = 3)
  m2 <- match_on(x, z = z)
  expect_equal(dim(m2), c(2, 2))
  expect_equal(m2, match_on(z ~ x, method = "euclidean",
                        data = data.frame(x,z)),
               check.attributes = FALSE)
  expect_identical(m1, m2)

  # if z is named same as x, but misordered, order and continue on.
  z <- c("b" = 0,  "a" = 0, "d" = 1, "c" = 1)
  x <- c("a" = NA, "a" = 1, "c" = 2, "d" =  3)
  m3 <- match_on(x, z = z)
  expect_equal(dim(m3), c(2, 2))
  expect_equal(m3, match_on(z ~ x, method = "euclidean",
                        data = data.frame(x,z)),
               check.attributes = FALSE)

  # if names of z and x are not the same, error
  z <- c("a" = 0,  "c" = 0, "b" = 1)
  x <- c("a" = NA, "b" = 1, "d" = 2)
  expect_error(match_on(x, z = z), "names")


})

test_that("matrix, ISM methods' within args",{
  scores <- rep(1:3, each = 4)
  z <- rep(c(0,1), 6)
  b <- rep(1:3, 4)
  names(z) <- names(scores) <- letters[1:12]
  ez <- exactMatch(z ~ b)
  sISM  <- match_on(scores, z=z)
  expect_s4_class(sISM, "DenseMatrix")
  expect_equivalent(match_on(sISM, within=ez), sISM+ez)
  sISM2  <- as.InfinitySparseMatrix(sISM)
  expect_equivalent(match_on(sISM2, within=ez), sISM+ez)
})

test_that("use of matrix, ISM, BISM methods w/ caliper arg", {
  scores <- rep(1:3, each = 4)
  z <- rep(c(0,1), 6)
  b <- rep(1:3, 4)
  names(z) <- names(scores) <- letters[1:12]
  sISM  <- match_on(scores, z=z)
  expect_s4_class(sISM, "DenseMatrix")
  sISM_cal  <- match_on(scores, z=z, caliper=1)
  expect_equivalent(sISM_cal, match_on(sISM, caliper=1))
  expect_equivalent(sISM_cal,
                  match_on(as.InfinitySparseMatrix(sISM), caliper=1)
                  )

  ez <- exactMatch(z ~ b)
  res <- match_on(scores, z = z, # checked in "Numeric: simple..."
                  caliper = 1, within = ez) # test above
  expect_equivalent(res, sISM_cal + ez)

  res1  <- match_on(sISM, caliper=1, within=ez)
  expect_equivalent(res, res1)
})

test_that("update() of match_on created objects", {
  Z <- rep(c(1,0), 10)
  X1 <- rep(c(1,3), 10)
  X2 <- rep(c(1,5), 10)
  B <- rep(c(1,2), each = 10)

  names(Z) <- names(X1) <- names(X2) <- names(B) <- letters[1:20]

  # first, with dense problems
  basic <- match_on(X1, z = Z)
  use.x2 <- update(basic, x = X2)
  expect_true(all(basic == 2))
  expect_true(all(use.x2 == 4))

  X1 <- X2
  expect_true(all((update(basic) - use.x2) == 0))

  # next, with sparse problems
  X1 <- rep(c(1,3), 10)
  names(X1) <- letters[1:20]

  basic <- match_on(X1, z = Z, within = exactMatch(Z ~ B))
  use.x2 <- update(basic, x = X2)
  expect_true(all(basic == 2))
  expect_true(all(use.x2 == 4))

  X1 <- X2
  expect_true(all((update(basic) - use.x2) == 0))

  # now create a distspec with addition, update should error
  X1 <- rep(c(1,3), 10)
  names(X1) <- letters[1:20]

  stratified <- match_on(X1, z = Z, within = exactMatch(Z ~ B))
  unstratified <- match_on(X1, z = Z)
  em <- exactMatch(Z ~ B)

  expect_equivalent(stratified, unstratified + em)

  expect_error(update(unstratified + em, x = X2))
})

test_that("Issue #44", {
  # problem: `within` negates proper caliper

  scores <- rep(c(1,2,2,3), each = 25)
  z <- rep(c(0,1), each = 50)

  names(scores) <- paste("A", 1:100, sep = "")

  # get the caliper only results
  res.cal <- match_on(scores, z = z, caliper = 1)
  expect_equal(max(res.cal), 1)

  # now make up a within arg
  names(z) <- names(scores)
  www  <- exactMatch(x = as.factor(rep(c(0,1), 50)),
                     treatment = z)

  res.w  <- match_on(scores, z = z, within = www)
  expect_true(max(res.w) > 1)

  # they should safely interact
  res.cal.w <- match_on(scores, z = z, within = www, caliper = 1)
  expect_equal(max(res.cal.w), 1)


})

test_that("Issue 48: caliper is a universal argument", {

  Z <- rep(c(1,0), 5)
  X <- -4:5
  names(X) <- names(Z) <- letters[1:10]

  res.num <- match_on(X, z = Z, caliper = 1)
  expect_true(all(res.num <= 1))

  res.glm <- match_on(glm(Z ~ X, family = binomial), caliper = 1)
  expect_true(all(res.glm <= 1))

  res.fmla <- match_on(Z ~ X, caliper = 1)
  expect_true(all(res.fmla <= 1))
})

# test_that("bayesglm, brglm", {
#
#   # packaages are added at top of file
#   data(nuclearplants)
#
#   by <- bayesglm(pr ~ cost, data=nuclearplants, family=binomial)
#   expect_true(all(class(by) == c("bayesglm", "glm", "lm")))
#   m1 <- match_on(by, data=nuclearplants)
#   expect_true(class(m1)[1] %in% c("InfinitySparseMatrix", "BlockedInfinitySparseMatrix", "DenseMatrix"))
#
#   br <- brglm(pr ~ cost, data=nuclearplants, family=binomial, method="glm.fit")
#   expect_true(all(class(br) == c("brglm", "glm", "lm")))
#   m2 <- match_on(br, data=nuclearplants)
#   expect_true(class(m2)[1] %in% c("InfinitySparseMatrix", "BlockedInfinitySparseMatrix", "DenseMatrix"))
# })

test_that("standardization scale from within match_on", {
  n <- 16
  Z <- numeric(n)
  Z[sample.int(n, n/2)] <- 1
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- rep(c(0,1), n/2)

  test.glm <- glm(Z ~ X1 + X2 + B, family = binomial()) # the coefs should be zero or so

  expect_silent(result.glm0 <- match_on(test.glm))
  expect_is(result.glm0, "matrix")

  expect_silent(result.glm1 <- match_on(test.glm, standardization.scale=mad))
  expect_equivalent(result.glm1, result.glm0)

  expect_silent(result.glm2 <- match_on(test.glm, standardization.scale=sd))
  expect_is(result.glm2, "matrix")

  expect_silent(result.glm4 <- match_on(test.glm, standardization.scale=1))
  expect_is(result.glm4, "matrix")
})

test_that("standardization_scale with svyglm",{
  if (requireNamespace("survey", quietly = TRUE)) {
    data(api, package = "survey")
    d <- survey::svydesign(id = ~ 1, probs = 1, data = apistrat)
    sglm <- survey::svyglm(sch.wide ~ ell + meals + mobility, design = d,
                           family=quasibinomial())
    aglm <- glm(sch.wide ~ ell + meals + mobility,data = apistrat,
                family = binomial)
    expect_equivalent(sglm$linear.predictors, aglm$linear.predictors)
    sd_sglm <- standardization_scale(sglm$linear.predictor,
                                     trtgrp=sglm$y,
                                     standardizer = svy_sd,
                                     svydesign_ = sglm$'survey.design')
    sd_aglm <- standardization_scale(aglm$linear.predictor,
                                     trtgrp=aglm$y,
                                     standardizer = stats::sd)
    expect_equivalent(sd_sglm, sd_aglm)

    d_w <- survey::svydesign(id = ~ 1, weights = ~ pw, data = apistrat)
    sglm_w <- survey::svyglm(sch.wide ~ ell + meals + mobility, design = d_w,
                             family = quasibinomial())
    mad_sglm_w <- standardization_scale(sglm_w$linear.predictor,
                                        trtgrp=sglm_w$y,
                                        standardizer = svy_mad,
                                        svydesign_ = sglm_w$'survey.design')
    expect_true(is.numeric(mad_sglm_w))
    expect_gt(mad_sglm_w, 0)
    sd_sglm_w <- standardization_scale(sglm_w$linear.predictor,
                                     trtgrp=sglm_w$y,
                                     standardizer = svy_sd,
                                     svydesign_ = sglm_w$'survey.design')
    expect_true(is.numeric(sd_sglm_w))
    expect_gt(sd_sglm_w, 0)
  }
  expect_true(TRUE) # avoiding empty test warning
})

test_that("Building exactMatch from formula with strata", {

  d <- data.frame(x = rnorm(8),
                  t = rep(0:1, 4),
                  z = rep(0:1, each=4))

  em <- exactMatch(t ~ z, data=d)

  nw <- makeWithinFromStrata(t ~ strata(z), d)

  expect_true(is(nw, "list"))
  expect_true(length(nw) == 2)
  expect_true(all(names(nw) == c("x", "within")))

  expect_true(is(nw$x, "formula"))
  expect_true(is(nw$within, "BlockedInfinitySparseMatrix"))
  expect_true(all(unlist(subdim(nw$within)) == 2))

  expect_equivalent(em, nw$within)
  expect_equivalent(t ~ 1, nw$x)

  nw1 <- makeWithinFromStrata(!!t ~ strata(z), d)
  expect_equivalent(nw1$within, nw$within)

  nw2 <- makeWithinFromStrata(t ~ x + strata(z), d)

  expect_equivalent(em, nw2$within)
  expect_equivalent(t ~ x, nw2$x)


})


test_that("Using strata instead of within arguments", {
  data(nuclearplants)

  m1 <- match_on(pr ~ cost, within=exactMatch(pr ~ pt, data=nuclearplants),
                  data=nuclearplants)
  m2 <- match_on(pr ~ cost + strata(pt), data=nuclearplants)
  m2b <- match_on(pr ~ cost, data=nuclearplants)

  expect_true(is(m1, "BlockedInfinitySparseMatrix"))

  expect_true(is(m1, "BlockedInfinitySparseMatrix"))

  expect_true(is(m2, "BlockedInfinitySparseMatrix"))
  expect_true(all.equal(m1, m2, check.attributes=FALSE))
  expect_true(!isTRUE(all.equal(m2, m2b, check.attributes=FALSE)))

  # handling more complicated strata calls
  m3 <- match_on(pr ~ cost, within=exactMatch(pr ~ pt + ct + ne,
                                               data=nuclearplants),
                  data=nuclearplants)
  m4 <- match_on(pr ~ cost + strata(pt) + strata(ct, ne), data=nuclearplants)

  expect_true(all.equal(m3, m4, check.attributes=FALSE))

  e1 <- exactMatch(pr ~ ne, data=nuclearplants)
  e2 <- exactMatch(pr ~ ne + ct, data=nuclearplants)

  m5 <- match_on(pr ~ cost, within=e2, data=nuclearplants)
  m6 <- match_on(pr ~ cost + strata(ct), within=e1, data=nuclearplants)
  m7 <- match_on(pr ~ cost + strata(ct, ne), data=nuclearplants)

  expect_true(all.equal(m5, m6, check.attributes=FALSE))
  expect_true(all.equal(m5, m7, check.attributes=FALSE))

})

test_that("strata in GLMs", {

  data(nuclearplants)

  m1 <- match_on(glm(pr ~ t1 + ne, data=nuclearplants, family=binomial),
                 within=exactMatch(pr ~ ne, data=nuclearplants),
                 data=nuclearplants)
  m2 <- match_on(glm(pr ~ t1 + strata(ne), data=nuclearplants, family=binomial),
                 data=nuclearplants)

  expect_true(is(m1, "BlockedInfinitySparseMatrix"))
  expect_true(is(m2, "BlockedInfinitySparseMatrix"))
  expect_true(all.equal(m1, m2, check.attributes=FALSE))

  m3 <- match_on(glm(pr ~ t1 + ne + interaction(ct,pt), data=nuclearplants,
                     family=binomial),
                 within=exactMatch(pr ~ ne + ct*pt, data=nuclearplants),
                 data=nuclearplants)
  m4 <- match_on(glm(pr ~ t1 + strata(ne) + strata(ct, pt),
                     data=nuclearplants, family=binomial), data=nuclearplants)

  expect_true(all.equal(m3, m4, check.attributes=FALSE))

  e1 <- exactMatch(pr ~ ne, data=nuclearplants)
  e2 <- exactMatch(pr ~ ne + ct, data=nuclearplants)

  m5 <- match_on(glm(pr ~ cost + ne + ct, data=nuclearplants, family=binomial),
                 within=e2, data=nuclearplants)
  m6 <- match_on(glm(pr ~ cost + ne + strata(ct), data=nuclearplants,
                     family=binomial),
                 within=e1, data=nuclearplants)
  m7 <- match_on(glm(pr ~ cost + strata(ct) + strata(ne), data=nuclearplants,
                     family=binomial),
                 data=nuclearplants)

  expect_true(all.equal(m5, m6, check.attributes=FALSE))
  expect_true(all.equal(m5, m7, check.attributes=FALSE))


  # strata(a,b) is equivalent to interaction(a,b)
  m8 <- match_on(glm(pr ~ cost + strata(ct,ne), data=nuclearplants,
                     family=binomial),
                 data=nuclearplants)
  suppressWarnings(
    m9 <- match_on(glm(pr ~ cost + interaction(ne, ct) + strata(ct),
                       data=nuclearplants, family=binomial),
                   within=e1, data=nuclearplants)
  )
  m10 <- match_on(glm(pr ~ cost + interaction(ne,ct), data=nuclearplants,
                      family=binomial),
                  within=e2, data=nuclearplants)
  m10b <- match_on(glm(pr ~ cost + ne*ct, data=nuclearplants, family=binomial),
                  within=e2, data=nuclearplants)
  # m9 is a bit weight because of the double inclusion of ct, and is an unlikely
  # way for users to enter code, but the extra ct is of course ignored.

  expect_true(all.equal(m8, m9, check.attributes=FALSE))
  expect_true(all.equal(m8, m10, check.attributes=FALSE))

  # Fixed effects should be added when exactMatching
  m11 <- match_on(glm(pr ~ cost, data=nuclearplants, family=binomial),
                  within=e2, data=nuclearplants)

  expect_false(isTRUE(all.equal(m8, m11, check.attributes=FALSE)))


})

test_that("Subsetting an ISM by passing a new data object to match_on", {

  data <- data.frame(z = rep(c(0,1), 13), x = rnorm(26))
  rownames(data) <- letters

  x <- match_on(z ~ x, data = data)
  expect_equal(dim(x), c(13, 13))

  d2 <- data.frame(w = 10:15)
  rownames(d2) <- letters[10:15]

  y <- match_on(x, data = d2)
  expect_equal(dim(y), c(3, 3))

  y2 <- match_on(optmatch:::as.InfinitySparseMatrix(x), data = d2)
  expect_equal(dim(y2), c(3, 3))
})

test_that("#114 informative error if caliper in formula", {

  data(nuclearplants)

  m <- match_on(pr ~ cost, data = nuclearplants)
  expect_error(match_on(pr ~ t1 + caliper(m), data = nuclearplants),
               "be applied via")
  expect_error(match_on(pr ~ t1 + caliper(m), data = nuclearplants),
               "within\\ =\\ caliper\\(m\\)")
  expect_error(match_on(pr ~ t1 + caliper(m) + caliper(n), data = nuclearplants),
               "within\\ =\\ caliper\\(m\\)\\ \\+\\ caliper\\(n\\)")

  em <- exactMatch(pr ~ pt, data = nuclearplants)
  expect_error(match_on(pr ~ t1 + caliper(m) + caliper(n),
                        data = nuclearplants, within = em),
               "within\\ =\\ em\\ \\+\\ caliper\\(m\\)\\ \\+\\ caliper\\(n\\)")

})

test_that("variables named strata or caliper are allowed", {

  data <- data.frame(z = rep(0:1, each = 10),
                     x = rnorm(20),
                     strata = rep(1:4, times = 5),
                     caliper = rnorm(20))
  # Try some valid albeit confusing inputs, along with some oddly formed input
  expect_silent(match_on(z ~ x + strata, data = data))
  expect_silent(match_on(z ~ x + caliper, data = data))
  expect_silent(match_on(z ~ x + strata + strata(strata), data = data))
  expect_silent(match_on(z ~ x + caliper + strata(strata), data = data))

  data(nuclearplants)
  m1 <- match_on(pr ~ t1 + t2 + strata(pt), data = nuclearplants)

  names(nuclearplants)[11] <- "strata"

  m2 <- match_on(pr ~ t1 + t2 + strata(strata), data = nuclearplants)
  expect_identical(as.matrix(m1), as.matrix(m2))

})

test_that("dot and strata in formula", {
  data <- data.frame(z = rep(0:1, each = 5),
                     x = rnorm(10),
                     s = rep(0:1, times = 5))

  expect_error(m <- match_on(z ~ . - s + strata(s), data = data),
               "Cannot use . expansion", fixed = TRUE)

})

test_that("#123: Supporting NA's in treatment, match_on.formula", {
  data <- data.frame(z = rep(0:1, each = 5),
                     x = rnorm(10))
  m <- match_on(z ~ x, data = data)
  expect_equal(dim(m), c(5, 5))
  expect_equal(rownames(m), rownames(data)[!is.na(data$z) & data$z == 1])
  expect_equal(colnames(m), rownames(data)[!is.na(data$z) & data$z == 0])

  # Now add an NA

  data$z[1] <- NA
  m <- match_on(z ~ x, data = data)
  expect_equal(dim(m), c(5, 4))
  expect_equal(rownames(m), rownames(data)[!is.na(data$z) & data$z == 1])
  expect_equal(colnames(m), rownames(data)[!is.na(data$z) & data$z == 0])

  data$z[c(2,5,6,7)] <- NA
  m <- match_on(z ~ x, data = data)
  expect_equal(dim(m), c(3, 2))
  expect_equal(rownames(m), rownames(data)[!is.na(data$z) & data$z == 1])
  expect_equal(colnames(m), rownames(data)[!is.na(data$z) & data$z == 0])


})

test_that("#123: Supporting NA's in treatment, match_on.numeric", {
  z <- rep(0:1, each = 5)
  x <- rnorm(10)
  names(z) <- names(x) <- 1:10
  m <- match_on(x, z = z)
  expect_equal(dim(m), c(5, 5))
  expect_equal(colnames(m), names(z[1:5])[!is.na(z[1:5])])
  expect_equal(rownames(m), names(z[6:10])[!is.na(z[6:10])])

  data <- data.frame(z, x)
  m2 <- match_on(x, z = z, data = data)
  expect_equivalent(m, m2)

  # Now add an NA

  z[1] <- NA
  m <- match_on(x, z = z)
  expect_equal(dim(m), c(5, 4))
  expect_equal(colnames(m), names(z[1:5])[!is.na(z[1:5])])
  expect_equal(rownames(m), names(z[6:10])[!is.na(z[6:10])])

  data <- data.frame(z, x)
  m2 <- match_on(x, z = z, data = data)
  expect_equivalent(m, m2)

  z[c(2,5,6,7)] <- NA
  m <- match_on(x, z = z)
  expect_equal(dim(m), c(3, 2))
  expect_equal(colnames(m), names(z[1:5])[!is.na(z[1:5])])
  expect_equal(rownames(m), names(z[6:10])[!is.na(z[6:10])])

  data <- data.frame(z, x)
  m2 <- match_on(x, z = z, data = data)
  expect_equivalent(m, m2)
})

test_that("#123: Supporting NA's in treatment, match_on.function", {

  data <- data.frame(z = rep(0:1, each = 5),
                     x = rnorm(10))

  sdiffs <- function(index, data, z) {
    abs(data[index[,1], "x"] - data[index[,2], "x"])
  }

  m <- match_on(sdiffs, z = data$z, data = data)
  expect_equal(dim(m), c(5, 5))
  expect_equal(colnames(m), rownames(data)[1:5][!is.na(data$z[1:5])])
  expect_equal(rownames(m), rownames(data)[6:10][!is.na(data$z[6:10])])

  data$z[1] <- NA

  m <- match_on(sdiffs, z = data$z, data = data)
  expect_equal(dim(m), c(5, 4))
  expect_equal(colnames(m), rownames(data)[1:5][!is.na(data$z[1:5])])
  expect_equal(rownames(m), rownames(data)[6:10][!is.na(data$z[6:10])])

  data$z[c(2,5,6,7)] <- NA

  m <- match_on(sdiffs, z = data$z, data = data)
  expect_equal(dim(m), c(3, 2))
  expect_equal(colnames(m), rownames(data)[1:5][!is.na(data$z[1:5])])
  expect_equal(rownames(m), rownames(data)[6:10][!is.na(data$z[6:10])])

})

test_that("#123: Supporting NA's in treatment, match_on.glm/bigglm", {

  data <- data.frame(z = rep(0:1, each = 5),
                     x = rnorm(10))

  mod <- glm(z ~ x, data = data, family = binomial)

  m <- match_on(mod)
  expect_equal(dim(m), c(5, 5))
  expect_equal(colnames(m), rownames(data)[1:5][!is.na(data$z[1:5])])
  expect_equal(rownames(m), rownames(data)[6:10][!is.na(data$z[6:10])])

  m2 <- match_on(mod, data = data)
  expect_equivalent(m, m2)

  data$z[1] <- NA

  mod <- glm(z ~ x, data = data, family = binomial)

  m <- match_on(mod)
  expect_equal(dim(m), c(5, 4))
  expect_equal(colnames(m), rownames(data)[1:5][!is.na(data$z[1:5])])
  expect_equal(rownames(m), rownames(data)[6:10][!is.na(data$z[6:10])])

  m2 <- match_on(mod, data = data)
  expect_equivalent(m, m2)

  data$z[c(2,5,6,7)] <- NA

  mod <- glm(z ~ x, data = data, family = binomial)

  m <- match_on(mod)
  expect_equal(dim(m), c(3, 2))
  expect_equal(colnames(m), rownames(data)[1:5][!is.na(data$z[1:5])])
  expect_equal(rownames(m), rownames(data)[6:10][!is.na(data$z[6:10])])

  m2 <- match_on(mod, data = data)
  expect_equivalent(m, m2)
})

test_that("147: within=caliper with NA's", {
  data(nuclearplants)
  #testing with missing data
  nuclearplants$cost[2] <- NA
  m.missing <- match_on(pr ~ cost, within=exactMatch(pr ~ pt, data=nuclearplants),
                        data=nuclearplants)
  expect_true(is(m.missing, "BlockedInfinitySparseMatrix"))
  m.drop <- match_on(pr ~ cost, within=exactMatch(pr ~ pt, data=nuclearplants[-2, ]),
                     data=nuclearplants[-2,])
  expect_true(all.equal(m.drop@.Data, m.missing@.Data))

  np <- nuclearplants[1:6,]
  np$cost[2] <- NA
  c <- caliper(match_on(pr ~ cost, data = np, method = "euclidean"), width = 100)
  m.1 <- match_on(pr ~ cap, within = c, data = np)

  m <- match_on(pr ~ cap, data = np)
  m.2 <- m + c
  expect_true(all.equal(dim(m.1),dim(m.2)))
  expect_true(is(m.1, "InfinitySparseMatrix"))
  expect_true(is(m.2, "InfinitySparseMatrix"))
  expect_true(all.equal(sort(m.2@.Data), sort(m.1@.Data)))
})

test_that("Exclude argument for match_on with caliper arg", {
    set.seed(10303920)
    n <- 20
    x <- rnorm(n, sd = 4)
    names(x) <- letters[1:20]
    z <- c(rep(0, 10), rep(1, 10))
    mint <- names(which.min(x[z == 1]))
    maxc <- names(which.min(x[z == 0]))
    cal <- sd(x) / 5

    m <- match_on(x, z = z, caliper = cal, exclude = c(mint, maxc))

    mm <- as.matrix(m)
    expect_true(sum(is.finite(mm)) < 100)
    expect_equal(sum(is.finite(mm[mint, ])), 10)
    expect_equal(sum(is.finite(mm[, maxc])), 10)

    ## exclude only treatment
    mt <- match_on(x, z = z, caliper = cal, exclude = mint)
    mm <- as.matrix(mt)
    expect_true(sum(is.finite(mm)) < 100)
    expect_equal(sum(is.finite(mm[mint, ])), 10)
    expect_true(sum(is.finite(mm[, maxc])) < 10)

    ## exclude only control
    mc <- match_on(x, z = z, caliper = cal, exclude = maxc)
    mm <- as.matrix(mc)
    expect_true(sum(is.finite(mm)) < 100)
    expect_equal(sum(is.finite(mm[, maxc])), 10)
    expect_true(sum(is.finite(mm[mint, ])) < 10)

    ## repeating with other versions of match_on
    df <- data.frame(z = z, x = x)
    mce <- match_on(z ~ x, data = df, caliper = cal, exclude = maxc, method = "euclidean")
    mme <- as.matrix(mce)
    expect_true(sum(is.finite(mme)) < 100)
    expect_equal(sum(is.finite(mme[, maxc])), 10)
    expect_true(sum(is.finite(mme[mint, ])) < 10)

    mg <- match_on(glm(z ~ x, data = df, family = binomial),
                    caliper = cal, exclude = c(mint, maxc))
    mm <- as.matrix(mg)
    expect_true(sum(is.finite(mm)) < 100)
    expect_equal(sum(is.finite(mm[mint, ])), 10)
    expect_equal(sum(is.finite(mm[, maxc])), 10)

    exmat <- matrix(1:9, nrow = 3, dimnames = list(c('a', 'b', 'c'), c('x', 'y', 'z')))
    mm <- as.matrix(match_on(exmat, caliper = 4, exclude = 'z'))
    expect_equivalent(mm[, 'y'], c(4, Inf, Inf))
    expect_equivalent(mm[, 'z'], 7:9)

    mm <- as.matrix(match_on(as.InfinitySparseMatrix(exmat), caliper = 4, exclude = 'z'))
    expect_equivalent(mm[, 'y'], c(4, Inf, Inf))
    expect_equivalent(mm[, 'z'], 7:9)
})


test_that("No longer support user-defined distances in match_on.formula", {
  data(nuclearplants)
  expect_warning(match_on(pr ~ cost, data = nuclearplants, method = optmatch:::compute_euclidean), "not supported")
})
