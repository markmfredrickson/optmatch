################################################################################
# Mdist: methods to create distance specifications, possibliy sparse
################################################################################

library(testthat)
library(biglm)
library(survival)

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

  res.factor <- match_on(as.factor(Z) ~ X1)
  expect_equivalent(res.one, res.factor)

  # euclidean distances
  # first, compute what the distances should be for the data.
  euclid <- as.matrix(dist(test.data[,-1], method = "euclidean", upper = T))
  z <- as.logical(Z)
  euclid <- euclid[z, !z]
  expect_true(all(abs(match_on(Z ~ X1 + X2 + B, method = "euclidean") - euclid) <
    .00001)) # there is some rounding error, but it is small

  # passing a function name for method
  expect_true(all(abs(match_on(Z ~ X1 + X2 + B, method = "compute_euclidean") - euclid) <
    .00001)) # there is some rounding error, but it is small


  # excluding matches combined with a formula
  stratify <- exactMatch(Z ~ B)
  res.strat <- match_on(Z ~ X1 + X2, within = stratify)
  expect_is(res.strat, "InfinitySparseMatrix")
  expect_equal(length(res.strat), 2 * (n/4)^2)

})

# while the formula method often handles mahalanobis distances, separating the tests for clarity
test_that("Mahalanobis distance calcualtions", {
  badData <- data.frame(Z = as.factor(rep(c(0,1), 10)),
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
test_that("Bigglm distances", {
    n <- 16
    Z <- c(rep(0, n/2), rep(1, n/2))
    X1 <- rnorm(n, mean = 5)
    X2 <- rnorm(n, mean = -2, sd = 2)
    B <- rep(c(0,1), n/2)
    test.data <- data.frame(Z, X1, X2, B)

    bgps <- bigglm(Z ~ X1 + X2, data = test.data, family = binomial())
    res.bg <- match_on(bgps, data = test.data)

    # compare to glm
    res.glm <- match_on(glm(Z ~ X1 + X2, data = test.data, family = binomial()))
    expect_equivalent(res.bg, res.glm)
})

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

  # this is the test case from:
  # https://github.com/markmfredrickson/optmatch/issues/44

  coxps <- predict(coxph(Surv(start, stop, event) ~ age + year + transplant + cluster(id), data=heart))
  names(coxps) <- row.names(heart)
  coxmoA <- match_on(coxps, z = heart$event, caliper = 1)
  expect_true(max(coxmoA) <= 1)

  coxmoC <- match_on(coxps, within = exactMatch(event ~ transplant, data = heart), z = heart$event, caliper = 1)
  expect_true(max(coxmoC) <= 1)

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

test_that("numeric standardization scale", {
  n <- 16
  Z <- numeric(n)
  Z[sample.int(n, n/2)] <- 1
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- rep(c(0,1), n/2)

  test.glm <- glm(Z ~ X1 + X2 + B, family = binomial()) # the coefs should be zero or so

  result.glm <- match_on(test.glm, standardization.scale=1)
})

