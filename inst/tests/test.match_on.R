################################################################################
# Mdist: methods to create distance specifications, possibliy sparse
################################################################################

library(testthat)

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
  
  expect_true(is(result.glm, "DistanceSpecification"))
  expect_equal(length(result.glm), (n/2)^2)
  
})

test_that("Distances from formulas", {
  n <- 16
  Z <- rep(c(0,1), n/2)
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- as.factor(rep(c(0,1), each = n/2))

  test.data <- data.frame(Z, X1, X2, B)

  result.fmla <- match_on(Z ~ X1 + X2 + B, data = test.data)
  expect_true(is(result.fmla, "DistanceSpecification"))

  # test pulling from the environment, like lm does
  result.envir <- match_on(Z ~ X1 + X2 + B)
  expect_equal(result.fmla, result.envir)

  expect_error(match_on(~ X1 + X2, data = test.data))
  expect_error(match_on(Z ~ 1, data = test.data))

  # checking diferent classes of responses
  res.one <- match_on(Z ~ X1) 
  res.logical <- match_on(as.logical(Z) ~ X1)
  expect_identical(res.one, res.logical)

  res.factor <- match_on(as.factor(Z) ~ X1)
  expect_identical(res.one, res.factor)

  # specifying distances
  euclid <- as.matrix(dist(test.data[,-1], method = "euclidean", upper = T))
  z <- as.logical(Z)
  euclid <- euclid[z, !z]
  # there are 3 columns x1, x2, b, so diag(3) is the identity matrix
  # also, match_on returns euclidean distance squared, so square the reference
  # result
  expect_true(all(abs(match_on(Z ~ X1 + X2 + B, inv.scale.matrix = diag(3)) - euclid^2) <
    .00001)) # there is some rounding error, but it is small


  # excluding matches combined with a formula
  stratify <- exactMatch(Z ~ B)
  res.strat <- match_on(Z ~ X1 + X2, exclusions = stratify)
  expect_is(res.strat, "InfinitySparseMatrix")
  expect_equal(length(res.strat), 2 * (n/4)^2)

})

test_that("Distances from functions", {
  n <- 16
  Z <- c(rep(0, n/2), rep(1, n/2))
  X1 <- rep(c(1,2,3,4), each = n/4)
  B <- rep(c(0,1), n/2)
  test.data <- data.frame(Z, X1, B)
  
  sdiffs <- function(t, c) {
    abs(t$X1 - c$X1)
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
  if (require('biglm')) {
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
    expect_equal(res.bg, res.glm) 
  }
})

test_that("Jake found a bug 2010-06-14", {
  ### Issue appears to be a missing row.names/class

  jb.sdiffs <- function(treatments, controls) {
    abs(treatments$X1 - controls$X2)
  }
  
  n <- 16
  Z <- c(rep(0, n/2), rep(1, n/2))
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- rep(c(0,1), n/2)
  test.data <- data.frame(Z, X1, X2, B)

  absdist1 <- match_on(jb.sdiffs, z = Z, data = test.data)
  # failing because fmatch is in transition, commentb back in later
  # expect_true(length(pairmatch(absdist1)) > 0)
 
})

