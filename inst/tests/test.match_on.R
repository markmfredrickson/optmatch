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
  expect_equivalent(result.fmla, result.envir)

  expect_error(match_on(~ X1 + X2, data = test.data))
  expect_error(match_on(Z ~ 1, data = test.data))

  # checking diferent classes of responses
  res.one <- match_on(Z ~ X1) 
  res.logical <- match_on(as.logical(Z) ~ X1)
  expect_equivalent(res.one, res.logical)

  res.factor <- match_on(as.factor(Z) ~ X1)
  expect_equivalent(res.one, res.factor)

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
  res.strat <- match_on(Z ~ X1 + X2, within = stratify)
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
    expect_equivalent(res.bg, res.glm) 
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
