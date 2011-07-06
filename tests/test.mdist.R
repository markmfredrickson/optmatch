################################################################################
# Mdist: methods to create distance specifications, possibliy sparse
################################################################################

library(testthat)

context("mdist function")

test_that("Distances from glms", {
  
  n <- 16
  Z <- numeric(n)
  Z[sample.int(n, n/2)] <- 1
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- rep(c(0,1), n/2)

  test.glm <- glm(Z ~ X1 + X2 + B, family = binomial()) # the coefs should be zero or so
   
  result.glm <- mdist(test.glm)
  
  expect_true(is(result.glm, "DistanceSpecification"))
  expect_equal(length(result.glm), (n/2)^2)
  
})

test_that("Distances from formulas", {
  n <- 16
  Z <- numeric(n)
  Z[sample.int(n, n/2)] <- 1
  X1 <- rnorm(n, mean = 5)
  X2 <- rnorm(n, mean = -2, sd = 2)
  B <- rep(c(0,1), n/2)

  test.data <- data.frame(Z, X1, X2, B)

  result.fmla <- mdist(Z ~ X1 + X2 + B, data = test.data)
  expect_true(is(result.fmla, "DistanceSpecification"))
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
  
  result.function <- mdist(sdiffs, z = Z, data = test.data)
  expect_equal(dim(result.function), c(8,8))

  # the result is a blocked matrix:
  # | 2 1 |
  # | 3 2 |

  expect_equal(mean(result.function), 2)
 
  # no treatment indicator
  expect_error(mdist(sdiffs, data = test.data))

  # no data
  expect_error(mdist(sdiffs, z = Z))
  
})

###### Using mad() instead of sd() for GLM distances
###
###result <- mdist(glm(pr ~ t1 + t2 + cost, data = nuclearplants, family = binomial()))
###
#### this is an odd test, but a simple way to make sure mad is running, not SD(). 
#### I would like a better test of the actual values, but it works
###test(mean(result$m) > 2)
###
###

test_that("Errors for numeric vectors", {
  expect_error(mdist(1:10))
})

###### Stratifying by a pipe (|) character in formulas
###
###main.fmla <- pr ~ t1 + t2
###strat.fmla <- ~ pt
###combined.fmla <- pr ~ t1 + t2 | pt
###
###result.main <- mdist(main.fmla, structure.fmla = strat.fmla, data = nuclearplants)
###result.combined <- mdist(combined.fmla, data = nuclearplants)
###
###test(identical(result.main, result.combined))
###
###### Informatively insist that one of formulas specify the treatment group
###shouldError(mdist(~t1+t2, structure.fmla=~pt, data=nuclearplants))
###test(identical(mdist(pr~t1+t2, structure.fmla=~pt, data=nuclearplants),
###               mdist(~t1+t2, structure.fmla=pr~pt, data=nuclearplants))
###     )
###### Finding "data" when it isn't given as an argument
###### Caveats:
###### * data's row.names get lost when you don't pass data as explicit argument;
###### thus testing with 'all.equal(unlist(<...>),unlist(<...>))' rather than 'identical(<...>,<...>)'.
###### * with(nuclearplants, mdist(fmla)) bombs for obscure scoping-related reasons,
###### namely that the environment of fmla is the globalenv rather than that created by 'with'.
###### This despite the facts that identical(fmla,pr ~ t1 + t2 + pt) is TRUE and that
###### with(nuclearplants, mdist(pr ~ t1 + t2 + pt)) runs fine.
###### But then with(nuclearplants, lm(fmla)) bombs too, for same reason, so don't worry be happy.
###attach(nuclearplants)
###test(all.equal(unlist(result.fmla),unlist(mdist(fmla))))
###test(all.equal(unlist(result.main),unlist(mdist(main.fmla, structure.fmla=strat.fmla))))
###test(all.equal(unlist(result.combined),unlist(mdist(combined.fmla))) )
###detach("nuclearplants")
###test(identical(fmla,pr ~ t1 + t2 + pt))
###test(all.equal(unlist(result.fmla),unlist(with(nuclearplants, mdist(pr ~ t1 + t2 + pt)))))
###test(identical(combined.fmla, pr ~ t1 + t2 | pt))
###test(all.equal(unlist(result.combined), unlist(with(nuclearplants, mdist(pr ~ t1 + t2 | pt)))))
###test(all.equal(unlist(result.fmla), unlist(with(nuclearplants[-which(names(nuclearplants)=="pt")],
###                                           mdist(update(pr ~ t1 + t2 + pt,.~.-pt + nuclearplants$pt))
###                                           )
###                                      )
###          )
###     )
###test(all.equal(unlist(result.combined), unlist(with(nuclearplants, mdist(pr ~ t1 + t2, structure.fmla=strat.fmla)))))
###
###### bigglm method
###if (require('biglm')) {
###bgps <- bigglm(fmla, data=nuclearplants, family=binomial() )
###shouldError(mdist(bgps, structure.fmla=pr ~ 1))
###shouldError(mdist(bgps, data=nuclearplants))
###result.bigglm1 <- mdist(bgps, structure.fmla=pr ~ 1, data=nuclearplants)
###result.bigglm2 <- mdist(bgps, structure.fmla=pr ~ 1, data=nuclearplants,
###                        standardization.scale=sd)
###result.bigglm2 <- mdist(bgps, structure.fmla=pr ~ 1, data=nuclearplants,
###                        standardization.scale=NULL)
###}
###
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

  absdist1 <- mdist(jb.sdiffs, z = Z, data = test.data)
  # failing because fmatch is in transition, commentb back in later
  # expect_true(length(pairmatch(absdist1)) > 0)
 
})
