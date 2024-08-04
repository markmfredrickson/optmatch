################################################################################
# Testing resolution / tolerance specification
################################################################################

context("test subproblem resolution specification functionality")

test_that("get_subproblem_info basics", {
  d <- data.frame(position = rep(1:4, each = 4),
                  z = rep(0:1, 8),
                  rownames=letters[1:16])
  dist <- match_on(z ~ position, inv.scale.matrix = diag(1), data=d)

  res.mat <- fullmatch(dist, data=d)
  expect_true(identical(length(unique(get_subproblem_info(res.mat, type = "subproblem"))), 1L))
  expect_true(identical(length(get_subproblem_info(res.mat, type = "resolution")), 1L))

  data(nuclearplants)
  f2 <- fullmatch(pr ~ cost + strata(pt), data=nuclearplants)

  expect_true(identical(length(unique(get_subproblem_info(f2, type = "subproblem"))), 2L))
  expect_true(identical(length(get_subproblem_info(f2, type = "resolution")), 2L))
})

test_that("fullmatch: basic error checks", {
  data(nuclearplants)
  res.custom <- list("0" = 0.001142857, "1" = 0.001142857)
  expect_silent(fullmatch(pr ~ cost + strata(pt), data=nuclearplants, resolution = res.custom))
  bad.res <- list("2" = 0.001142857, "1" = 0.001142857)
  # subproblem is not in the original specification
  expect_error(fullmatch(pr ~ cost + strata(pt), data=nuclearplants, resolution = bad.res))
  # no hint + subset of epsilons should fail
  bad.res <- list("1" = 0.001142857)
  expect_error(fullmatch(pr ~ cost + strata(pt), data=nuclearplants, resolution = bad.res))
  #unnamed arguments should fail, length 1 list is a special case
  bad.res <- list(0.001142857)
  expect_error(fullmatch(pr ~ cost + strata(pt), data=nuclearplants, resolution = bad.res))
  #unnamed arguments should fail
  bad.res <- list(0.001142857, "1" =.35)
  expect_error(fullmatch(pr ~ cost + strata(pt), data=nuclearplants, resolution = bad.res))

  #unnamed arguments should fail
  bad.res <- list(0.001142857, .35)
  expect_error(fullmatch(pr ~ cost + strata(pt), data=nuclearplants, resolution = bad.res))

})

test_that("pairmatch: basic error checks", {
  data(nuclearplants)

  psm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
  em <- exactMatch(pr ~ pt, data = nuclearplants)
  res.pm <- pairmatch(match_on(psm) + em, data=nuclearplants)
  res.custom <- list("0" = 0.001142857, "1" = 0.001142857)
  expect_silent(pairmatch(match_on(psm) + em, data=nuclearplants, resolution = res.custom))

  bad.res <- list("2" = 0.001142857, "1" = 0.001142857)
  expect_error(pairmatch(match_on(psm) + em, data=nuclearplants, resolution = bad.res))

  bad.res <- list("1" = 0.001142857)
  expect_error(pairmatch(match_on(psm) + em, data=nuclearplants, resolution = bad.res))


  bad.res <- list(0.001142857, 0.001142857)
  expect_error(pairmatch(match_on(psm) + em, data=nuclearplants, resolution = bad.res))

  bad.res <- list("1" = 0.001142857,
                  0.001142857)
  expect_error(pairmatch(match_on(psm) + em, data=nuclearplants, resolution = bad.res))

  bad.res <- list(0.001142857)
  expect_error(pairmatch(match_on(psm) + em, data=nuclearplants, resolution = bad.res))

})

test_that("fullmatch: error checking with hints", {
  set.seed(201905)
  data <- data.frame(z = rep(0:1, each = 5),
                     x = rnorm(10), fac=rep(c(rep("a",2), rep("b",3)),2) )
  mo  <- match_on(z ~ x, data=data)
  mos <- match_on(z ~ x + strata(fac), data=data)
  f1b <- fullmatch(mos, min.c=.5, max.c=2, data = data, tol=0.1)
  res.custom <- list("b" = 0.0001666667)
  expect_silent(fullmatch(mos, min.c=.5, max.c=2, data = data, tol=0.0001, hint=f1b, resolution = res.custom))

  res.custom <- list("b" = 0.0001666667, "a" = 0.0001666667) #out of order should be ok
  expect_silent(fullmatch(mos, min.c=.5, max.c=2, data = data, tol=0.0001, hint=f1b, resolution = res.custom))

  res.custom <- list("b" = 0.01)
  expect_silent(fullmatch(mos, min.c=.5, max.c=2, data = data, hint=f1b, resolution = res.custom)  )
})


test_that("fullmatch: with hint", {
  set.seed(201905)
  data <- data.frame(z = rep(0:1, each = 5),
                     x = rnorm(10), fac=rep(c(rep("a",2), rep("b",3)),2) )
  mo  <- match_on(z ~ x, data=data)

  mos <- match_on(z ~ x + strata(fac), data=data)
  f1b <- fullmatch(mos, min.c=.5, max.c=2, data = data, tol=0.1)

  f1c <- fullmatch(mos, min.c=.5, max.c=2, data = data, tol=0.0001, hint=f1b)
  res.custom <- list("b" = 0.0001666667)
  f1c <- fullmatch(mos, min.c=.5, max.c=2, data = data, hint=f1b, resolution = res.custom)
  resos <- get_subproblem_info(f1c, "resolution")
  expect_true(resos["b"] == 0.0001666667)

  res.custom <- list("b" = 0.01)
  f1c <- fullmatch(mos, min.c=.5, max.c=2, data = data, hint=f1b, resolution = res.custom)
  resos <- get_subproblem_info(f1c, "resolution")
  expect_true(resos["b"] == 0.01)

  res.custom <- list("b" = 0.01, "a" = .001)
  f1c <- fullmatch(mos, min.c=.5, max.c=2, data = data, hint=f1b, resolution = res.custom)
  resos <- get_subproblem_info(f1c, "resolution")
  expect_true(resos["b"] == 0.01)
  expect_true(resos["a"] == 0.001)
})

test_that("fullmatch: subproblem specification", {
  data(nuclearplants)
  f1 <- fullmatch(pr ~ cost + strata(pt), data=nuclearplants)
  res.custom <- list("0" = 0.001, "1" = 0.01)
  f2 <- fullmatch(pr ~ cost + strata(pt), data=nuclearplants, resolution = res.custom)
  resos.old <- get_subproblem_info(f1, "resolution")
  resos <- get_subproblem_info(f2, "resolution")
  expect_true(resos["0"] == 0.001)
  expect_true(resos["1"] == 0.01)
  expect_true(all(resos.old != resos))
})


test_that("pairmatch: error checking with hints", {
  set.seed(201905)
  data <- data.frame(z = rep(0:1, each = 5),
                     x = rnorm(10), fac=rep(c(rep("a",2), rep("b",3)),2) )
  mo  <- match_on(z ~ x, data=data)
  mos <- match_on(z ~ x + strata(fac), data=data)
  p1b <- pairmatch(mos, data = data, tol=0.1)

  res.custom <- list("b" = 0.001, "a" = 0.01) #out of order should be ok
  expect_silent(pairmatch(mos, data = data, hint=p1b, resolution = res.custom))
  res.custom <- list("b" = 0.001)
  expect_silent(pairmatch(mos, data = data, hint=p1b, resolution = res.custom))
})

test_that("pairmatch: with hint", {
  set.seed(201905)
  data <- data.frame(z = rep(0:1, each = 5),
                     x = rnorm(10), fac=rep(c(rep("a",2), rep("b",3)),2) )
  mo  <- match_on(z ~ x, data=data)
  mos <- match_on(z ~ x + strata(fac), data=data)
  p1b <- pairmatch(mos, data = data, tol=0.1)

  res.custom <- list("b" = 0.001, "a" = 0.01) #out of order should be ok
  pm1 <- pairmatch(mos, data = data, hint=p1b, resolution = res.custom)
  resos <- get_subproblem_info(pm1, "resolution")
  expect_true(resos["b"] == 0.001)
  expect_true(resos["a"] == 0.01)

  res.custom <- list("b" = 0.001) #should be able to specify subset here
  pm1 <- pairmatch(mos, data = data, hint=p1b, resolution = res.custom)
  resos <- get_subproblem_info(pm1, "resolution")
  expect_true(resos["b"] == 0.001)

})

test_that("pairmatch: subproblem specification", {
  data(nuclearplants)
  psm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
  em <- exactMatch(pr ~ pt, data = nuclearplants)
  res.pm <- pairmatch(match_on(psm) + em, data=nuclearplants)
  res.custom <- list("0" = 0.001,
                     "1" = 0.01)
  pm1 <- pairmatch(match_on(psm) + em, data=nuclearplants, resolution = res.custom)
  resos.old <- get_subproblem_info(res.pm, "resolution")
  resos <- get_subproblem_info(pm1, "resolution")
  expect_true(resos["0"] == 0.001)
  expect_true(resos["1"] == 0.01)
  expect_true(all(resos.old != resos))
})


test_that("fullmatch: no subproblems", {
  d <- data.frame(position = rep(1:4, each = 4),
                  z = rep(0:1, 8),
                  rownames=letters[1:16])
  dist <- match_on(z ~ position, inv.scale.matrix = diag(1), data=d)
  res.mat <- fullmatch(dist, data=d, resolution = list(.1))
  resos <- get_subproblem_info(res.mat, "resolution")
  expect_true(length(resos) == 1)
  expect_true(all(resos == .1))
})

test_that("pairmatch: no subproblems", {
  data(nuclearplants)

  psm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
  res.pm <- pairmatch(match_on(psm), data=nuclearplants, resolution = list(.1))
  resos <- get_subproblem_info(res.pm, "resolution")
  expect_true(length(resos) == 1)
  expect_true(all(resos == .1))
})
