################################################################################
# Tests for the optmatch object and basic methods
################################################################################

library(testthat)
context("Optmatch object")

test_that("Object creation", {
  dist <- diag(5)
  dimnames(dist) <- list(letters[1:5], letters[6:10])

  # recreate the result of running fullmatch. Must have err and cells fields
  ms <- list(
    list(err = 0, cells = c(a = 1, f = 1, b = 2, g = 2)),
    list(err = 0, cells = c(c = 1, h = 1, d = 2, i = 2, e = 3, j = 3)))

  res.opt <- makeOptmatch(dist, ms, NULL)

  expect_equal(length(res.opt), 10)
  expect_is(res.opt, "factor")
  expect_is(res.opt, "optmatch")

  # two levels of matches shouldn't be 1.NA, 2.NA, should just be NA
  ms2 <- list(
    list(err = 0, cells = c(a = 1, f = 1, b = 1, g = NA)),
    list(err = 0, cells = c(c = 1, h = 1, d = 2, i = 2, e = 2, j = NA)))

  res.opt2 <- makeOptmatch(dist, ms2, NULL)
  expect_true(all(is.na(res.opt2[c("g", "j")])))

})

test_that("Object subsetting", {
  dist <- diag(5)
  dimnames(dist) <- list(letters[1:5], letters[6:10])

  ms <- list(list(err = 0, cells = c(a = 1, f = 1, b = 2, g = 2)),
             list(err = 0, cells = c(c = 1, h = 1, d = 2, i = 2, e = 3, j = 3)))

  res.opt <- makeOptmatch(dist, ms, NULL)

  expect_equal(names(res.opt[1:4]), c("a", "f", "b", "g"))
  expect_equal(length(res.opt[c("a", "b")]), 2)

})

test_that("Subsetting preserves subproblem", {
  data(nuclearplants)

  # 1 subproblem
  f <- fullmatch(pr ~ cost, data=nuclearplants)

  ssf <- f[25:28]
  spssf <- attr(ssf, "subproblem")

  expect_true(all(spssf ==  attr(f, "subproblem")[25:28]))
  expect_true(all.equal(names(spssf),names(ssf)))


  # 2 subproblems
  f <- fullmatch(pr ~ cost, within=exactMatch(pr ~ pt, data=nuclearplants), data=nuclearplants)

  ssf <- f[25:28]
  spssf <- attr(ssf, "subproblem")

  expect_true(all(spssf ==  attr(f, "subproblem")[25:28]))
  expect_true(all.equal(names(spssf),names(ssf)))

  # no subproblems
  f <- fullmatch(pr ~ cost, data=nuclearplants)
  attr(f, "subproblem") <- NULL

  ssf <- f[25:28]
  spssf <- attr(ssf, "subproblem")

  expect_true(is.null(spssf))

})

test_that("Matched distances", {
  # see R/matched.distances.R for the function
  # it is only called by makeOptmatch internally, so putting the tests here
  # start with an easy case:
  dist <- matrix(Inf, nrow = 5, ncol = 5)
  diag(dist) <- 1:5
  dimnames(dist) <- list(letters[1:5], letters[6:10])
  dist.match <- as.factor(c(1.1,1.1,1.2,1.2,2.1,2.1,2.2,2.2,2.3,2.3))
  names(dist.match) <- c("a","f","b","g","c","h","d","i","e","j")
  class(dist.match) <- c("optmatch", "factor")

  res.md <- matched.distances(dist.match, dist)
  expect_equivalent(as.vector(res.md), 1:5)

  # now an ISM version
  dist.i <- as.InfinitySparseMatrix(dist)
  res.mdi <- matched.distances(dist.match, dist.i)
  expect_equivalent(as.vector(res.mdi), 1:5)

  # proper names
  res.names <- matched.distances(dist.match, dist, preserve.unit.names = TRUE)
  expect_equal(names(res.names), c("1.1", "1.2", "2.1", "2.2", "2.3"))

  res.names.i <- matched.distances(dist.match, dist.i, preserve.unit.names = TRUE)
  expect_equal(names(res.names.i), c("1.1", "1.2", "2.1", "2.2", "2.3"))

  # matches with more than one item in a strata
  match.multiple <- as.factor(c(1.1,1.1,NA,1.1,2.1,2.1,2.2,2.2,2.3,2.3))
  names(match.multiple) <- c("a","f","b","g","c","h","d","i","e","j")
  class(match.multiple) <- c("optmatch", "factor")

  dist.multiple <- dist
  dist.multiple["a", "g"] <- 99

  res.multiple <- matched.distances(match.multiple, dist.multiple, preserve.unit.names = T)
  expect_equal(length(res.multiple), 4) # 4 matches, four item list
  expect_equal(as.vector(unlist(res.multiple)), c(1, 99, 3, 4, 5))
  expect_equal(as.vector(unlist(sapply(res.multiple, names))), c("f", "g", "h", "i", "j"))


})

test_that("Match carries info about subproblems", {
  Z <- rep(c(0,1), 8)
  B <- as.factor(rep(c(1,2), each = 8))
  names(Z) <- names(B) <- letters[1:16]
  match <- pairmatch(exactMatch(Z ~ B), data = Z) # assure data order by passing Z

  # subproblem attribute should be a factor indicating which group each item maps to
  expect_equal(class(attr(match, "subproblem")), "factor")
  expect_equal(length(match), length(attr(match, "subproblem")))
  expect_equivalent(B, attr(match, "subproblem"))

})

test_that("Indicating failing subproblems", {
  Z <- rep(c(0,1), 8)
  B <- as.factor(rep(c(1,2), each = 8))
  names(Z) <- names(B) <- letters[1:16]
  match <- pairmatch(exactMatch(Z ~ B), data = Z) # assure data order by passing Z

  expect_equal(sum(subproblemSuccess(match)), 2)
  expect_true(all(names(subproblemSuccess(match)) %in%  c("1", "2")))
})

test_that("optmatch_restrictions", {
  Z <- c(1,0,0,0,0,1,0,0)
  B <- c(rep('a', 5), rep('b', 3))
  d <- as.data.frame(cbind(Z,B))

  res.b <- exactMatch(Z ~ B, data=d)

  f <- fullmatch(res.b, data=d)
  o <- optmatch_restrictions(f)
  expect_true(all(o$min.controls == 0))
  expect_true(all(o$max.controls == Inf))
  expect_true(all(is.na(o$omit.fraction)))
  expect_true(all(is.null(o$mean.controls)))

  f <- fullmatch(res.b, data=d, mean.controls = 1)
  o <- optmatch_restrictions(f)
  expect_true(all(o$min.controls == 0))
  expect_true(all(o$max.controls == Inf))
  expect_true(all(is.null(o$omit.fraction)))
  expect_true(all(o$mean.controls == 1))

  f <- fullmatch(res.b, data=d, mean.controls = 1, max.controls=c(1,2), min.controls=c(1, 1/2))
  o <- optmatch_restrictions(f)
  expect_true(all(o$min.controls == c(1, 1/2)))
  expect_true(all(o$max.controls == c(1,2)))
  expect_true(all(is.null(o$omit.fraction)))
  expect_true(all(o$mean.controls == 1))

  expect_true(all(names(o$min.controls) == c('a','b')))
  expect_true(all(names(o$max.controls) == c('a','b')))
  expect_true(all(names(o$mean.controls) == c('a','b')))

  f <- fullmatch(res.b, data=d, max.controls=1)
  o <- optmatch_restrictions(f)
  expect_true(all(o$min.controls == 0))
  expect_true(all(o$max.controls == 1))
  expect_true(all(o$omit.fraction == c(3/4, 1/2)))
  expect_true(all(is.null(o$mean.controls)))
})

test_that("optmatch_same_distance", {
  Z <- c(1,0,0,0,0,1,0,0)
  B <- c(rep('a', 5), rep('b', 3))
  d <- as.data.frame(cbind(Z,B))

  res.b <- exactMatch(Z ~ B, data=d)
  res.b2 <- res.b
  res.b2@.Data[1] <- 1

  f1 <- fullmatch(res.b, data=d)
  f2 <- fullmatch(res.b2, data=d)
  f3 <- fullmatch(res.b, data=d, max.controls=1)

  expect_true(optmatch_same_distance(f1, res.b))
  expect_true(optmatch_same_distance(f2, res.b2))
  expect_true(optmatch_same_distance(f3, res.b))

  expect_true(!optmatch_same_distance(f1, res.b2))
  expect_true(!optmatch_same_distance(f2, res.b))
  expect_true(!optmatch_same_distance(f3, res.b2))

  expect_error(optmatch_same_distance(res.b, res.b), "obj must be an optmatch object")
  expect_error(optmatch_same_distance(f1, as.matrix(res.b)), "newdist must be a valid distance")
})


## test_that("update.optmatch", {
  Z <- c(1,0,0,0,0,1,0,0)
  B <- c(rep('a', 5), rep('b', 3))
  d <- as.data.frame(cbind(Z,B))
  rm(Z)
  rm(B)

##   res.b <- exactMatch(Z ~ B, data=d)

##   f1 <- fullmatch(res.b, data=d)
##   f2 <- fullmatch(res.b, data=d, max.controls = 2)
##   f3 <- fullmatch(res.b, data=d, max.controls = 1)
##   f4 <- fullmatch(res.b, data=d, max.controls = 1, min.controls = 1)
##   f5 <- fullmatch(res.b, data=d, omit.fraction = 1/7)
##   f6 <- fullmatch(res.b, data=d, mean.controls = 1)
##   f7 <- fullmatch(res.b, data=d, tol = .00001)

##   u2 <- update(f1, max.controls=2)
##   u3 <- update(u2, max.controls=1)
##   u4 <- update(u3, min.controls=1)
##   u5 <- update(f1, omit.fraction = 1/7)
##   u6 <- update(f1, mean.controls = 1)
##   u7 <- update(f1, tol = .00001)

##   expect_true(identical(f2, u2))
##   expect_true(identical(f3, u3))
##   expect_true(identical(f4, u4))
##   expect_true(identical(f5, u5))
##   expect_true(identical(f6, u6))
##   expect_true(identical(f7, u7))

##   # update without arguments shouldn't change anything
##   f1 <- fullmatch(res.b, data=d)
##   u1 <- update(f1)
##   expect_true(identical(f1, u1))

##   f1 <- fullmatch(res.b, data=d)
##   u1 <- update(f1,data=d)
##   expect_true(identical(f1, u1))


##   # passing a difference distance
##  set.seed(9876)
##  x <- rnorm(10)
##  y <- runif(10)
##  z <- c(rep(0,6), rep(1,4))
##  d1 <- as.data.frame(cbind(x,y,z))
##  rm(x)
##  rm(y)
##  rm(z)

##   res.b1 <- match_on(z ~ x, data=d1)
##   res.b2 <- match_on(z ~ y, data=d1)

##   f1 <- fullmatch(res.b1, data = d1)
##   f2 <- fullmatch(res.b2, data = d1)

##   expect_true(!identical(f1,f2))

##   expect_warning(u1 <- update(f2, x=res.b1))
##   expect_warning(u2 <- update(f1, x=res.b2))
##   expect_true(identical(f1,u1))
##   expect_true(identical(f2,u2))
##   expect_true(!identical(f2,u1))

##   f3 <- fullmatch(res.b1, data = d1, max.controls = 2)
##   u3a <- update(f1, max.controls = 2)
##   expect_warning(u3b <- update(f2, x = res.b1, max.controls = 2))

##   expect_true(identical(f3, u3a))
##   expect_true(identical(f3, u3b))

##   # change distance between calls

##   res.c <- match_on(z ~ x, data = d2)

##   fc <- fullmatch(res.c, data=d1)

##   res.c <- match_on(z ~ y, data = d1)

##   expect_warning(uc <- update(fc, distance = res.c))

##   expect_true(!identical(fc, uc))

##   # odd ordering of parameters
##   fo <- fullmatch(data = d1, x = res.c)
##   uo <- update(fo, max.controls=2)
##   fo <- fullmatch(data = d1, x = res.c, max.controls=2)

##   expect_true(identical(fo, uo))

##   # two updates, first changing data, only one warning

##   ftu <- fullmatch(res.c, data=d1)
##   expect_warning(utu1 <- update(ftu, x=res.b))
##   expect_warning(utu2 <- update(utu1, max.controls=2), "The problem is infeasible with the given constraints; some units were omitted to allow a match.")
##   expect_warning(ftu2 <- fullmatch(res.b, data=d1, max.controls=2))
##   attr(ftu2, "call") <- NULL
##   attr(utu2, "call") <- NULL
##   expect_true(identical(ftu2, utu2))
## })

## test_that("update.optmatch with fullmatch ui simplications", {
##   set.seed(9876)
##   x <- rnorm(10)
##   y <- runif(10)
##   z <- c(rep(0,6), rep(1,4))
##   d1 <- as.data.frame(cbind(x,y,z))
##   rm(x)
##   rm(y)
##   rm(z)

##   f1 <- fullmatch(z~y+x, data=d1)
##   a <- update(f1, x=z~y, max.controls=2)

## })

test_that("num_eligible_matches", {
  x <- rnorm(10)
  y <- as.factor(rep(c("a", "b"), 5))
  y2 <- as.factor(rep(1:2, 5))
  z <- c(rep(0,6), rep(1,4))
  d1 <- as.data.frame(cbind(x,z))
  d1$y <- y
  d1$y2 <- y2
  rm(x)
  rm(y)
  rm(y2)
  rm(z)

  a <- match_on(z ~ x, data = d1)
  expect_true(num_eligible_matches(a) == 24)

  b <- caliper(a, 1e-5)
  expect_true(num_eligible_matches(b) == 0)

  c <- exactMatch(z ~ y, data=d1)
  nemc <- num_eligible_matches(c)
  expect_true(identical(nemc, list(a=as.integer(6), b=as.integer(6))))

  c2 <- exactMatch(z ~ y2, data=d1)
  nemc2 <- num_eligible_matches(c2)
  expect_true(identical(nemc2, list(`1`=as.integer(6), `2`=as.integer(6))))

  expect_true(num_eligible_matches(optmatch:::as.InfinitySparseMatrix(c)) == 12)

  d <- matrix(rep(1:2, 10), 10, 2)
  d <- caliper(d, 1.5)
  expect_true(num_eligible_matches(d) == 10)



})
