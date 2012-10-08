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
