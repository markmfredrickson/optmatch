################################################################################
# Tests for the optmatch object and basic methods
################################################################################
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

test_that("Subsetting drops any matched.distances attributes", {
  data(nuclearplants)

  f1 <- fullmatch(glm(pr ~ t1 + ne, data=nuclearplants, family=binomial),
                  within=exactMatch(pr ~ ne, data=nuclearplants),
                  data=nuclearplants)

  expect_true(is.null(attr(f1, "matched.distances")))

  # Add the attribute (because it is no longer created by default)
  attr(f1, "matched.distances") <- runif(length(levels(f1)))

  expect_true(!is.null(attr(f1, "matched.distances")))

  f2 <- f1[1:10]
  f3 <- f1[1:10, drop=TRUE]

  expect_true(is.null(attr(f2, "matched.distances")))
  expect_true(is.null(attr(f3, "matched.distances")))
})

test_that("Summary properly handles matched.distances #106", {
  data(nuclearplants)
  dist <- mdist(glm(pr~.-(pr+cost), family=binomial(),
                    data=nuclearplants))

  pm <- pairmatch(dist, data=nuclearplants)

  s1 <- summary(pm)

  expect_true(is.null(s1$total.distance))

  # if we add matched.distances in manually, should re-appear in
  # summary.
  attr(pm, "matched.distances") <- matched.distances(pm, dist)

  s2 <- summary(pm)

  expect_true(!is.null(s2$total.distance))
  expect_true(!is.null(s2$total.tolerance))
  expect_true(!is.null(s2$matched.dist.quantiles))

  # Double check that the match isn't getting affected.
  expect_identical(s1$thematch[sort(names(s1$thematch))], s2$thematch[sort(names(s2$thematch))])
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
  d <- data.frame(Z = c(1,0,0,0,0,1,0,0),
                  B = rep(c('a', 'b'), times=c(5, 3)))

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

  expect_warning(f <- fullmatch(res.b, data=d, max.controls=1),
                 "infeasible")
  o <- optmatch_restrictions(f)
  expect_true(all(o$min.controls == 0))
  expect_true(all(o$max.controls == 1))
  expect_true(all(o$omit.fraction == c(3/4, 1/2)))
  expect_true(all(is.null(o$mean.controls)))
})

test_that("optmatch_same_distance", {
  d <- data.frame(Z = c(1,0,0,0,0,1,0,0),
                  B = rep(c('a', 'b'), times=c(5, 3)))

  res.b <- exactMatch(Z ~ B, data=d)
  res.b2 <- res.b
  res.b2@.Data[1] <- 1

  f1 <- fullmatch(res.b, data=d)
  f2 <- fullmatch(res.b2, data=d)
  expect_warning(f3 <- fullmatch(res.b, data=d, max.controls=1),
                 "infeasible")

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
##  Z <- c(1,0,0,0,0,1,0,0)
##  B <- c(rep('a', 5), rep('b', 3))
##  d <- as.data.frame(cbind(Z,B))
##  rm(Z)
##  rm(B)

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

  a <- matrix(rep(0,9), nrow=3)
  class(a) <- c("DenseMatrix", class(a))
  expect_true(num_eligible_matches(a) == 9)

  a[1] <- Inf
  expect_true(num_eligible_matches(a) == 8)

  b <- makeInfinitySparseMatrix(1:4,
                                rows=c(1L,1L,2L,3L),
                                cols=c(1L,2L,3L,3L),
                                dimension=c(3L,3L),
                                colnames=letters[1:3],
                                rownames=LETTERS[1:3])

  expect_true(num_eligible_matches(b) == 4)
  c <- b
  c[1] <- Inf
  expect_true(num_eligible_matches(c) == 3)

  d <- as(b, "BlockedInfinitySparseMatrix")
  d@groups <- factor(c("cat","cat","dog","cat","cat","dog"))
  names(d@groups) <- c(LETTERS[1:3], letters[1:3])

  nem <- num_eligible_matches(d)
  expect_equal(names(nem), c("cat", "dog"))
  expect_equal(nem[[1]], 3)
  expect_equal(nem[[2]], 1)
  expect_true(num_eligible_matches.InfinitySparseMatrix(d) == 4)
  expect_true(num_eligible_matches(as.InfinitySparseMatrix(d)) == 4)

  e <- d
  e[1] <- Inf

  nem <- num_eligible_matches(e)
  expect_equal(names(nem), c("cat", "dog"))
  expect_equal(nem[[1]], 2)
  expect_equal(nem[[2]], 1)
  expect_true(num_eligible_matches.InfinitySparseMatrix(e) == 3)
  expect_true(num_eligible_matches(as.InfinitySparseMatrix(e)) == 3)


  d <- matrix(rep(1:2, 10), 10, 2)
  d <- caliper(d, 1.5)
  expect_true(num_eligible_matches(d) == 10)

})

test_that("equality of matches", {
  data(nuclearplants)

  # Truly identical matches
  f1 <- fullmatch(pr ~ cost, data=nuclearplants)
  f2 <- fullmatch(pr ~ cost, data=nuclearplants)

  expect_true(compare_optmatch(f1, f2))

  # Completely different matched
  data(plantdist)
  expect_warning(p1 <- fullmatch(plantdist), "Without \'data\'")

  expect_false(compare_optmatch(p1,f1))

  # Same match, different call
  f3 <- fullmatch(pr ~ cost, data=nuclearplants, max=100)

  expect_true(compare_optmatch(f1, f3))

  # Matches with unmatched objects
  expect_warning(f4 <- fullmatch(pr ~ cost, data=nuclearplants, max=1), "infeasible")
  expect_warning(f5 <- fullmatch(pr ~ cost, data=nuclearplants, max=1, min=1), "infeasible")

  expect_true(compare_optmatch(f4,f5))

  # Make sure its not returning true for everything!
  expect_false(compare_optmatch(f1,f4))
  expect_false(compare_optmatch(f3,f5))

  # Re-ordering
  nuclearplants2 <- nuclearplants[sample(1:nrow(nuclearplants)),]

  f6 <- fullmatch(pr ~ cost, data=nuclearplants2)
  # f1, f2, f6 are all the same, but f6 has a different order
  expect_true(all(f1 == f2))
  expect_false(all(f1 == f6))
  # But compare_optmatch doesn't care!
  expect_true(compare_optmatch(f1, f6))

  # Try with blocked
  b1 <- fullmatch(pr ~ cost, data=nuclearplants, within=exactMatch(pr ~ ne, data=nuclearplants))
  nuclearplants$ne2 <- 1 - nuclearplants$ne
  b2 <- fullmatch(pr ~ cost, data=nuclearplants, within=exactMatch(pr ~ ne2, data=nuclearplants))

  expect_error(all(b1 == b2), "sets of factors are different")
  # But compare_optmatch doesn't care!
  expect_true(compare_optmatch(b1, b2))

  # Make some wonky observation names
  row.names(nuclearplants) <- sapply(1:nrow(nuclearplants), function(x)
                                     paste0(sample(strsplit("!@#$%^&*()_+1234567890asdfghjkl", "")[[1]], 10, TRUE), collapse=""))

  w1 <- fullmatch(pr ~ cost, data=nuclearplants)
  w2 <- fullmatch(pr ~ cost, data=nuclearplants, max=10)
  expect_true(compare_optmatch(w1,w2))

  wb1 <- fullmatch(pr ~ cost, data=nuclearplants, within=exactMatch(pr ~ ne, data=nuclearplants))
  wb2 <- fullmatch(pr ~ cost, data=nuclearplants, within=exactMatch(pr ~ ne2, data=nuclearplants))

  expect_error(all(wb1 == wb2), "sets of factors are different")
  # But compare_optmatch doesn't care!
  expect_true(compare_optmatch(wb1, wb2))

  # If we drop NA members, should be the same match
  f4_dropna<- f4[!is.na(f4)]
  expect_true(compare_optmatch(f4, f4_dropna))

  # The problem that motivated this function: Two matches are identical, except one has an extra NA
  f4b <- f4
  f4b[1] <- NA

  # This doesn't catch it!
  expect_true(all(f4 == f4b, na.rm=TRUE))

  # This does!
  expect_false(compare_optmatch(f4, f4b))

  # Differing names should always be false.
  f1b <- f1
  names(f1b)[1] <- "Z"
  expect_false(compare_optmatch(f1, f1b))

  f1c <- f1
  names(f1c)[1] <- "A"
  expect_false(compare_optmatch(f1, f1c))

  ## # Saving this to test time.
  ## n <- 20000
  ## s1 <- as.factor(sample(letters, n, TRUE))
  ## names(s1) <- sample(LETTERS, n, TRUE)
  ## s2 <- s1[sample(seq_along(s1), n, TRUE)]
  ## system.time(compare_optmatch(s1,s2))
  ## # Taking about .3sec on laptop.

})
