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
  dist <- match_on(glm(pr~.-(pr+cost), family=binomial(),
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

  options("optmatch_verbose_messaging" = TRUE)
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

  options("optmatch_verbose_messaging" = TRUE)
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


test_that("update.optmatch basics", {
  d <- data.frame(z = rep(0:1, each = 50),
                  b = rnorm(100))

  # update without arguments shouldn't change anything
  f1 <- fullmatch(z ~ b, data = d)
  expect_is(update(f1), "optmatch")
  expect_true(identical(f1, update(f1)))
})

test_that("update without changing distance", {
  options("optmatch_verbose_messaging" = FALSE)

  d <- data.frame(z = rep(0:1, each = 50),
                  b = rnorm(100))

  f1 <- fullmatch(z ~ b, data = d)
  f2 <- fullmatch(z ~ b, data = d, max.controls = 2)
  f3 <- fullmatch(z ~ b, data = d, max.controls = 1)
  f4 <- fullmatch(z ~ b, data = d, max.controls = 1,
                  min.controls = 1)
  f5 <- fullmatch(z ~ b, data = d, omit.fraction = 1/7)
  f6 <- fullmatch(z ~ b, data = d, mean.controls = 1)
  f7 <- fullmatch(z ~ b, data = d, tol = .00001)

  expect_true(identical(f2, update(f1, max.controls = 2)))
  expect_true(identical(f3, update(f1, max.controls = 1)))
  expect_true(identical(f4, update(f1, max.controls = 1,
                                   min.controls = 1)))
  expect_true(identical(f5, update(f1, omit.fraction = 1/7)))
  expect_true(identical(f6, update(f1, mean.controls = 1)))
  expect_true(identical(f7, update(f1, tol = .00001)))
})

test_that("upadate passing a different distance as x argument", {
  options("optmatch_verbose_messaging" = FALSE)

  # passing a difference distance
  set.seed(9876)
  d1 <- data.frame(x = rnorm(10),
                   y = runif(10),
                   z = c(rep(0,6), rep(1,4)))

  res.b1 <- match_on(z ~ x, data = d1)
  res.b2 <- match_on(z ~ y, data = d1)

  f1 <- fullmatch(res.b1, data = d1)
  f2 <- fullmatch(res.b2, data = d1)

  expect_true(!identical(as.vector(f1),as.vector(f2)))

  # When verbose messaging is off, this should produce no distance warning
  options("optmatch_verbose_messaging" = FALSE)
  u1 <- update(f2, x = res.b1)
  u2 <- update(f1, x = res.b2)
  expect_true(identical(f1,u1))
  expect_true(identical(f2,u2))
  expect_true(!identical(f2,u1))
  expect_true(!identical(as.vector(f2),as.vector(u1)))


  # If verbose messaing is enabled, should produce warning
  options("optmatch_verbose_messaging" = TRUE)
  expect_warning(update(f2, x = res.b1), "different than distance")
  expect_warning(update(f1, x = res.b2), "different than distance")
  options("optmatch_verbose_messaging" = FALSE)

  # ensure changing distance + other arguments works
  f3 <- fullmatch(res.b1, data = d1, max.controls = 2)
  u3a <- update(f1, max.controls = 2)
  u3b <- update(f2, x = res.b1, max.controls = 2)

  expect_true(identical(f3, u3a))
  expect_true(identical(f3, u3b))

})

test_that("update when distance is changed outside of update", {
  options("optmatch_verbose_messaging" = FALSE)

  set.seed(9876)
  d1 <- data.frame(x = rnorm(10),
                   y = runif(10),
                   z = c(rep(0,6), rep(1,4)))

  res.c <- match_on(z ~ x, data = d1)

  fc <- fullmatch(res.c, data = d1)

  res.c <- match_on(z ~ y, data = d1)

  uc <- update(fc, x = res.c)
  expect_true(!identical(as.vector(fc), as.vector(uc)))

  # verbose should produce warning

  options("optmatch_verbose_messaging" = TRUE)
  expect_warning(update(fc, x = res.c), "different than distance")
})

test_that("Update arguments change be ordered differently", {
  options("optmatch_verbose_messaging" = FALSE)

  set.seed(9876)
  d1 <- data.frame(x = rnorm(10),
                   y = runif(10),
                   z = c(rep(0,6), rep(1,4)))
  res.c <- match_on(z ~ y, data = d1)
  # odd ordering of parameters

  fo <- fullmatch(data = d1, x = res.c)
  uo <- update(fo, max.controls = 2)
  fo <- fullmatch(data = d1, x = res.c, max.controls = 2)

  expect_true(identical(fo, uo))
})

test_that("Update supporting new formula", {
  data(nuclearplants)

  f1 <- fullmatch(pr ~ cost, data = nuclearplants)
  f2 <- fullmatch(pr ~ t1, data = nuclearplants)

  options("optmatch_verbose_messaging" = FALSE)
  expect_error(update(f2, pr ~ cost), "must be named")
  f3 <- update(f2, x = pr ~ cost)
  expect_identical(f1, f3)
  expect_identical(update(f1, x = pr ~ cost + t1),
                   update(f2, x = pr ~ cost + t1))
})

test_that("update warning for implicit distance changes", {
  data("nuclearplants")
  p <- pairmatch(pr ~ cap, data = nuclearplants)

  # Calipering
  expect_warning(expect_is(up <- update(p, caliper = 1.5),
                           "optmatch"),
                 "different than distance")

  pcal <- pairmatch(pr ~ cap, data = nuclearplants, caliper = 1.5)
  expect_identical(up, pcal)

  # Within
  em <- exactMatch(pr ~ pt, data = nuclearplants)

  expect_warning(uem <- update(p, within = em),
                 "different than distance")

  pe <- pairmatch(pr ~ cap, data = nuclearplants, within = em)

  expect_identical(pe, uem)
})

test_that("update producing errors properly", {
  data(nuclearplants)
  f <- fullmatch(pr ~ cost, data = nuclearplants)
  call <- attr(f, "call")
  attr(f, "call") <- NULL
  expect_error(update(f), "must have a call")
  attr(f, "call") <- 7
  expect_error(update(f), "not a valid")
  attr(f, "call") <- list(call, call)
  expect_error(update(f), "combined optmatch")
})

test_that("num_eligible_matches", {

  options("optmatch_verbose_messaging" = TRUE)
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

test_that("combining optmatch objects", {
  data(nuclearplants)
  f1 <- fullmatch(pr ~ t1, data = nuclearplants[nuclearplants$pt == 0,])

  expect_is(c(f1), "optmatch")

  f2 <- fullmatch(pr ~ t1, data = nuclearplants[nuclearplants$pt == 1,])

  fc <- c(f1, f2)

  expect_equal(length(fc), length(f1) + length(f2))
  for (a in c("subproblem", "contrast.group", "levels")) {
    expect_equal(length(attr(fc, a)),
                 length(attr(f1, a)) + length(attr(f2, a)))
  }

  expect_is(attr(fc, "hashed.distance"), "list")
  expect_length(attr(fc, "hashed.distance"), 2)
  expect_is(attr(fc, "call"), "list")
  expect_length(attr(fc, "call"), 2)

  for (a in c("min.controls", "max.controls", "omit.fraction", "exceedances")) {
    expect_is(attr(fc, a), "numeric")
    expect_length(attr(fc, a), 2)
    expect_equivalent(attr(fc, a)[1], attr(f1, a))
    expect_equivalent(attr(fc, a)[2], attr(f2, a))
  }


  expect_error(c(f1, f1), "duplicated")

  full <- fullmatch(pr ~ t1, data = nuclearplants,
                    within = exactMatch(pr ~ pt, data = nuclearplants))

  expect_true(compare_optmatch(fc, full))

  levels(full) <- levels(fc)
  expect_equivalent(full, fc)

  p1 <- pairmatch(pr ~ t1, data = nuclearplants[nuclearplants$pt == 0,])

  expect_is(c(p1), "optmatch")

  p2 <- pairmatch(pr ~ t1, data = nuclearplants[nuclearplants$pt == 1,])

  pc <- c(p1, p2)

  expect_equal(length(pc), length(p1) + length(p2))
  for (a in c("subproblem", "contrast.group", "levels")) {
    expect_equal(length(attr(pc, a)),
                 length(attr(p1, a)) + length(attr(p2, a)))
  }

  expect_error(c(p1, p1), "duplicated")

  expect_identical(is.na(p1), is.na(pc)[1:26])
  expect_identical(is.na(p2), is.na(pc)[27:32])

  f1 <- fullmatch(pr ~ t1, data = nuclearplants[1:10,])
  f2 <- fullmatch(pr ~ t1, data = nuclearplants[11:25,])
  f3 <- fullmatch(pr ~ t1, data = nuclearplants[26:32,])

  fc <- c(f1, f2, f3)
  expect_is(fc, "optmatch")

  expect_equal(length(fc), length(f1) + length(f2) + length(f3))
  for (a in c("subproblem", "contrast.group", "levels")) {
    expect_equal(length(attr(fc, a)),
                 length(attr(f1, a)) + length(attr(f2, a)) +
                   length(attr(f3, a)))
  }

  expect_is(attr(fc, "hashed.distance"), "list")
  expect_length(attr(fc, "hashed.distance"), 3)
  expect_is(attr(fc, "call"), "list")
  expect_length(attr(fc, "call"), 3)

  for (a in c("min.controls", "max.controls", "omit.fraction", "exceedances")) {
    expect_is(attr(fc, a), "numeric")
    expect_length(attr(fc, a), 3)
    expect_equivalent(attr(fc, a)[1], attr(f1, a))
    expect_equivalent(attr(fc, a)[2], attr(f2, a))
    expect_equivalent(attr(fc, a)[3], attr(f3, a))
  }

  # Min, Max, etc carry forward properly
  options("optmatch_verbose_messaging" = FALSE)
  f1 <- fullmatch(pr ~ t1, data = nuclearplants[1:25,],
                  min = 1, max = 2)

  f2 <- fullmatch(pr ~ t1, data = nuclearplants[26:32,],
                  max = 3, omit.fraction = .1)

  fc <- c(f1, f2)

  expect_equivalent(attr(fc, "max.controls"),
                    c(attr(f1, "max.controls"),
                      attr(f2, "max.controls")))
  expect_equivalent(attr(fc, "min.controls"),
                    c(attr(f1, "min.controls"),
                      attr(f2, "min.controls")))
  expect_equivalent(attr(fc, "omit.fraction"),
                    c(attr(f1, "omit.fraction"),
                      attr(f2, "omit.fraction")))

  # Functions taking optmatch objects

  f1 <- fullmatch(pr ~ t1, data = nuclearplants[1:25,],
                  min = 1, max = 2)
  f2 <- fullmatch(pr ~ t1, data = nuclearplants[26:32,],
                  min = 1, max = 2)
  fc <- c(f1, f2)
  nuclearplants$treat <- rep(0:1, times = c(25, 7))
  full <- fullmatch(pr ~ t1, data = nuclearplants, min = 1, max = 2,
                    within = exactMatch(pr ~ treat, data = nuclearplants))

  expect_true(compare_optmatch(fc, full))
  expect_identical(matched(fc), matched(full))

  expect_identical(optmatch_restrictions(fc), optmatch_restrictions(full))
  expect_identical(stratumStructure(fc), stratumStructure(full))
  expect_identical(summary(fc)$effective.sample.size,
                   summary(full)$effective.sample.size)
  expect_identical(summary(fc)$matched.set.structures,
                   summary(full)$matched.set.structures)
})
