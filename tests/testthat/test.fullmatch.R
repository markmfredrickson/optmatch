################################################################################
# Fullmatch tests
################################################################################

context("fullmatch function")

# test whether two matches are the same. Uses all.equal exceedances to
# ignore errors below some tolerance. After checking those, strips
# attributes that may differ but not break `identical` status.
match_compare <- function(match1, match2) {
  expect_true(all.equal(attr(match1, "exceedances"),
                        attr(match2, "exceedances")))

  attr(match1, "hashed.distance") <- NULL
  attr(match2, "hashed.distance") <- NULL
  attr(match1, "exceedances") <- NULL
  attr(match2, "exceedances") <- NULL
  attr(match1, "call") <- NULL
  attr(match2, "call") <- NULL

  expect_true(identical(match1, match2))
}

test_that("No cross strata matches", {
  # test data
  d <- data.frame(Z = rep(c(0,1), 4),
                  B = rep(c(0,1), each = 4))
  distances <- 1 + exactMatch(Z ~ B, data=d)

  res <- fullmatch(distances, data=d)
  expect_false(any(res[1:4] %in% res[5:8]))
})

test_that("Basic Matches", {
  d <- data.frame(position = rep(1:4, each = 4),
                  z = rep(0:1, 8),
                  rownames=letters[1:16])
  dist <- match_on(z ~ position, inv.scale.matrix = diag(1), data=d)

  res.mat <- fullmatch(dist, data=d)
  res.ism <- fullmatch(as.InfinitySparseMatrix(dist), data=d)
  expect_equivalent(res.mat, res.ism)

  allin <- exactMatch(rep(1, 16), structure(d$z, names=rownames(d)))
  expect_equivalent(fullmatch(dist + allin, data=d), res.mat)
  expect_equivalent(fullmatch(as.InfinitySparseMatrix(dist) + allin, data=d), res.mat)

  # Now that we know they are all the same, check that we got what we
  # want.  While this is not explicitly blocked, the position vector
  # should completely determine who matches whom (though for the same
  # level, this test is agnostic about which pairs should be matched.
  for (i in 1:4) {
    in.level <- rownames(d)[d$position == i] # possibly reorders factor
    not.in.level <- rownames(d)[d$position != i]
    expect_false(any(res.mat[in.level] %in% res.mat[not.in.level]))
  }
})

test_that("Checks input", {
  # no dimnames, bad!
  m <- matrix(1:8, nrow = 2, ncol = 4)
  # Throw both an error and a warning about lack of data
  expect_warning(expect_error(fullmatch(m)))
  expect_warning(expect_error(fullmatch(as.InfinitySparseMatrix(m))))
  # then to make sure it was the names
  dimnames(m) <- list(treated = c("A", "B"),
                      control = c("C", "D", "E", "F"))
  # Still should warn about lack of data
  expect_warning(expect_is(fullmatch(m), "factor")) # no error expected
  expect_warning(expect_true(all(names(fullmatch(m)) %in% LETTERS[1:6])))
  expect_warning(expect_false(any(is.na(fullmatch(m)))))

  # add only colnames
  m <- matrix(1:8, nrow = 2, ncol = 4)
  colnames(m) <- LETTERS[3:6]
  expect_warning(expect_error(fullmatch(m)))
  expect_warning(expect_error(fullmatch(as.InfinitySparseMatrix(m))))

  # repeat for rownames
  m <- matrix(1:8, nrow = 2, ncol = 4)
  rownames(m) <- LETTERS[1:2]
  expect_warning(expect_error(fullmatch(m)))
  expect_warning(expect_error(fullmatch(as.InfinitySparseMatrix(m))))

  # a logical matrix should case an error
  ml <- matrix(rep(c(T,F), 4), nrow = 2, ncol = 2, dimnames =
    list(letters[1:2], letters[3:4]))

  expect_warning(expect_error(fullmatch(ml)))
  ml <- replace(ml, 1:4, as.numeric(ml))
  expect_warning(expect_is(fullmatch(ml), "factor"))

  # row and columns share names
  dimnames(ml) <- list(letters[1:2], letters[2:3])
  expect_warning(expect_error(fullmatch(ml)))

  # the min, max, and omit must be same length as the number of
  # subproblems, which might be more than 1 if using exactMatch, e.g.

  m <- matrix(1, nrow = 2, ncol = 2, dimnames = list(c("a", "b"),
                                                     c('d', 'e')))

  expect_warning(expect_error(fullmatch(m, min.controls = c(0,0))))
  expect_warning(expect_error(fullmatch(m, max.controls = c(Inf,Inf))))
  expect_warning(expect_error(fullmatch(m, omit.fraction = c(1, 1))))

  B <- rep(1:5, each = 2)
  names(B) <- letters[1:10]
  em <- exactMatch(B, rep(c(0,1), 5))
  expect_warning(expect_error(fullmatch(em, min.controls = c(0,0))))
  expect_warning(expect_error(fullmatch(em, max.controls = c(Inf,Inf))))
  expect_warning(expect_error(fullmatch(em, omit.fraction = c(1, 1))))

})

test_that("fullmatch warns when given a 'within' arg that it's going to ignore", {
    m <- matrix(1, nrow = 2, ncol = 3,
                dimnames = list(c("a", "b"), c('d', 'e', 'f')))
    B <- rep(1:3, each = 2)
    names(B) <- letters[1:6]
    em <- exactMatch(B, rep(c(0,1), 3))
    expect_warning(fullmatch(m, within=em), "gnor")
    expect_warning(fullmatch(as.InfinitySparseMatrix(m), within=em), "gnor")
})

test_that("Reversion Test: Fullmatch handles omit.fraction for matrices", {
  # this bug was discovered while working on pairmatch, but it would appear to be
  # a fullmatch bug, though it might actually be in in subdivstrat or fmatch.

  A <- matrix(c(1,1,Inf,1,1,Inf,1,1,Inf,1,1,Inf), nrow = 3)
  dimnames(A) <- list(1:3, 4:7)
  Ai <- as.InfinitySparseMatrix(A)

  # the omit.fraction values as computed by pairmatch and was the same
  # for both A and Ai

  res.a <- fullmatch(A, min.controls = 1, max.controls = 1, omit.fraction = 0.5, data=data.frame(1:7))
  res.ai <- fullmatch(Ai, min.controls = 1, max.controls = 1, omit.fraction = 0.5, data=data.frame(1:7))

  expect_equivalent(res.a, res.ai)

})

test_that("Reversion Test: Inf entries in matrix", {
  # this was handled in the previous version just fine
  d <- matrix(c(1,2, 3,4, Inf, Inf), nrow = 2, dimnames = list(c(1,2), c(3,4, "U")))
  expect_warning(expect_equal(length(fullmatch(d)), 5))

  # the previous version also returned all 5 entries, not just the matched ones.
  expect_warning(expect_equal(length(fullmatch(as.InfinitySparseMatrix(d))), 5))
})

test_that("Reversion Test: Proper labeling of NAs", {
  # NA was being labeled as 1.NA
  A <- matrix(c(1,1,Inf,1,1,Inf,1,1,Inf,1,1,Inf), nrow = 3)
  dimnames(A) <- list(1:3, 4:7)

  Ai <- as.InfinitySparseMatrix(A)
  res <- fullmatch(Ai, data=data.frame(1:7))

  expect_true(is.na(res[[3]])) # improperly labeled as "1.NA"
  expect_true(!all(is.na(res[-3])))

  # also, entirely failing problems should be labeled with NAs, not subgroup.NA
  # https://github.com/markmfredrickson/optmatch/issues/22
  scores <- c(10,1,10,10,10,10,10,10,10, 1,1,1,1,1,10,10,10)
  B <- c(rep(1, 9), rep(2, 8))
  Z <- c(rep(c(0,1), each = 4), 0, rep(c(0,1), each = 4))

  names(scores) <- names(B) <- names(Z) <- letters[1:17]

  d <- match_on(scores, z = Z, within = exactMatch(Z ~ B))
  expect_warning(res <- pairmatch(caliper(d, 2)))

  expect_equal(sum(is.na(res)), 9)

  # repeat using an optmatch.dlist object, just in case...
  od <- list(matrix(c(0,0,0,0, Inf,Inf,Inf,Inf, rep(0, 12)), nrow = 4, ncol = 5, dimnames = list(letters[1:4], letters[5:9])),
             matrix(c(0,0,0,0, rep(Inf, 12)), byrow = T, nrow = 4, ncol = 4, dimnames = list(letters[10:13], letters[14:17])))

  class(od) <- c("optmatch.dlist", "list")

  expect_warning(expect_equal(sum(is.na(pairmatch(od))), 9))

  # while we're at it, check that match failed only indicates that 1 level failed
  expect_warning(expect_equal(sum(matchfailed(pairmatch(od))), 8))
})

test_that("Results are in 'data order'", {
  # https://github.com/markmfredrickson/optmatch/issues/14

  df <- data.frame(z = rep(c(0,1), 5), x = 1:10, y = rnorm(10))
  df$w <- df$y + rnorm(10)
  rownames(df) <- letters[1:10][sample(1:10)]

  # add some NAs to the df:
  df[3, "y"] <- NA

  # mahal based ISM object
  m <- as.matrix(match_on(z ~ x + y + w, data = df))

  # make the first control unit unmatchable
  m[, 1] <- Inf

  res <- fullmatch(m, data = df)

  expect_equal(names(res), rownames(df))

  # shuffle the order of the distance matrix
  m2 <- m[sample(1:5),]
  res <- fullmatch(m2, data = df)

  expect_equal(names(res), rownames(df))

  mm <- as.InfinitySparseMatrix(m)

  res <- fullmatch(mm, data = df)

  expect_equal(names(res), rownames(df))

  # not supplying a data argument is grounds for a warning
  expect_warning(fullmatch(mm), "data")

  # data argument should have useful names attached
  tmp <- as.matrix(df)
  rownames(tmp) <- colnames(tmp) <- NULL
  expect_error(fullmatch(mm, data = tmp), "are not found")

  # Catching issue #56 rather than returning nonsense
  expect_error(fullmatch(mm, data=mm), "are not found")

})

test_that("Complete Inf matrices/ISMs => all NA optmatch object", {
  # Issue 15: https://github.com/markmfredrickson/optmatch/issues/15

  m <- matrix(Inf, nrow = 3, ncol = 4)
  rownames(m) <- LETTERS[1:3]
  colnames(m) <- letters[23:26]

  expect_warning(res.m <- fullmatch(m))

  expect_true(all(is.na(res.m)))

  ism <- as.InfinitySparseMatrix(m)

  expect_warning(res.ism <- fullmatch(ism))

  expect_true(all(is.na(res.ism)))

})

test_that("Both mdist and match_on objects accepted", {
  # this test depends on mdist and match_on tests passing
  # it will probably fail if those fail

  n <- 14
  test.data <- data.frame(Z = c(rep(0, n/2), rep(1, n/2)),
                          X1 = rnorm(n, mean = 5),
                          X2 = rnorm(n, mean = -2, sd = 2),
                          B = rep(c(0,1), n/2))

  model <- glm(Z ~ X1 + X2, data = test.data, family = binomial())
  tmp <- mdist(model)
  names(tmp) <- c(1) # mdist adds an 'm' to the front by default
  res.mdist <- fullmatch(tmp, data=test.data)
  res.mon <- fullmatch(match_on(model), data=test.data)

  expect_equivalent(res.mdist, res.mon)

})

test_that("full() and pair() are alises to _match functions", {

  n <- 14
  test.data <- data.frame(Z = c(rep(0, n/2), rep(1, n/2)),
                          X1 = rnorm(n, mean = 5),
                          X2 = rnorm(n, mean = -2, sd = 2),
                          B = rep(c(0,1), n/2))

  model <- glm(Z ~ X1 + X2, data = test.data, family = binomial())
  dists <- match_on(model)
  expect_equivalent(fullmatch(dists, data=test.data),
                    full(dists, data=test.data))
  expect_equivalent(pairmatch(dists, data=test.data),
                    pair(dists, data=test.data))
})

test_that("fullmatch UI cleanup", {
  n <- 14
  test.data <- data.frame(Z = c(rep(0, n/2), rep(1, n/2)),
                          X1 = rnorm(n, mean = 5),
                          X2 = rnorm(n, mean = -2, sd = 2),
                          B = rep(c(0,1), n/2))

  m <- match_on(Z~X1 + X2, within=exactMatch(Z~B, data=test.data), data=test.data, caliper=2)

  fm.dist <- fullmatch(m, data=test.data)

  fm.form <- fullmatch(Z~X1 + X2, within=exactMatch(Z~B, data=test.data), data=test.data, caliper=2)

  match_compare(fm.dist, fm.form)

  # with "with()"

  fm.with <- with(data=test.data, fullmatch(Z~X1 + X2, within=exactMatch(Z~B), caliper=2))

  match_compare(fm.dist, fm.with)

  # passing a glm
  ps <- glm(Z~X1+X2, data=test.data, family=binomial)

  fm.ps <- fullmatch(ps, data=test.data, caliper=2)

  fm.glm <- fullmatch(glm(Z~X1+X2, data=test.data, family=binomial), data=test.data, caliper=2)
  fm.glm2 <- fullmatch(glm(Z~X1+X2, data=test.data, family=binomial), caliper=2)

  match_compare(fm.ps, fm.glm)

  # passing inherited from glm

  class(ps) <- c("foo", class(ps))

  fm.foo <- fullmatch(ps, data=test.data, caliper=2)

  match_compare(fm.ps, fm.foo)

  # with scores

  ps <- glm(Z~X2, data=test.data, family=binomial)

  m <- match_on(Z ~ X1 + scores(ps), within=exactMatch(Z~B, data=test.data), data=test.data)

  fm.dist <- fullmatch(m, data=test.data)

  fm.form <- fullmatch(Z~ X1 + scores(ps), within=exactMatch(Z~B, data=test.data), data=test.data)

  match_compare(fm.dist, fm.form)

  # passing numeric

  X1 <- test.data$X1
  Z <- test.data$Z

  names(X1) <- row.names(test.data)
  names(Z) <- row.names(test.data)
  fm.vector <- fullmatch(X1,z=Z, data=test.data, caliper=1)
  expect_warning(fm.vector2 <- fullmatch(X1,z=Z, caliper=1))

  m <- match_on(X1, z=Z, caliper=1)
  fm.mi <- fullmatch(m, data=test.data)

  match_compare(fm.vector, fm.mi)

  # function

  n <- 16
  test.data <- data.frame(Z = c(rep(0, n/2), rep(1, n/2)),
                          X1 = rep(c(1,2,3,4), each = n/4),
                          B = rep(c(0,1), n/2))

  sdiffs <- function(index, data, z) {
    abs(data[index[,1], "X1"] - data[index[,2], "X1"])
  }

  result.function <- match_on(sdiffs, z = test.data$Z, data = test.data)

  fm.funcres <- fullmatch(result.function, data=test.data)

  fm.func <- fullmatch(sdiffs, z = test.data$Z, data=test.data)
  expect_error(fullmatch(sdiffs, z = Z), "A data argument must be given when passing a function")

  match_compare(fm.funcres, fm.func)

  # passing bad arguments

  expect_error(fullmatch(test.data), "Invalid input, must be a potential argument to match_on")
  expect_error(fullmatch(TRUE), "Invalid input, must be a potential argument to match_on")


})

test_that("NAs in irrelevant data slots don't trip us up", {
  n <- 16
  test.data <- data.frame(Z = c(rep(0, n/2), rep(1, n/2)),
                          X1 = rep(c(1,2,3,4), each = n/4),
                          B = rep(c(0,1), n/2))
  test.data$B[1] <- NA

  expect_equal(length(fullmatch(Z~X1, data=test.data)), n)
})

test_that("Using strata instead of within arguments", {
  data(nuclearplants)

  f1 <- fullmatch(pr ~ cost, within=exactMatch(pr ~ pt, data=nuclearplants),
                  data=nuclearplants)
  f2 <- fullmatch(pr ~ cost + strata(pt), data=nuclearplants)
  f2b <- fullmatch(pr ~ cost, data=nuclearplants)

  expect_true(is(f1, "optmatch"))
  expect_true(is(f2, "optmatch"))
  expect_true(compare_optmatch(f1, f2))
  expect_false(compare_optmatch(f2, f2b))

  # handling more complicated strata calls
  f3 <- fullmatch(pr ~ cost, within=exactMatch(pr ~ pt + ct + ne,
                                               data=nuclearplants),
                  data=nuclearplants)
  f4 <- fullmatch(pr ~ cost + strata(pt) + strata(ct, ne), data=nuclearplants)

  expect_true(compare_optmatch(f3, f4))

  e1 <- exactMatch(pr ~ ne, data=nuclearplants)
  e2 <- exactMatch(pr ~ ne + ct, data=nuclearplants)

  f5 <- fullmatch(pr ~ cost, within=e2, data=nuclearplants)
  f6 <- fullmatch(pr ~ cost + strata(ct), within=e1, data=nuclearplants)
  f7 <- fullmatch(pr ~ cost + strata(ct, ne), data=nuclearplants)

  expect_true(compare_optmatch(f5, f6))
  expect_true(compare_optmatch(f5, f7))

})

test_that("strata in GLMs", {

  data(nuclearplants)

  f1 <- fullmatch(glm(pr ~ t1 + ne, data=nuclearplants, family=binomial),
                  within=exactMatch(pr ~ ne, data=nuclearplants),
                  data=nuclearplants)
  f2 <- fullmatch(glm(pr ~ t1 + strata(ne),
                      data=nuclearplants, family=binomial),
                  data=nuclearplants)

  expect_true(is(f1, "optmatch"))
  expect_true(is(f2, "optmatch"))
  expect_true(compare_optmatch(f1, f2))

  f3 <- fullmatch(glm(pr ~ t1 + ne + interaction(ct,pt),
                      data=nuclearplants, family=binomial),
                  within=exactMatch(pr ~ ne + ct*pt, data=nuclearplants),
                  data=nuclearplants)
  f4 <- fullmatch(glm(pr ~ t1 + strata(ne) + strata(ct, pt),
                      data=nuclearplants, family=binomial),
                  data=nuclearplants)

  expect_true(compare_optmatch(f3, f4))

  e1 <- exactMatch(pr ~ ne, data=nuclearplants)
  e2 <- exactMatch(pr ~ ne + ct, data=nuclearplants)

  f5 <- fullmatch(glm(pr ~ cost + ne + ct, data=nuclearplants),
                  within=e2, data=nuclearplants)
  f6 <- fullmatch(glm(pr ~ cost + ne + strata(ct), data=nuclearplants),
                  within=e1, data=nuclearplants)
  f7 <- fullmatch(glm(pr ~ cost + strata(ct) + strata(ne), data=nuclearplants),
                  data=nuclearplants)

  expect_true(compare_optmatch(f5, f6))
  expect_true(compare_optmatch(f5, f7))

  # strata(a,b) is equivalent to interaction(a,b)
  f8 <- fullmatch(glm(pr ~ cost + strata(ct,ne), data=nuclearplants),
                  data=nuclearplants)
  suppressWarnings(
    f9 <- fullmatch(glm(pr ~ cost + interaction(ne, ct) + strata(ct),
                        data=nuclearplants),
                    within=e1, data=nuclearplants)
  )
  f10 <- fullmatch(glm(pr ~ cost + interaction(ne,ct), data=nuclearplants),
                   within=e2, data=nuclearplants)
  f11 <- fullmatch(glm(pr ~ cost + ne*ct, data=nuclearplants),
                   within=e2, data=nuclearplants)
  # f9 is a bit weird because of the double inclusion of ct, and is an
  # unlikely way for users to enter code, but the extra ct is of
  # course ignored.

  expect_true(compare_optmatch(f8, f9))
  expect_true(compare_optmatch(f8, f10))
  expect_true(compare_optmatch(f10, f11))

})

test_that("matched.distances attr removed per #57", {
  data(nuclearplants)

  f1 <- fullmatch(glm(pr ~ t1 + ne, data=nuclearplants, family=binomial),
                  within=exactMatch(pr ~ ne, data=nuclearplants),
                  data=nuclearplants)

  expect_true(is.null(attr(f1, "matched.distances")))
})
