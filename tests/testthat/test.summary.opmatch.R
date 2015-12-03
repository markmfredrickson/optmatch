################################################################################
# Tests for utility functions
################################################################################

context("Matching summaries")

test_that("Failing subgroups", {
  # good case:
  m <- matrix(1, nrow = 3, ncol = 4, dimnames = list(letters[1:3], LETTERS[23:26]))
  res.good <- summary(fullmatch(m))
  expect_true(all(res.good$matching.failed == 0))

  # good case, but one unit unmatched
  m[, 4] <- Inf

  res.one.unmatched <- summary(fullmatch(m))
  expect_true(all(res.one.unmatched$matching.failed == 0))

  # matching fails for all
  f <- matrix(Inf, nrow = 3, ncol = 4, dimnames = list(letters[1:3], LETTERS[23:26]))
  res.all.fail <- summary(fullmatch(f))
  expect_true(all(res.all.fail$matching.failed > 0))

  # blocked good, 2 2x2 groups
  b <- rep(c("A", "B"), each = 8)
  z <- rep(c(0, 1), 8)
  names(z) <- names(b) <- letters[1:16]
  em <- exactMatch(z ~ b)

  res.bg <- summary(fullmatch(em))
  expect_true(all(res.bg$matching.failed == 0))

  # blocked group as above but with one un-matched unit
  x <- matrix(0, nrow = 8, ncol = 8, dimnames = list(letters[c(2,4,6,8,10,12,14,16)],
                                                     letters[c(1,3,5,7,9,11,13,15)]))
  x[, "a"] <- Inf
  res.b.one.unmatched <- summary(fullmatch(em + x))
  expect_true(all(res.b.one.unmatched$matching.failed == 0))

  # now mock up a match in which one group failed (and there is also one unmatched unit)
  tmp <- fullmatch(em + x)
  tmp[c(letters[9:16])] <- NA
  res.b.subgrp.fail <- summary(tmp)
  expect_true(all(res.b.subgrp.fail$matching.failed ==  c(4,4)))
})

test_that("New matching.failed", {
  data(nuclearplants)
  # one subproblem, good
  # should be NULL

  f <- fullmatch(pt ~ cost, data=nuclearplants)

  expect_true(is.null(summary(f)$matching.failed))

  # one subproblem, bad
  # should be row matrix

  f <- fullmatch(pt ~ cost, data=nuclearplants, caliper=1e-8)

  expect_true(all(summary(f)$matching.failed ==  c(26,6)))

  # many subproblems, all good
  # should be NULL

  np <- nuclearplants[nuclearplants$pt==0,]

  frame <- exactMatch(pr ~ ne + ct, data=np)
  frame@.Data[frame@rows==2] <- Inf

  m <- match_on(pr ~ cost, within=frame, data=np)

  f <- fullmatch(m, data=np)

  expect_true(is.null(summary(f)$matching.failed))

  # many subproblems, all good, some excluded individuals
  # should be empty matrix

  f <- fullmatch(pt ~ cost, data=nuclearplants, within=exactMatch(pt ~ ne, data=nuclearplants))

  expect_true(is.null(summary(f)$matching.failed))

  # many subproblems, 1 bad
  # should be row matrix

  f <- fullmatch(m, data=np)
  f[attr(f, "subproblem") == "0.0"] <- NA

  expect_true(all(summary(f)$matching.failed ==  c(7,3)))

  # many subproblems, many bad
  # should be matrix of 2 rows

  f[attr(f, "subproblem") == "0.0"] <- NA
  f[attr(f, "subproblem") == "1.0"] <- NA

  mf <- summary(f)$matching.failed
  expect_true(all(row.names(mf) == c("0.0", "1.0")))
  expect_true(all(as.numeric(mf) == c(7,3,3,1)))

  # many subproblems, all bad
  # should be table of all z's

  f[1:26] <- NA

  mf <- summary(f)$matching.failed
  expect_true(all(row.names(mf) == c("0.0", "0.1", "1.0", "1.1")))
  expect_true(all(as.numeric(mf) == c(7,6,3,3,3,2,1,1)))

  # recovered
  data(nuclearplants)
  m <- match_on(pr ~ cost, data=nuclearplants, within=exactMatch(pr ~ ct + ne, data=nuclearplants))
  m@.Data[m@rows==2] <- Inf

  f <- fullmatch(m, data=nuclearplants)
  f[1] <- NA

  # there are 5 NA's, but matching.failed only reports the 4 in the bad subgroup
  expect_equal(sum(is.na(f)), 5)
  expect_true(all(summary(f)$matching.failed ==  c(3,1)))
})
