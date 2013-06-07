################################################################################
# Tests for utility functions
################################################################################

context("Matching summaries")

test_that("Failing subgroups", {
  # good case:
  m <- matrix(1, nrow = 3, ncol = 4, dimnames = list(letters[1:3], LETTERS[23:26]))
  res.good <- summary(fullmatch(m))
  expect_true(all(res.good$matching.failed == FALSE))

  # good case, but one unit unmatched
  m[, 4] <- Inf

  res.one.unmatched <- summary(fullmatch(m))
  expect_true(all(res.one.unmatched$matching.failed == FALSE))

  # matching fails for all
  f <- matrix(Inf, nrow = 3, ncol = 4, dimnames = list(letters[1:3], LETTERS[23:26]))
  res.all.fail <- summary(fullmatch(f))
  expect_true(all(res.all.fail$matching.failed == TRUE))

  # blocked good, 2 2x2 groups
  b <- rep(c("A", "B"), each = 8)
  z <- rep(c(0, 1), 8)
  names(z) <- names(b) <- letters[1:16]
  em <- exactMatch(z ~ b)

  res.bg <- summary(fullmatch(em))
  expect_true(all(res.bg$matching.failed == FALSE))

  # blocked group as above but with one un-matched unit
  x <- matrix(0, nrow = 8, ncol = 8, dimnames = list(letters[c(2,4,6,8,10,12,14,16)],
                                                     letters[c(1,3,5,7,9,11,13,15)]))
  x[, "a"] <- Inf
  res.b.one.unmatched <- summary(fullmatch(em + x))
  expect_true(all(res.b.one.unmatched$matching.failed == FALSE))

  # now mock up a match in which one group failed (and there is also one unmatched unit)
  tmp <- fullmatch(em + x)
  tmp[c(letters[9:16])] <- NA
  res.b.subgrp.fail <- summary(tmp)
  expect_equal(sum(res.b.subgrp.fail$matching.failed), 4)
})
  


