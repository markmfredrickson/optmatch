################################################################################
# Fullmatch tests
################################################################################

library(testthat)
library(optmatch)

context("fullmatch function")

test_that("No cross strata matches", {
  # test data
  Z <- rep(c(0,1), 4)
  B <- rep(c(0,1), each = 4)
  distances <- 1 + exactMatch(Z ~ B)

  res <- fullmatch(distances)
  expect_false(any(res[1:4] %in% res[5:8]))
})

test_that("Basic Matches", {
  position <- rep(1:4, each = 4)  
  z <- rep(0:1, 8)
  names(z) <- letters[1:16]
  dist <- mdist(z ~ position, inv.scale.matrix = diag(1))
  
  res.mat <- fullmatch(dist)
  res.ism <- fullmatch(as.InfinitySparseMatrix(dist))
  expect_equivalent(res.mat, res.ism)

  allin <- exactMatch(rep(1, 16), z)
  expect_equivalent(fullmatch(dist + allin), res.mat)
  expect_equivalent(fullmatch(as.InfinitySparseMatrix(dist) + allin), res.mat)

  # now that we know they are all the same, check that we got what we want
  # while this is not explicitly blocked, the position vector should
  # completely determine who matches whom (though for the same level, this
  # test is agnostic about which pairs should be matched.
  lapply(1:4, function(i) {
    in.level <- names(z)[position == i] # possibly reorders factor
    not.in.level <- names(z)[position != i]
    expect_false(any(res.mat[in.level] %in% res.mat[not.in.level]))
  })
})

test_that("Checks input", {
  # no dimnames, bad!
  m <- matrix(1:8, nrow = 2, ncol = 4)
  expect_error(fullmatch(m))
  expect_error(fullmatch(as.InfinitySparseMatrix(m)))
  # then to make sure it was the names
  dimnames(m) <- list(treated = c("A", "B"),
                      control = c("C", "D", "E", "F"))
  expect_is(fullmatch(m), "factor") # no error expected
  expect_true(all(names(fullmatch(m)) %in% LETTERS[1:6]))
  expect_false(any(is.na(fullmatch(m))))

  # add only colnames
  m <- matrix(1:8, nrow = 2, ncol = 4)
  colnames(m) <- LETTERS[3:6]  
  expect_error(fullmatch(m))
  expect_error(fullmatch(as.InfinitySparseMatrix(m)))

  # repeat for rownames
  m <- matrix(1:8, nrow = 2, ncol = 4)
  rownames(m) <- LETTERS[1:2]  
  expect_error(fullmatch(m))
  expect_error(fullmatch(as.InfinitySparseMatrix(m)))

  # a logical matrix should case an error
  ml <- matrix(rep(c(T,F), 4), nrow = 2, ncol = 2, dimnames =
    list(letters[1:2], letters[3:4]))

  expect_error(fullmatch(ml))
  ml <- replace(ml, 1:4, as.numeric(ml))
  expect_is(fullmatch(ml), "factor")

  # row and columns share names
  dimnames(ml) <- list(letters[1:2], letters[2:3])
  expect_error(fullmatch(ml))

  # the min, max, and omit must be same length as the number of subproblems,
  # which might be more than 1 if using exactMatch, e.g.

  m <- matrix(1, nrow = 2, ncol = 2, dimnames = list(c("a", "b"), c('d',
  'e')))  
  
  expect_error(fullmatch(m, min.controls = c(0,0)))
  expect_error(fullmatch(m, max.controls = c(Inf,Inf)))
  expect_error(fullmatch(m, omit.fraction = c(1, 1)))

  B <- rep(1:5, each = 2)
  names(B) <- letters[1:10]
  em <- exactMatch(B, rep(c(0,1), 5))
  expect_error(fullmatch(em, min.controls = c(0,0)))
  expect_error(fullmatch(em, max.controls = c(Inf,Inf)))
  expect_error(fullmatch(em, omit.fraction = c(1, 1)))

})

test_that("Reversion Test: Fullmatch handles omit.fraction for matrices", {
  # this bug was discovered while working on pairmatch, but it would appear to be
  # a fullmatch bug, though it might actually be in in subdivstrat or fmatch.

  A <- matrix(c(1,1,Inf,1,1,Inf,1,1,Inf,1,1,Inf), nrow = 3)
  dimnames(A) <- list(1:3, 4:7)
  Ai <- as.InfinitySparseMatrix(A)
  
  # the omit.fraction values as computed by pairmatch and was the same
  # for both A and Ai

  res.a <- fullmatch(A, min.controls = 1, max.controls = 1, omit.fraction = 0.5)
  res.ai <- fullmatch(Ai, min.controls = 1, max.controls = 1, omit.fraction = 0.5)

  expect_equivalent(res.a, res.ai)
  
})

test_that("Reversion Test: Inf entries in matrix", {
  # this was handled in the previous version just fine
  d <- matrix(c(1,2, 3,4, Inf, Inf), nrow = 2, dimnames = list(c(1,2), c(3,4, "U")))
  expect_equal(length(fullmatch(d)), 5)
  
  # the previous version also returned all 5 entries, not just the matched ones.
  expect_equal(length(fullmatch(as.InfinitySparseMatrix(d))), 5)
})

test_that("Reversion Test: Proper labeling of NAs", {
  # NA was being labeled as 1.NA
  A <- matrix(c(1,1,Inf,1,1,Inf,1,1,Inf,1,1,Inf), nrow = 3)
  dimnames(A) <- list(1:3, 4:7)

  Ai <- as.InfinitySparseMatrix(A)
  res <- fullmatch(Ai)

  expect_true(is.na(res[[3]])) # improperly labeled as "1.NA"
  expect_true(!all(is.na(res[-3])))

})

test_that("Results are in 'data order'", {
  df <- data.frame(z = rep(c(0,1), 5), x = 1:10, y = rnorm(10))
  df$w <- df$y + rnorm(10)
  rownames(df) <- letters[1:10]

  # mahal based ISM object
  m <- mdist(z ~ x + y + w, data = df)

  # make g unmatchable
  m[, "g"] <- Inf

  res <- fullmatch(m)

  expect_equal(names(res), rownames(df))

  mm <- as.InfinitySparseMatrix(m)

  res <- fullmatch(mm)

  expect_equal(names(res), rownames(df))
  
})
