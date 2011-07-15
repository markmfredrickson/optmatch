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
    in.level <- position == i
    expect_false(any(res.mat[in.level] %in% res.mat[!in.level]))
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
})


