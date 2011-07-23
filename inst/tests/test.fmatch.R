################################################################################
### R/Fortran Interface Tests
################################################################################

library(testthat)
context("R/Fortran Interface")

test_that("fmatch accepts DistanceSpecifications", {
  v <- c(1, Inf, 2,
         2, 1, Inf,
         3, 2, 1)

  # and doesn't accept other things...
  expect_error(fmatch(v, 2, 2))

  # the goal of this matrix is that there is a clear match to make
  # A:D, B:E, C:F
  m <- matrix(v, nrow = 3, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E", "F")

  res <- fmatch(m, 2, 2)

  expect_equal(length(res), 6)
  expect_equal(length(unique(res)), 3)
  expect_equivalent(res["A"], res["D"])
  expect_false(as.logical(res["A"] == res["B"])) # as.logical throws away the names

  M <- as.InfinitySparseMatrix(m)
  res.ism <- fmatch(M, 2, 2)
  expect_identical(res, res.ism)
})

test_that("Solutions -> factor helper", {
  v <- c(1, Inf, 2,
         2, 1, Inf,
         3, 2, 1)

  m <- matrix(v, nrow = 3, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E", "F")

  skeleton <- prepareMatching(m)

  pairs <- cbind(skeleton, solution = c(1,0,0,1,0,0,1))
  pairs.expected <- c(1,2,3,1,2,3)
  names(pairs.expected) <- c("D", "E", "F", "A", "B", "C")

  expect_equal(solution2factor(pairs), pairs.expected)

  groupOfFour <- cbind(skeleton, solution = c(1,1,0,1,1,0,1))
  gof.expected <- c(1,2,1,1,2,1)
  names(gof.expected) <- c("D", "E", "F", "A", "B", "C")
  
  expect_equal(solution2factor(groupOfFour), gof.expected)

  treatedNotMatched <- cbind(skeleton, solution = c(1,0,0,1,1,0,0))
  tnm.expected <- c(1,2, NA, 1,2,1)
  names(tnm.expected) <- c("D", "E", "F", "A", "B", "C")

  expect_equal(solution2factor(treatedNotMatched), tnm.expected)

  controlNotMatched <- cbind(skeleton, solution = c(0,0,1,1,0,0,1))
  cnm.expected <- c(1, 1, 3, NA, 1, 3)
  names(cnm.expected) <- c("D", "E", "F", "A", "B", "C")
  
  expect_equal(solution2factor(controlNotMatched), cnm.expected)
})


