################################################################################
# Fullmatch tests
################################################################################

library(testthat)
library(optmatch)

context("feasibility")

test_that("Problems bigger than the max are not feasible", {
  # the options('optmatch_max_problem_size') should exist and be less than 1e100
  expect_true(options("optmatch_max_problem_size")[[1]] < 1e100)

  # the helper function getMaxProblemSize should make the above easier
  expect_equal(options("optmatch_max_problem_size")[[1]], getMaxProblemSize())

  # now we use the options() function to set this limit(!) much lower and create a problem that should pass
  options('optmatch_max_problem_size' = 256)
  maxprob <- getMaxProblemSize()
  
  feasible.n <- floor(sqrt(maxprob)) - 1

  largeButFeasible <- matrix(0, nrow = feasible.n, ncol = feasible.n,
    dimnames = list(1:feasible.n, (feasible.n + 1):(2 * feasible.n)))
  
  fullmatch(largeButFeasible)

  # now we make too big a problem and try again
  infeasible <- matrix(0, nrow = feasible.n, ncol = 2 * feasible.n,
    dimnames = list(1:feasible.n, (feasible.n + 1):(3 * feasible.n)))

  expect_error(fullmatch(infeasible), "too many")

  setFeasibilityConstants() # reset the values to make sure that other tests pass
})

test_that("minExactMatch creates minimal exact match", {
  # no subproblem can be bigger than 8x8
  oldopts <- options("optmatch_max_problem_size" = 37) 

  df <- data.frame(Z = rep(c(1,0), 16),
                   E1 = rep(c(1,1,0,0,0,0,0,0), each = 4), # cuts size in 1/2, too big still
                   E2 = rep(c(1,1,0,0), 8),
                   E3 = rep(c(1,1,1,1,0,0,0,0), 4))
  
  res <- minExactMatch(Z ~ E1 + E2 + E3, data = df)

  expect_equal(length(findSubproblems(res)), 3) # uses E1 and partial E2, not E3
  expect_true(all(table(res@groups) %in% c(8, 12)))

  # the formula must have both a  left and right side
  expect_error(minExactMatch(~ E1 + E2), "Formula")

  setFeasibilityConstants() # reset the values to make sure that other tests pass
})

test_that("find size of caliper result", {
  
  # start with the helper function that computes the maximum number of comparisons
  scores <- c(1:5, seq(6, 22, by = 2)) 
  z <- rep(c(1,0), 7)

  # shuffle them so they are not in any useful order prior
  rndorder <- sample(1:14)
  scores <- scores[rndorder]
  z <- z[rndorder]

  # treated: controls within caliper
  # 1: 2
  # 3: 2, 4
  # 5: 4,6
  # 8: 6,10
  # 12: 10, 14
  # 16: 14, 18
  # 20: 18, 22
  # => Total: 13
  expect_equal(caliperSize(scores, z, 2), 13)
  
  # treated: controls within caliper
  # 1: 2, 4
  # 3: 2, 4, 6
  # 5: 2, 4, 6
  # 8: 6, 10
  # 12: 10, 14
  # 16: 14, 18
  # 20: 18, 22
  # => Total: 16
  expect_equal(caliperSize(scores, z, 3), 16)
  
  # include every one! (7 * 7 = 49)
  expect_equal(caliperSize(scores, z, 100), 49)

  # include no one
  expect_error(caliperSize(scores, z, 0), "Invalid caliper width")
  
  # a quicker upper bound test of the caliper size.
  # goal of tests: the upper bound should be less than the max and more than the true value

  res <- caliperUpperBound(scores, z, 3)
  expect_true(res >= 16 & res <= 49)

  res <- caliperUpperBound(scores, z, 2)
  expect_true(res >= 13 & res <= 49)

  # include every one! (7 * 7 = 49)
  expect_equal(caliperUpperBound(scores, z, 100), 49)

  # include no one
  expect_error(caliperUpperBound(scores, z, 0), "Invalid caliper width")

})
