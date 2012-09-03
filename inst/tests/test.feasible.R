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


