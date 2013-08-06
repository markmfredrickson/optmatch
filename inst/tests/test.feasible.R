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

  expect_equal(length(levels(res)), 3) # uses E1 and partial E2, not E3
  expect_true(all(table(res) %in% c(8, 12)))

  # the formula must have both a  left and right side
  expect_error(minExactMatch(~ E1 + E2), "Formula")

  scores <- rep(1:8,4)
  # minExactMatch can also take a caliper width and a set of scores
  minExactMatch(Z ~ E1 + E2, data = df, scores = scores, width = 1)

  # if you pass one, you must pass both arguments
  expect_error(minExactMatch(Z ~ E1 + E2 + E3, data = df, scores = scores), "width")
  expect_error(minExactMatch(Z ~ E1 + E2 + E3, data = df, width = 1), "scores")

  # the caliper whould allow the problem to be feasible, without using E2
  res <- minExactMatch(Z ~ E1 + E2, data = df, scores = scores, width = 0.5) # very narrow caliper
  expect_equal(length(levels(res)), 2) # goal: only split on E1

  setFeasibilityConstants() # reset the values to make sure that other tests pass

  # don't oversplit: e.g. I(E1 + E2)
  res <- minExactMatch(Z ~ I(E1 + E2) + E3, data = df)
  expect_equal(length(levels(res)), 3)

})

test_that("find size of caliper result", {
  
  # start with the helper function that computes the maximum number of comparisons
  scores <- c(1:5, seq(6, 22, by = 2)) 
  z <- rep(c(1,0), 7)
  b <- rep(c(1,0), each = 7)

  # shuffle them so they are not in any useful order prior
  rndorder <- sample(1:14)
  scores <- scores[rndorder]
  z <- z[rndorder]
  b <- b[rndorder]

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

  ### structure argument defines which units are matchable
  # NB: the groups defined by b are: 1,2,3,4,5,6,8 ; 10,12,14,16,18,20,22
  # so when the caliper width is 2, there are 12 possible matches
  # when caliper = 3, 15 possible
  # when caliper is very large, 3 * 4 + 4 * 2 = 24`
  expect_equal(sum(caliperSize(scores, z, 2, structure = b)), 12)
  expect_equal(sum(caliperSize(scores, z, 3, structure = b)), 15)
  expect_equal(sum(caliperSize(scores, z, 100, structure = b)), 24)

  # likewise for caliperUpperBound, structure argument
  # however, the structure now suggests that the max is 12 per level
  res <- caliperUpperBound(scores, z, 2)
  expect_true(res >= 12 & res <= 24)

  res <- caliperUpperBound(scores, z, 3)
  expect_true(res >= 15 & res <= 24)

  expect_equal(caliperUpperBound(scores, z, 100, structure = b), 24)

  # minCaliper: finds the caliper first caliiper size (from a left-to-right seq) that
  # will be feasible, ie. use fewer arcs than required

  # a caliper with width 2 has 13 arcs, 3 has 16 and would be too wide
  oldopts <- options("optmatch_max_problem_size" = 15)  
  expect_equal(maxCaliper(scores, z, 5:1), 2)

  # if 2 is missing, pick the next best
  expect_equal(maxCaliper(scores, z, c(5,4,3,1, 0.25)), 1)

  # since 2 isn't included, an error should be generated
  expect_error(maxCaliper(scores, z, 5:3), "caliper size")

  # introduce a structure argument, a factor indicating groups
  # even a very wide caliper is helped by the structure. without b, this would take a caliper of 3 
  expect_equal(maxCaliper(scores, z, 5:1, structure = b), 5)

  # tighten down the problem size to require a smaller caliper
  options("optmatch_max_problem_size" = 10)
  expect_equal(maxCaliper(scores, z, 5:1, structure = b), 4)

  # move this down so that the upper bound for width = 2 is too high (15)
  options("optmatch_max_problem_size" = 14)
  # use the upper bound, rather than the exact computation method, to get a caliper value
  expect_equal(maxCaliper(scores, z, 5:1, exact = FALSE), 1)
  
  # play nice with other tests
  setFeasibilityConstants()
})

test_that("match_on does not allow too large problems (via makedist fn)", {
  X <- rnorm(100)
  Z <- rep(c(0,1), 50)
  B <- rep(c(0,1), each = 50)
  
  # expected behavior:
  # exactMatch should create BlockedISMs of any size, as they can be strung
  # together to form smaller problems.
  # match_on, on the other hand, should give a warning when creating a match that
  # is too large, with a hint to use the within argument
  oldopts <- options(warn = 2, "optmatch_max_problem_size" = 25 * 25 + 1)

  # expect no error/warning
  blocking <- exactMatch(Z ~ B)
  # this should be ok
  match_on(Z ~ X, within = blocking)
  options(warn = 0) # back to normal warning behavior

  # give a warning that suggests the within argument
  expect_warning(match_on(Z ~ X), "within")

  # make the max problem smaller, and the warnining should pop up for blocked
  # problems
  options("optmatch_max_problem_size" = 25 * 25 - 1)
  expect_warning(match_on(Z ~ X))
  expect_warning(match_on(Z ~ X, within = blocking), "within")

  # now turn off the warning via the option optmatch_warn_on_big_problem
  options(warn = 2, "optmatch_warn_on_big_problem" = FALSE)
  match_on(Z ~ X)
  options(warn = 0)

  setFeasibilityConstants() 
  

})

