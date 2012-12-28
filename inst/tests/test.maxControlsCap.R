################################################################################
# maxControlsCap: finding the best input for fullmatch
################################################################################

library(testthat)

context("maxControlsCap function")

test_that("basics", {
  # needs to take both simple and stratified problems.
  n <- 16
  Z <- rep(c(0,1), n/2)
  my.names <- paste(rep(c("C", "T"), n/2), 1:16, sep = "")
  names(Z) <- my.names
  B <- c(rep("A", n/2), rep("B", n/2))

  em <- exactMatch(B, treatment = Z) # factor, factor implementation

  res <- maxControlsCap(em)
  
  expect_equal(length(res$strictest.feasible.max.controls), 2) # two level problem
  expect_true(all(!is.na(res$strictest.feasible.max.controls)))
})

test_that("Testing input", {
  # must pass a dist spec
  expect_error(maxControlsCap(1:10),"Invalid distance")

  # if min.controls is a vector, it must be named the same as the name of the subproblems
  Z <- rep(c(0,1), 8)
  B <- rep(letters[1:4], each = 4)
  res.em <- exactMatch(Z ~ B)

  expect_error(maxControlsCap(res.em, min.controls = c(x = 1, y = 2, z = 3, w = 4)),
              "Names of 'min.controls' must match the subproblems. See 'findSubproblems' and 'exactMatch'.")
  # should not expect an error for this:
  maxControlsCap(res.em, min.controls = c(a = 1, b = 1, c = 1, d = 1))

})

test_that("Correct output", {
  # borrowing from the R CMD Check tests to nail down proper output.
  data(nuclearplants, env = parent.env()) 

  mhd2a <- t(match_on(pr ~ date + cum.n, data = nuclearplants) + exactMatch(pr ~ pt, data = nuclearplants))
  res.mxcc <- maxControlsCap(mhd2a + caliper(mhd2a, sqrt(3)))

  expect_equivalent(res.mxcc$strictest.feasible.max.controls, c(1,1))
  
})
