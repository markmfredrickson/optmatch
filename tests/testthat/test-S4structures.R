################################################################################
# Tests for new data structures
################################################################################
context("problem and node data structure checks")

test_that("node.data: Single subproblem with double precision", {
  data(nuclearplants)
  f <- fullmatch(pr ~ cost, data=nuclearplants)
  nd <- attr(f, "node.data")
  expect_true(nrow(nd) == nrow(nuclearplants) + 2)
  expect_true(ncol(nd) == 4)
  expect_true(length(unique(nd$group)) == 1)
})


test_that("prob.data: Single subproblem with double precision", {
  data(nuclearplants)
  f <- fullmatch(pr ~ cost, data=nuclearplants)
  pd <- attr(f, "prob.data")
  expect_true(nrow(pd) == 1)
  expect_true(ncol(pd) == 8)
  expect_true(length(unique(pd$group)) == 1)
  expect_true(all(!is.na(pd$reso)))
  expect_true(all(!is.na(pd$tol)))
  expect_true(all((pd$exceedance) > 0))
})


test_that("node.data: two subproblems with double precision", {
  data(nuclearplants)
  f <- fullmatch(pr ~ cost, within=exactMatch(pr ~ pt, data=nuclearplants), data=nuclearplants)
  nd <- attr(f, "node.data")
  expect_true(nrow(nd) == nrow(nuclearplants) + 4)
  expect_true(ncol(nd) == 4)
  expect_true(length(unique(nd$group)) == 2)
})

test_that("prob.data: two subproblems with double precision", {
  data(nuclearplants)
  f <- fullmatch(pr ~ cost, within=exactMatch(pr ~ pt, data=nuclearplants), data=nuclearplants)
  pd <- attr(f, "prob.data")
  expect_true(nrow(pd) == 2)
  expect_true(ncol(pd) == 8)
  expect_true(length(unique(pd$group)) == 2)
  expect_true(all(!is.na(pd$reso)))
  expect_true(all(!is.na(pd$tol)))
  expect_true(all((pd$exceedance) > 0))
})


test_that("prob.data: two subproblems with integer precision", {
  data(nuclearplants)
  np <- nuclearplants
  np$cost <- round(np$cost)
  mfx <- round(match_on(pr ~ cost, within=exactMatch(pr ~ pt, data=np), data=np) * 500)
  f <- fullmatch(mfx, data=np)
  pd <- attr(f, "prob.data")
  expect_true(nrow(pd) == 2)
  expect_true(ncol(pd) == 8)
  expect_true(length(unique(pd$group)) == 2)
  expect_true(all(is.na(pd$reso)))
  expect_true(all(is.na(pd$tol)))
  expect_true(all((pd$exceedance) == 0))
})

test_that("node.data: Single subproblem with integer precision", {
  data(nuclearplants)
  np <- nuclearplants
  np$cost <- round(np$cost)
  mfx <- round(match_on(pr ~ cost, data=np) * 500)
  f <- fullmatch(mfx, data=np)
  nd <- attr(f, "node.data")
  expect_true(nrow(nd) == nrow(nuclearplants) + 2)
  expect_true(ncol(nd) == 4)
  expect_true(length(unique(nd$group)) == 1)
})


test_that("prob.data: Single subproblem with integer precision", {
  data(nuclearplants)
  np <- nuclearplants
  np$cost <- round(np$cost)
  mfx <- round(match_on(pr ~ cost, data=np) * 500)
  f <- fullmatch(mfx, data=np)
  pd <- attr(f, "prob.data")
  expect_true(nrow(pd) == 1)
  expect_true(ncol(pd) == 8)
  expect_true(length(unique(pd$group)) == 1)
  expect_true(all(is.na(pd$reso)))
  expect_true(all(is.na(pd$tol)))
  expect_true(all((pd$exceedance) == 0))
})


test_that("node.data: two subproblems with integer precision", {
  data(nuclearplants)
  np <- nuclearplants
  np$cost <- round(np$cost)
  mfx <- round(match_on(pr ~ cost, within=exactMatch(pr ~ pt, data=np), data=np) * 500)
  f <- fullmatch(mfx, data=np)
  nd <- attr(f, "node.data")
  expect_true(nrow(nd) == nrow(nuclearplants) + 4)
  expect_true(ncol(nd) == 4)
  expect_true(length(unique(nd$group)) == 2)
})


test_that("prob.data: two subproblems with integer precision", {
  data(nuclearplants)
  np <- nuclearplants
  np$cost <- round(np$cost)
  mfx <- round(match_on(pr ~ cost, within=exactMatch(pr ~ pt, data=np), data=np) * 500)
  f <- fullmatch(mfx, data=np)
  pd <- attr(f, "prob.data")
  expect_true(nrow(pd) == 2)
  expect_true(ncol(pd) == 8)
  expect_true(length(unique(pd$group)) == 2)
  expect_true(all(is.na(pd$reso))) # NA for integer distance
  expect_true(all(is.na(pd$tol)))
  expect_true(all((pd$exceedance) == 0)) #exceedance should be zero for integer data
})
