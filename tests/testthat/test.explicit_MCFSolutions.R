################################################################################
# Testing some experimental work on returning feasibility explicitly via MCFSolutions
################################################################################

context("Test Explicit MCFSolutions")

test_that("simple case, all feasible", {
  d <- data.frame(position = rep(1:4, each = 4),
                  z = rep(0:1, 8),
                  rownames=letters[1:16])
  dist <- match_on(z ~ position, inv.scale.matrix = diag(1), data=d)
  res.mat <- fullmatch(dist, data=d)
  expect_true(identical((get_subproblem_info(res.mat, type = "feasible")), TRUE))


  data(nuclearplants)
  f2 <- fullmatch(pr ~ cost + strata(pt), data=nuclearplants)
  expect_true(identical((get_subproblem_info(f2, type = "feasible")), c(TRUE, TRUE)))
  dt <- get_subproblem_info(f2, type = c("resolution", "group"))
  expect_true(all(dim(dt) == c(2,2)))

})


test_that("Feasibility with subproblems + recovery",{
  data(nuclearplants)

  mm <- match_on(pr ~ cost + t1 + t2, data=nuclearplants)
  # infeasible as given
  options("fullmatch_try_recovery" = FALSE)
  expect_warning(f <- fullmatch(mm, data=nuclearplants, max.controls=2))
  expect_true(identical((get_subproblem_info(f, type = "feasible")), c(FALSE)))
  options("fullmatch_try_recovery" = TRUE)
  f <- fullmatch(mm, data=nuclearplants, max.controls=2)
  expect_true(identical((get_subproblem_info(f, type = "feasible")), c(TRUE)))


  # not infeasible as mean.controls provided
  f <- fullmatch(mm, data=nuclearplants, max.controls=2, mean.controls=2)
  expect_true(identical((get_subproblem_info(f, type = "feasible")), c(TRUE)))


  options("fullmatch_try_recovery" = FALSE)
  mm <- match_on(pr ~ cost + t1 + t2, data=nuclearplants, within=exactMatch(pr ~ pt, data=nuclearplants))
  # infeasible as given for subproblem 1, feasible for subproblem 2
  expect_warning(f <- fullmatch(mm, data=nuclearplants, max.controls=2))
  expect_true(identical((get_subproblem_info(f, type = "feasible")), c(FALSE, TRUE)))

  options("fullmatch_try_recovery" = TRUE)
  f <- fullmatch(mm, data=nuclearplants, max.controls=2)
  expect_true(identical((get_subproblem_info(f, type = "feasible")), c(TRUE, TRUE)))

})


test_that("completely empty problem", {
  m <- matrix(Inf, nrow = 3, ncol = 4)
  rownames(m) <- LETTERS[1:3]
  colnames(m) <- letters[23:26]
  expect_warning(res.m <- fullmatch(m))
  expect_true(get_subproblem_info(res.m, type = "feasible") == FALSE)

})
