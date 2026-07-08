################################################################################
# Max-flow feasibility tests (issue #200 prototype)
################################################################################

context("maxflow feasibility")

## A can match any control; B can only match c1
simple.dist <- matrix(c(1, 1, 1, 1,
                        1, Inf, Inf, Inf),
                      nrow = 2, byrow = TRUE,
                      dimnames = list(c("A", "B"),
                                      c("c1", "c2", "c3", "c4")))

test_that("pair matching: max usable controls and minimal omit.fraction", {
  r <- maxflow_feasibility(simple.dist, min.controls = 1, max.controls = 1)
  expect_true(r$feasible)
  expect_equal(r$max.controls.matchable, 2L) # B takes c1, A takes one other
  expect_equal(r$min.omit.fraction, 0.5)
  expect_equal(r$max.mean.controls, 1)
  expect_equal(r$deficient.row.units, character(0))
})

test_that("1:k matching: larger max.controls raises the max flow", {
  r <- maxflow_feasibility(simple.dist, min.controls = 1, max.controls = 2)
  expect_true(r$feasible)
  expect_equal(r$max.controls.matchable, 3L) # A: 2 of c2..c4, B: c1
  expect_equal(r$min.omit.fraction, 0.25)
  expect_equal(r$max.mean.controls, 1.5)

  r <- maxflow_feasibility(simple.dist, min.controls = 1, max.controls = 4)
  expect_equal(r$max.controls.matchable, 4L) # all controls usable
  expect_equal(r$min.omit.fraction, 0)
})

test_that("infeasible lower bounds are detected and attributed", {
  ## A and B compete for the single control both can accept
  contested <- matrix(c(1, 1, Inf,
                        1, Inf, Inf),
                      nrow = 2, byrow = TRUE,
                      dimnames = list(c("A", "B"), c("c1", "c2", "c3")))
  r <- maxflow_feasibility(contested, min.controls = 2, max.controls = 2)
  expect_false(r$feasible)
  expect_true(is.na(r$min.omit.fraction))
  expect_true(is.na(r$max.mean.controls))
  ## Max flows are not unique, so attribution can vary: A may or may not
  ## appear, but B (who can reach only c1) is deficient in every max flow
  expect_true("B" %in% r$deficient.row.units)

  ## a row unit with no finite distances at all
  isolated <- matrix(c(1, 1, Inf, Inf),
                     nrow = 2, byrow = TRUE,
                     dimnames = list(c("A", "B"), c("c1", "c2")))
  r <- maxflow_feasibility(isolated)
  expect_false(r$feasible)
  expect_equal(r$deficient.row.units, "B")
})

test_that("input validation", {
  expect_error(maxflow_feasibility(simple.dist, min.controls = 0),
               "min.controls must be at least 1")
  expect_error(maxflow_feasibility(simple.dist, min.controls = 2,
                                   max.controls = 1),
               "min.controls may not exceed max.controls")
})

test_that("fullmatch agrees: min.omit.fraction is exactly feasible", {
  old <- options(fullmatch_try_recovery = FALSE)
  on.exit(options(old))

  r <- maxflow_feasibility(simple.dist, min.controls = 1, max.controls = 1)

  ## at the max-flow-derived omit.fraction the problem solves ...
  f.at <- suppressWarnings(
    fullmatch(simple.dist, min.controls = 1, max.controls = 1,
              omit.fraction = r$min.omit.fraction))
  expect_true(all(subproblemSuccess(f.at)))
  expect_equal(sum(!is.na(f.at)[c("c1", "c2", "c3", "c4")]),
               r$max.controls.matchable)

  ## ... but asking to keep even one more control is infeasible
  n.keep <- r$max.controls.matchable + 1L
  f.past <- suppressWarnings(
    fullmatch(simple.dist, min.controls = 1, max.controls = 1,
              omit.fraction = 1 - n.keep / ncol(simple.dist)))
  expect_false(any(subproblemSuccess(f.past)))
})

test_that("fullmatch agrees when min.controls < max.controls", {
  old <- options(fullmatch_try_recovery = FALSE)
  on.exit(options(old))

  r <- maxflow_feasibility(simple.dist, min.controls = 1, max.controls = 2)

  ## design assumption: the max flow under max.controls caps is jointly
  ## attainable with the min.controls lower bounds
  f.at <- suppressWarnings(
    fullmatch(simple.dist, min.controls = 1, max.controls = 2,
              omit.fraction = r$min.omit.fraction))
  expect_true(all(subproblemSuccess(f.at)))
  expect_equal(sum(!is.na(f.at)[c("c1", "c2", "c3", "c4")]),
               r$max.controls.matchable)
})
