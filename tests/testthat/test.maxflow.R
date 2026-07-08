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

test_that("drop.isolated.rows exempts unmatchable treated units", {
  isolated <- matrix(c(1, 1, Inf, Inf),
                     nrow = 2, byrow = TRUE,
                     dimnames = list(c("A", "B"), c("c1", "c2")))
  r <- maxflow_feasibility(isolated, drop.isolated.rows = TRUE)
  expect_true(r$feasible)
  expect_equal(r$unmatchable.row.units, "B")
  expect_equal(r$deficient.row.units, character(0))
  expect_equal(r$max.controls.matchable, 1L) # A alone, capped at 1
  expect_equal(r$min.omit.fraction, 0.5)
  expect_equal(r$max.mean.controls, 1) # per required (matchable) row unit
})

test_that("recovery via max flow omits as few controls as possible", {
  old <- options(fullmatch_try_recovery = TRUE)
  on.exit(options(old))

  ## Cost-guided recovery (the pre-#200 heuristic) would solve without
  ## max.controls, giving B one control (cost 5) and A four (cost 4), and
  ## after capping at max.controls conclude only 3 controls are usable
  ## (omit.fraction 0.4). The max flow answer: all of B on c1, c2 and two
  ## of c3..c5 on A, so 4 controls (omit.fraction 0.2).
  d <- matrix(c(1, 1, 1, 1, 1,
                5, 5, Inf, Inf, Inf),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("A", "B"), paste0("c", 1:5)))
  fm <- suppressWarnings(
    fullmatch(d, min.controls = 1, max.controls = 2))
  expect_true(all(subproblemSuccess(fm)))
  expect_equal(unname(attr(fm, "omit.fraction")), 0.2)
  expect_equal(sum(!is.na(fm[paste0("c", 1:5)])), 4)
  ## B's matched set keeps both of its permissible controls
  expect_equal(fm[["B"]], fm[["c1"]])
  expect_equal(fm[["B"]], fm[["c2"]])
})

test_that("recovery reports why omission cannot help, when verbose", {
  old <- options(fullmatch_try_recovery = TRUE,
                 optmatch_verbose_messaging = TRUE)
  on.exit(options(old))

  ## A and B each need 2 controls; B can only ever reach c1
  contested <- matrix(c(1, 1, Inf,
                        1, Inf, Inf),
                      nrow = 2, byrow = TRUE,
                      dimnames = list(c("A", "B"), c("c1", "c2", "c3")))
  w <- capture_warnings(
    fm <- fullmatch(contested, min.controls = 2, max.controls = 2))
  expect_true(all(is.na(fm)))
  expect_true(any(grepl("cannot be given min.controls", w)))

  ## a treated unit with no eligible controls at all is named, alongside
  ## the deficiency that actually makes the subproblem unrecoverable
  isolated <- matrix(c(1, 1,
                       1, Inf,
                       Inf, Inf),
                     nrow = 3, byrow = TRUE,
                     dimnames = list(c("A", "B", "C"), c("c1", "c2")))
  w <- capture_warnings(
    fm <- fullmatch(isolated, min.controls = 2, max.controls = 2))
  expect_true(any(grepl("no permissible control: C", w)))
  expect_true(any(grepl("cannot be given min.controls", w)))
})

test_that("#226: pairmatch says plainly when controls demand is impossible", {
  m <- matrix(1, nrow = 2, ncol = 3,
              dimnames = list(c("t1", "t2"), c("c1", "c2", "c3")))
  expect_error(pairmatch(m, controls = 2),
               "not enough controls in some subclasses")
  expect_error(pairmatch(m, controls = 2),
               "needing 4 controls \\(controls = 2\\), but only 3 eligible")
  ## controls = 1 with more treated than controls remains allowed (#116)
  expect_error(suppressWarnings(pairmatch(t(m))), NA)
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
