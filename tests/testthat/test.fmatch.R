################################################################################
### R/Fortran Interface Tests
################################################################################

context("R/Fortran Interface")

if (requireNamespace("rrelaxiv", quietly = TRUE)) {
  slvr <- "RELAX-IV"
} else {
  slvr <- "LEMON"
}


test_that("fmatch accepts DistanceSpecifications", {
  v <- c(1, Inf, 2,
         2, 1, Inf,
         3, 2, 1)

  # and doesn't accept other things...
  expect_error(fmatch(v, 2, 2, solver = slvr))

  # the goal of this matrix is that there is a clear match to make
  # A:D, B:E, C:F
  m <- matrix(v, nrow = 3, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E", "F")
  pm <- prepareMatching(m)

  res <- fmatch(pm, 2, 2, solver = slvr)

  expect_equal(dim(res), c(7,4)) # seven non-Inf entries

  # check that A-D is a pair and A-B is not a match
  expect_equal(res[res$control == "A" & res$treated == "D", "solution"], 1)
  expect_equal(res[res$control == "A" & res$treated == "B",
    "solution"], numeric(0))

  M <- as.InfinitySparseMatrix(m)
  pM <- prepareMatching(M)
  res.ism <- fmatch(pM, 2, 2, solver = slvr)
  expect_identical(res$solution, res.ism$solution)
})
#}

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

  # handles failed matchings by returning NULL
  noMatches <- cbind(skeleton, solution = -1)

  expect_true(is.null(solution2factor(noMatches)))
})

test_that("LEMON solvers", {
  v <- c(1, Inf, 2,
         2, 1, Inf,
         3, 2, 1)
  m <- matrix(v, nrow = 3, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E", "F")
  pm <- prepareMatching(m)

  expect_error(fmatch(pm, 2, 2))


  f_lemon <- fmatch(pm, 2, 2, solver = slvr)
  f_cycle <- fmatch(pm, 2, 2, solver = LEMON("CycleCancelling"))
  f_capac <- fmatch(pm, 2, 2, solver = LEMON("CapacityScaling"))
  f_costs <- fmatch(pm, 2, 2, solver = LEMON("CostScaling"))
  f_netwo <- fmatch(pm, 2, 2, solver = LEMON("NetworkSimplex"))

  expect_identical(f_lemon, f_cycle)
  expect_identical(f_lemon, f_capac)
  expect_identical(f_lemon, f_costs)
  expect_identical(f_lemon, f_netwo)

  if (requireNamespace("rrelaxiv", quietly = TRUE)) {
    f_relax <- fmatch(pm, 2, 2, solver = "RELAX-IV")
    expect_identical(f_relax, f_lemon)
  }

})
