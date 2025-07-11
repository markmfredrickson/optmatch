################################################################################
### Tests for objects implementing the distance specification protocol
################################################################################

context("Distance Specification Protocol")

test_that("Matrix => nodes and arcs", {
  # tests:
  expect_true(isGeneric("prepareMatching"))

  # matrix and matrix.csr implement prepareMatching
  # each returns the proper number of nodes and arcs

  # test data: 4 arcs (2 pairs unmatchable)
  m <- matrix(c(1, Inf, 1, 2, 2, Inf), nrow = 2, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E")

  m.result <- prepareMatching(m)
  expect_equal(dim(m.result), c(4, 3))
  expect_equal(unique(m.result$treatment), as.factor(c("D", "E")))
  expect_equal(unique(m.result$control), as.factor(c("A", "B", "C")))

  # NA's and NaN's the same as Inf's
  mm  <- m
  mm[is.infinite(m)]  <- c(NA_real_, NaN)
  expect_equal(prepareMatching(mm), m.result)
})

test_that("ISM => nodes and arcs", {
  m <- matrix(c(1, Inf, 1, 2, 2, Inf), nrow = 2, ncol = 3)
  A.nonames <- as.InfinitySparseMatrix(m)

  expect_error(prepareMatching(A.nonames))

  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E")
  A <- as.InfinitySparseMatrix(m)

  res.ISM <- prepareMatching(A)

  expect_equal(dim(res.ISM), c(4, 3))
  expect_equal(unique(res.ISM$treatment), as.factor(c("D", "E")))
  expect_equal(unique(res.ISM$control), as.factor(c("A", "B", "C")))

  # NA's and NaN's the same as Inf's
  mm  <- m
  mm[is.infinite(m)]  <- c(NA_real_, NaN)
  expect_equal(prepareMatching(as.InfinitySparseMatrix(mm)), res.ISM)

})

test_that("Subproblems", {
  m <- matrix(1, nrow = 2, ncol = 2)
  expect_false(subproblems(m))
  expect_false(subproblems(as.InfinitySparseMatrix(m)))

  B <- rep(c(0,1), each = 5)
  names(B) <- letters[1:10]
  em <- exactMatch(B, rep(c(0,1), 5))
  res.em <- subproblems(em)
  expect_equal(length(res.em), 2)

  expect_true(all(vapply(res.em,
                         validDistanceSpecification,
                         logical(1))))

  m1 <- matrix(0, nrow = 2, ncol = 3,
    dimnames = list(treatment = c("b", "d"), control = c("a", "c", "e")))

  m2 <- matrix(0, nrow = 3, ncol = 2,
    dimnames = list(treatment = c("f", "h", "j"), control = c("g", "i")))

  expect_equal(as.matrix(res.em[[1]]), m1)
  expect_equal(as.matrix(res.em[[2]]), m2)

  # findSubproblems should always return a list, of length # of probs
  expect_equal(length(findSubproblems(m)), 1)
  expect_equal(length(findSubproblems(as.InfinitySparseMatrix(m))), 1)

  B2 <- rep(1:5, each = 2)
  names(B2) <- letters[1:10]
  em2 <- exactMatch(B2, rep(c(0,1), 5))
  em.subps <- findSubproblems(em2)
  expect_equal(length(em.subps), 5)
  expect_is(em.subps, "list")
  lapply(em.subps, function(i) {
    expect_true(validDistanceSpecification(i))
  })


  # a dense problem that looks likes a sparse problem
  position <- rep(1:4, each = 4)
  z <- rep(0:1, 8)
  names(z) <- letters[1:16]
  dist <- match_on(z ~ position, inv.scale.matrix = diag(1))
  allin <- exactMatch(rep(1, 16), z)

  res.allin <- findSubproblems(dist + allin)
  expect_equal(length(res.allin), 1)
  expect_equal(dim(res.allin[[1]]), c(8,8))
  expect_equal(res.allin[[1]]@.Data, as.vector(dist))

  # factors that create unworkable subproblems (e.g. only control units)
  # use z from above
  bad.factor <- c(rep(1:4, each = 2), rep(c(4,5), 4))
  bad.match <- exactMatch(bad.factor, z)
  expect_equal(dim(bad.match), c(8,8))

  # don't return invalid subproblems, in this case, only return 1 subprob
  res.sp <- findSubproblems(bad.match)
  expect_equal(length(res.sp), 4) # there are 5 levels in the factor, but only 4 valid subprobs
})

test_that("Validating DistSpecs", {
  # nameless matrices are not valid:
  m <- matrix(1:4, nrow = 2)

  expect_error(validDistanceSpecification(m))
  expect_error(validDistanceSpecification(as.InfinitySparseMatrix(m)))

  dimnames(m) <- list(1:2, 3:4)
  expect_true(validDistanceSpecification(m))
  expect_true(validDistanceSpecification(as.InfinitySparseMatrix(m)))

  # matrices/isms must be numeric
  m2 <- matrix(letters[1:4], nrow = 2)
  dimnames(m2) <- list(1:2, 3:4)
  expect_equal(mode(m2), "character")

  # while character based ISMs are ok, they are not valid distance specifications
  expect_error(validDistanceSpecification(m2))
  expect_error(validDistanceSpecification(as.InfinitySparseMatrix(m2)))

  # all Inf is valid
  m3 <- matrix(Inf, nrow = 3, ncol = 4)
  rownames(m3) <- LETTERS[1:3]
  colnames(m3) <- letters[23:26]

  expect_true(validDistanceSpecification(m3))
  expect_true(validDistanceSpecification(as.InfinitySparseMatrix(m3)))

  # negative distances not OK
  m4  <- matrix(1:4 -2, nrow = 2)
  dimnames(m4) <- list(1:2, 3:4)
  expect_error(validDistanceSpecification(m4), "can't be negative")
  expect_error(validDistanceSpecification(as.InfinitySparseMatrix(m4)),
               "can't be negative")
  m5  <- matrix(1:4 -1, nrow = 2)
  dimnames(m5) <- list(1:2, 3:4)
  expect_true(validDistanceSpecification(m5))
  expect_true(validDistanceSpecification(as.InfinitySparseMatrix(m5)))
})
