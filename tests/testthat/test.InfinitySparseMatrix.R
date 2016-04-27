################################################################################
### Tests for the InfinitySparseMatrix class
###############################################################################

context("InfinitySparseMatrix tests")

test_that("ISM Basics", {
  A <- makeInfinitySparseMatrix(c(1,2,3), cols = c(1,2, 2), rows = c(1,1,2))
  expect_is(A, "InfinitySparseMatrix")
  expect_equal(dim(A), c(2,2))

  # converting to the equivalent matrix
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2)
  expect_equivalent(as.matrix(A), m)

  # converting from a matrix to a ISM
  expect_equivalent(as.InfinitySparseMatrix(m), A)
  # and back again
  expect_equivalent(as.matrix(as.InfinitySparseMatrix(m)), m)
  # expect_equal(as(m, "InfinitySparseMatrix"), A)

  # a more complicated examples, missing an entire row/col
  w <- matrix(c(1,Inf,2, 3, Inf, 4), nrow = 3)
  B <- as.InfinitySparseMatrix(w)
  expect_equivalent(as.matrix(B), w)

  y <- matrix(c(1,2,3,Inf, Inf, Inf), nrow = 3)
  D <- as.InfinitySparseMatrix(y)
  expect_equivalent(as.matrix(D), y)

  # the as() technique should be equivalent
  expect_equivalent(as(D, "matrix"), y)
  expect_equivalent(A, as(m, "InfinitySparseMatrix"))

})

test_that("ISM Handles Names", {
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2,
              dimnames = list(treated = c("A", "B"),
                              control = c("C", "D")))

  expect_equal(as.matrix(as(m, "InfinitySparseMatrix")), m)

  A <- makeInfinitySparseMatrix(c(1,2,3), rows = c(1,1,2), cols = c(1,2,2))
  expect_true(is.null(dimnames(A)))

  dms <- list(treated = c("A", "B"), control = c("x", "y"))
  dimnames(A) <- dms
  expect_equal(dimnames(A), dms)

  dimnames(m) <- dms
  expect_equal(as.matrix(A), m)

  dimnames(A) <- NULL
  expect_null(dimnames(A))
})

test_that("Math Ops", {
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2)
  A <- as.InfinitySparseMatrix(m)

  # scalar math
  expect_equivalent(as.matrix(A + 1), m + 1)
  expect_equivalent(as.matrix(A - 1), m - 1)
  expect_equivalent(as.matrix(A * 2), m * 2)
  expect_equivalent(as.matrix(A / 2), m / 2)

  # matrix element wise math
  expect_equivalent(as.matrix(A + A), m + m)

  # Inf - Inf or Inf / Inf gives NA
  mm <- m - m
  mm[is.na(mm)] <- Inf

  md <- m / m
  md[is.na(md)] <- Inf

  expect_equivalent(as.matrix(A - A), mm)
  expect_equivalent(as.matrix(A * A), m * m)
  expect_equivalent(as.matrix(A / A), md)

  # The harder case is when the matrix has non-identical row/col ids

  q <- matrix(c(1, 2, Inf, 4), nrow = 2, ncol = 2)
  B <- as.InfinitySparseMatrix(q)

  expect_equivalent(as.matrix(A + B), m + q)
  expect_equivalent(as.matrix(A * B), m * q)

  # TODO, make up temp matrices for sub and div

  # dense + sparse => sparse
  Aq = A + q
  expect_is(Aq, "InfinitySparseMatrix")
  expect_equivalent(as.matrix(Aq), m + q)

  # make sure it works the other direction (and with mult)
  qA = q * A
  expect_is(qA, "InfinitySparseMatrix")
  expect_equivalent(as.matrix(qA), q * m)

  # names should be sticky across arithmatic
  # TODO, math should reorder by names in case that changes things
  colnames(A) <- paste("C", 1:2, sep = "")
  rownames(A) <- paste("T", 1:2, sep = "")
  colnames(q) <- paste("C", 1:2, sep = "")
  rownames(q) <- paste("T", 1:2, sep = "")

  Aq = A + q
  expect_equal(colnames(Aq), c("C1", "C2"))
  expect_equal(rownames(Aq), c("T1", "T2"))

  # math ops over two matrices with same rows/names bu in different orders
  B <- as.InfinitySparseMatrix(q) # q got rownames later
  q.reordered <- q[,2:1]
  C <- as.InfinitySparseMatrix(q.reordered)
  expect_equal(colnames(C), rev(colnames(B)))
  expect_equal(A + C, A  + B)
})

test_that("Math ops with vectors", {

  # Small matrix with manual calculation
  m <- matrix(c(1, 4, 2, 3), nrow = 2, ncol = 2)
  A <- optmatch:::as.InfinitySparseMatrix(m)
  v <- 1:2

  expect_true(all.equal(attributes(A), attributes(A/v)))
  expect_true(all.equal(attributes(A), attributes(A*v)))
  expect_true(all.equal(attributes(A), attributes(A-v)))
  expect_true(all.equal(attributes(A), attributes(A+v)))
  expect_true(all.equal(attributes(A), attributes(A^v)))
  expect_true(all.equal(attributes(A), attributes(A%%v)))
  expect_true(all.equal(attributes(A), attributes(A%/%v)))

  expect_true(all.equal(attributes(A), attributes(v/A)))
  expect_true(all.equal(attributes(A), attributes(v*A)))
  expect_true(all.equal(attributes(A), attributes(v-A)))
  expect_true(all.equal(attributes(A), attributes(v+A)))
  expect_true(all.equal(attributes(A), attributes(v^A)))
  expect_true(all.equal(attributes(A), attributes(v%%A)))
  expect_true(all.equal(attributes(A), attributes(v%/%A)))

  expect_true(all(as.vector(A/v) == c(1,2,2,3/2)))
  expect_true(all(as.vector(A*v) == c(1,8,2,6)))
  expect_true(all(as.vector(A+v) == c(2,6,3,5)))
  expect_true(all(as.vector(A-v) == c(0,2,1,1)))
  expect_true(all(as.vector(A^v) == c(1,16,2,9)))
  expect_true(all(as.vector(A%%v) == c(0,0,0,1)))
  expect_true(all(as.vector(A%/%v) == c(1,2,2,1)))

  expect_true(all(as.vector(v/A) == c(1,1/2, 1/2, 2/3)))
  expect_true(all(as.vector(v*A) == c(1,8,2,6)))
  expect_true(all(as.vector(v+A) == c(2,6,3,5)))
  expect_true(all(as.vector(v-A) == c(0,-2,-1,-1)))
  expect_true(all(as.vector(v^A) == c(1,16,1,8)))
  expect_true(all(as.vector(v%%A) == c(0,2,1,2)))
  expect_true(all(as.vector(v%/%A) == c(1,0,0,0)))

  # Logical operations
  m2 <- m
  m2[1,2] <- Inf
  A2 <- optmatch:::as.InfinitySparseMatrix(m2)
  expect_is(A2 <= c(1,3), "InfinitySparseMatrix")
  expect_equal(as.vector(A2 <= c(1,3)), c(T, F, T))
  expect_equal(as.vector(c(1,3) >= A2), c(T, F, T))

  # BlockedInfinitySparseMatrix
  x <- c(rep(1,4), rep(2,2), rep(3,5))
  set.seed(1)
  y <- runif(11)
  z <- c(0,0,1,0,1,0,1,1,0,0,0)

  A <- match_on(z~y, within=exactMatch(z~x))
  m <- as.matrix(A)
  v <- 1:5

  # There's some dimensionality issues here, so we'll get lots of "not
  # a multiple" warnings.
  expect_warning({expect_true(all(as.matrix(A/v) == m/v))
    expect_true(all(as.matrix(A*v) == m*v))
    expect_true(all(as.matrix(A+v) == m+v))
    expect_true(all(as.matrix(A-v) == m-v))
    expect_true(all(as.matrix(A^v) == m^v))
    expect_true(all(as.matrix(A%%v) == m%%v, na.rm=TRUE))
    expect_true(all(as.matrix(A%/%v) == m%/%v, na.rm=TRUE))
    vm <- v/m
    vm[!is.finite(as.matrix(A))] <- Inf
    expect_true(all(as.matrix(v/A) == vm))
    expect_true(all(as.matrix(v*A) == v*m))
    expect_true(all(as.matrix(v+A) == v+m))
    vm <- v-m
    vm[!is.finite(as.matrix(A))] <- Inf
    expect_true(all(as.matrix(v-A) == vm))
    vm <- v^m
    vm[!is.finite(as.matrix(A))] <- Inf
    expect_true(all(as.matrix(v^A) == vm))
    expect_true(all(as.matrix(v%%A) == v%%m, na.rm=TRUE))
    expect_true(all(as.matrix(v%/%A) == v%/%m, na.rm=TRUE))},
    "not a multiple")


  # Error on non-numeric input
  expect_error("a"*A, "non-numeric")

})

test_that("Subsetting", {
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2)
  rownames(m) <- c("A", "B")
  colnames(m) <- c("C", "D")
  A <- as.InfinitySparseMatrix(m)

  res.sub <- subset(A, c(TRUE, FALSE))
  expect_equal(res.sub@.Data, c(1, 2))
  expect_equal(res.sub@cols, c(1,2))
  expect_equal(res.sub@rows, c(1,1))

})

test_that("cbinding ISMs and matrices", {
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2)
  rownames(m) <- c("A", "B")
  colnames(m) <- c("C", "D")
  A <- as.InfinitySparseMatrix(m)

  # Expect warnings for duplicate column names
  expect_warning(res.AA <- cbind(A, A))
  expect_equal(length(res.AA), 6)
  expect_equal(dim(res.AA), c(2, 4))

  # and the names should be uniquified (that's a word, really!)
  expect_equal(length(unique(colnames(res.AA))), 4)

  # same for matrices
  expect_warning(res.Am <- cbind(A, m))
  expect_equal(res.Am, res.AA)

  # flipped name order shouldn't matter
  m2 <- m
  rownames(m2) <- c("B", "A")
  expect_warning(res.Am2 <- cbind(A, m2))

  m4 <- matrix(1, nrow = 2, ncol = 3)
  rownames(m4) <- c("A", "C")
  colnames(m4) <- c("X", "Y", "Z")
  expect_error(cbind(A, m4))

  m5 <- matrix(1, nrow = 3, ncol = 2)
  rownames(m5) <- c("A", "B", "C")
  colnames(m5) <- c("X", "Y")
  expect_error(cbind(A, m5))

})

test_that("rbinding ISMs and matrices", {
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2)
  rownames(m) <- c("A", "B")
  colnames(m) <- c("C", "D")
  A <- as.InfinitySparseMatrix(m)

  # Expect warnings for duplicate row names
  expect_warning(res.AA <- rbind(A, A))
  expect_equal(length(res.AA), 6)
  expect_equal(dim(res.AA), c(4,2))

  # and the names should be uniquified (that's a word, really!)
  expect_equal(length(unique(rownames(res.AA))), 4)

  expect_warning(res.Am <- rbind(A, m), "share row names")
  expect_equal(res.Am, res.AA)

  # flipped column names should not matter
  m2 <- m
  colnames(m2) <- c("D", "C")
  expect_warning(res.Am2 <- rbind(A, m2))

  m4 <- matrix(1, nrow = 2, ncol = 2)
  rownames(m4) <- c("A", "B")
  colnames(m4) <- c("X", "Y")
  expect_error(rbind(A, m4))

  m5 <- matrix(1, nrow = 2, ncol = 3)
  rownames(m5) <- c("A", "B")
  colnames(m5) <- c("C", "D", "E")

  expect_error(rbind(A, m5))

})


test_that("t(ransform) function", {
  # set up the names on the dims backwards to that when
  # we call t(m), everything is labeled properly
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2,
              dimnames = list(control = c("A", "B"),
                              treated = c("C", "D")))
  A <- as.InfinitySparseMatrix(m)

  expect_equal(as.matrix(t(A)), t(m))

})

################################################################################
# Tests for the BlockedISM subclass
################################################################################

test_that("BlockedISM addition", {
  Z <- rep(c(0,1), 8)
  B1 <- rep(1:4, each = 4)
  B2 <- rep(c(0,1), each = 8)

  res.b1 <- exactMatch(Z ~ B1)
  res.b2 <- exactMatch(Z ~ B2)

  res.b1b1 <- res.b1 + res.b1
  expect_equal(res.b1b1@groups, res.b1@groups)

  # should use the smaller of the two's groups
  res.b2b1 <- res.b2 + res.b1
  expect_equal(res.b2b1@groups, res.b1@groups)

  expect_is(res.b2 + 1, "BlockedInfinitySparseMatrix")
  expect_is(res.b2 + matrix(1, nrow = 8, ncol = 8),
    "BlockedInfinitySparseMatrix")
  expect_is(matrix(1, nrow = 8, ncol = 8) + res.b2,
    "BlockedInfinitySparseMatrix")
})

test_that("Get subproblem size of each block", {
  Z <- rep(c(0,1), 8)
  B1 <- c(rep('a',3),rep('b', 3), rep('c', 6), rep('d', 4))
  B2 <- c(rep(0, 7), rep(1, 9))
  B3 <- c('a', rep('b', 15)) # group a has no treatment.

  res.b1 <- exactMatch(Z ~ B1)
  res.b2 <- exactMatch(Z ~ B2)
  res.b3 <- exactMatch(Z ~ B3)

  expect_equal(subdim(res.b1), list('a' = c(1, 2),'b' = c(2, 1),'c' = c(3, 3),'d' = c(2, 2)))
  expect_equal(subdim(res.b2), list('0' = c(3, 4),'1' = c(5, 4)))
  expect_equal(subdim(res.b3), list('b' = c(8, 7)))

  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2,
              dimnames = list(control = c("A", "B"),
                  treated = c("C", "D")))
  a <- as.InfinitySparseMatrix(m)

  # subdim on a matrix or non-blocked ISM is equivalent to calling dim
  expect_equal(subdim(m), list(dim(m)))
  expect_equal(subdim(a), list(dim(a)))

  # test on optmatch.dlist
  od <- list(matrix(c(0,0,0,0, Inf,Inf,Inf,Inf, rep(0, 12)), nrow = 2, ncol = 10, dimnames = list(letters[1:2], letters[3:12])),
             matrix(c(0,0,0,0, rep(Inf, 16)), byrow = T, nrow = 5, ncol = 4, dimnames = list(letters[10:14], letters[15:18])))

  class(od) <- c("optmatch.dlist", "list")

  expect_equal(subdim(od), list(c(2, 10), c(5, 4)))

  # test on DenseMatrix
  W <- rnorm(16)

  m <- match_on(Z ~ W)

  expect_equal(subdim(m), list(dim(m)))
})


test_that("ISM sorting", {
  X <- makeInfinitySparseMatrix(data = c(6,5,2,3,1),
                                cols = c(2,1,2,1,1),
                                rows = c(3,3,1,2,1))

  # Output should still be ISM
  expect_is(X, "InfinitySparseMatrix")
  expect_is(sort(X), "InfinitySparseMatrix")
  expect_is(sort(X, byCol=TRUE), "InfinitySparseMatrix")

  X.rows <- sort(X, byCol=FALSE)
  X.cols <- sort(X, byCol=TRUE)

  expect_identical(dim(X.cols), dim(X))
  expect_identical(dim(X.rows), dim(X))

  expect_identical(as.matrix(X.cols), as.matrix(X))
  expect_identical(as.matrix(X.rows), as.matrix(X))


  # pairwise coords should be sorted, e.g.
  # (1,1), (1,2), (2,1), (2,2)
  # In original X, this is not true.
  coordrc <- as.numeric(paste(attr(X, "rows"), attr(X, "cols"), sep=""))
  coordcr <- as.numeric(paste(attr(X, "cols"), attr(X, "rows"), sep=""))
  expect_true(is.unsorted(coordrc))
  expect_true(is.unsorted(coordcr))

  # When sorting by column, then when looking at column/row it should be
  # true.
  coordrc.sortcols <- as.numeric(paste(attr(X.cols, "rows"), attr(X.cols, "cols"), sep=""))
  coordcr.sortcols <- as.numeric(paste(attr(X.cols, "cols"), attr(X.cols, "rows"), sep=""))
  expect_true(is.unsorted(coordrc.sortcols))
  expect_true(!is.unsorted(coordcr.sortcols))

  # Ditto when sorting by row & looking at row first.
  coordrc.sortrows <- as.numeric(paste(attr(X.rows, "rows"), attr(X.rows, "cols"), sep=""))
  coordcr.sortrows <- as.numeric(paste(attr(X.rows, "cols"), attr(X.rows, "rows"), sep=""))
  expect_true(!is.unsorted(coordrc.sortrows))
  expect_true(is.unsorted(coordcr.sortrows))

  # Checking for bad input on byCol
  expect_silent(sort(X, byCol=1))
  expect_error(sort(X, byCol="a"))
  expect_warning(expect_error(sort(X, byCol=c(1,1))))

  # Checking argument `decreasing`
  X.rows <- sort(X, byCol=FALSE, decreasing=TRUE)
  X.cols <- sort(X, byCol=TRUE, decreasing=TRUE)

  expect_identical(as.matrix(X.rows), as.matrix(X))
  expect_identical(as.matrix(X.cols), as.matrix(X))

  coordrc.sortcols <- as.numeric(paste(attr(X.cols, "rows"), attr(X.cols, "cols"), sep=""))
  coordcr.sortcols <- as.numeric(paste(attr(X.cols, "cols"), attr(X.cols, "rows"), sep=""))
  expect_true(is.unsorted(coordrc.sortcols))
  # to check sorting, reverse the order.
  expect_true(!is.unsorted(rev(coordcr.sortcols)))

  data(nuclearplants)

  m <- match_on(pr ~ cost, data=nuclearplants, caliper=1)

  m.rows <- sort(m, byCol=FALSE)
  m.cols <- sort(m, byCol=TRUE)

  # by default, ISM's are row dominant, so resorting by row should not
  # have any impact.
  expect_identical(m, m.rows)

  # However, sorting by column should change the internals, but not
  # externals.
  expect_identical(as.matrix(m), as.matrix(m.cols))
  expect_false(identical(m, m.cols))

  # Double-sorting
  expect_identical(m, sort(m.cols))
})


test_that("BISM sorting", {

  b <- makeInfinitySparseMatrix(c(1,2,3,4,5,6),
                                cols=c(1,2,2,3,4,3),
                                rows=c(1,1,2,3,3,4),
                                colnames=c("1", "3", "5", "7"),
                                rownames=c("2", "4", "6", "8"))

  attr(b, "groups") <- factor(rep(c(1,2), each=4))
  names(attr(b, "groups")) <- 1:8
  class(b) <- "BlockedInfinitySparseMatrix"


  # Output should still be BISM
  expect_is(b, "BlockedInfinitySparseMatrix")
  expect_is(sort(b), "BlockedInfinitySparseMatrix")
  expect_is(sort(b, byCol=TRUE), "BlockedInfinitySparseMatrix")

  b.rows <- sort(b, byCol=FALSE)
  b.cols <- sort(b, byCol=TRUE)

  expect_identical(dim(b.cols), dim(b))
  expect_identical(dim(b.rows), dim(b))

  expect_identical(as.matrix(b.cols), as.matrix(b))
  expect_identical(as.matrix(b.rows), as.matrix(b))

  expect_identical(as.matrix(b), as.matrix(sort(b, decreasing=TRUE)))

  # When sorting by column, then when looking at column/row it should be
  # true.
  coordrc.sortcols <- as.numeric(paste(attr(b.cols, "rows"), attr(b.cols, "cols"), sep=""))
  coordcr.sortcols <- as.numeric(paste(attr(b.cols, "cols"), attr(b.cols, "rows"), sep=""))
  expect_true(is.unsorted(coordrc.sortcols))
  expect_true(!is.unsorted(coordcr.sortcols))

  # Ditto when sorting by row & looking at row first.
  coordrc.sortrows <- as.numeric(paste(attr(b.rows, "rows"), attr(b.rows, "cols"), sep=""))
  coordcr.sortrows <- as.numeric(paste(attr(b.rows, "cols"), attr(b.rows, "rows"), sep=""))
  expect_true(!is.unsorted(coordrc.sortrows))
  expect_true(is.unsorted(coordcr.sortrows))

  # Checking for bad input on byCol
  expect_silent(sort(b, byCol=1))
  expect_error(sort(b, byCol="a"))
  expect_warning(expect_error(sort(b, byCol=c(1,1))))

  data(nuclearplants)

  m <- match_on(pr ~ cost, data=nuclearplants,
                within=exactMatch(pr ~ ct, data=nuclearplants))

  m.rows <- sort(m, byCol=FALSE)
  m.cols <- sort(m, byCol=TRUE)

  # by default, ISM's are row dominant, so resorting by row should not
  # have any impact.
  expect_identical(m, m.rows)

  # However, sorting by column should change the internals, but not
  # externals.
  expect_identical(as.matrix(m), as.matrix(m.cols))
  expect_false(identical(m, m.cols))

  # Double-sorting
  expect_identical(m, sort(m.cols))

})
