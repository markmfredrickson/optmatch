################################################################################
# Tests for exactMatch function: a function to create InfinitySpareMatrices
################################################################################

context("exactMatch function")

test_that("Exact Match on Factors", {
  n <- 16
  Z <- rep(c(0,1), each = n/2)
  my.names <- c(LETTERS[1:(n/2)], letters[(26 - n/2 + 1):26])
  names(Z) <- my.names

  W <- rnorm(16)
  B <- rep(c(0,1), n/2)
  test.data <- data.frame(Z, W, B)

  res <- exactMatch(B, treatment = Z) # factor, factor implementation

  # the resulting matrix should be block diagonal with 32 non-inf entries

  expect_equal(dim(res), c(8,8))
  expect_equal(length(res), 32)

  expect_error(exactMatch(B, rep(1:(n/4), 4)))
  expect_error(exactMatch(B, c(Z, 0)))
  expect_error(exactMatch(c(B, 1), Z))

  # row and column names
  expect_equal(rownames(res), my.names[Z == 1])
  expect_equal(colnames(res), my.names[Z == 0])
})

test_that("Exact match on formula", {
  n <- 16
  Z <- rep(c(0,1), n/2)
  my.names <- paste(rep(c("C", "T"), n/2), 1:16, sep = "")
  names(Z) <- my.names

  W <- rnorm(16)
  B <- c(rep(0, n/2), rep(1, n/2))
  test.data <- data.frame(Z, W, B)

  res <- exactMatch(Z ~ B)

  # the resulting matrix should be block diagonal
  m0 <- matrix(0, nrow = n/4, ncol = n/4)
  mInf <- matrix(Inf, nrow = n/4, ncol = n/4)

  tmp1 <- cbind(m0, mInf)
  tmp2 <- cbind(mInf, m0)
  m <- rbind(tmp1, tmp2)

  expect_equivalent(as.matrix(res), m)
  expect_equal(dim(res), c(8,8))

  res.data <- exactMatch(Z ~ B, data = test.data)
  expect_equivalent(res.data, res)

  # combine mulitiple factors into a single factor
  B2 <- rep(c(0,1), 4, each = 2)

  # combine them by hand into a single factor
  BB <- B + 2 * B2
  res.bb <- exactMatch(BB, Z)

  res.multi <- exactMatch(Z ~ B + B2)

  expect_equal(as.matrix(res.bb), as.matrix(res.multi))

})

test_that("Use proper environment or data.frame", {
  n <- 16
  Z <- rep(c(0,1), n/2)
  W <- rnorm(16)
  B <- c(rep(0, n/2), rep(1, n/2))
  test.data <- data.frame(a = Z, x = W, c = B)

  names(Z) <- letters[1:n]
  rownames(test.data) <- letters[1:n]

  res.envir <- exactMatch(Z ~ B)
  res.df <- exactMatch(a ~ c, data = test.data)

  expect_equivalent(res.envir, res.df)

})

test_that("Makes correct mask", {
  # this data gave me problems with a makedist() test.
  # it should produces a matrix with a 2x3 0 matrix in
  # the upper left and a 3x2 0 m matrix in the lower right
  # it was producing a 3x3 and a 2x2 for some reason.

  set.seed(20110629)
  data <- data.frame(z = rep(c(1,0), 5),
                     y = rnorm(10),
                     b = rep(c(1,0), each = 5))
  rownames(data) <- letters[1:10]
  Y <- data$z
  A <- data$b
  names(Y) <- rownames(data)
  names(A) <- rownames(data)

  reference <- matrix(c(rep(c(0,0,0,Inf,Inf), 2),
                        rep(c(Inf, Inf, Inf, 0, 0), 3)),
                      nrow = 5, ncol = 5)

  mask.df <- exactMatch(z ~ b, data = data)
  expect_equal(length(mask.df), 3*2 + 2*3) # sizes of the 0 blocks

  mask.fac <- exactMatch(A, Y)
  expect_equal(length(mask.fac), 12)

  expect_equivalent(mask.df, mask.fac)

})

test_that("Must have names", {
  expect_error(exactMatch(rep(c(0,1), each = 5), rep(c(0,1), 5)))
  Z <- rep(c(0,1), 8)
  B <- rep(1:4, each = 4)
  names(B) <- letters[1:6]
  em <- exactMatch(B, Z)

  expect_false(is.null(em@colnames))
  expect_false(is.null(em@rownames))
  expect_false(is.null(names(em@groups)))

  position <- rep(1:4, each = 4)
  z <- rep(0:1, 8)
  names(z) <- letters[1:16]
  dist <- match_on(z ~ position, inv.scale.matrix = diag(1))
  allin <- exactMatch(rep(1, 16), z)

  expect_equal(names(allin@groups), letters[1:16])
})


test_that("Contains grouping information", {
  d <- data.frame(Z = rep(c(0,1), 8),
                  B = rep(letters[1:4], each = 4))

  res.em <- exactMatch(Z ~ B, data=d)
  expect_is(res.em, "BlockedInfinitySparseMatrix")

  # the grouping factor must have names
  expect_equal(length(names(res.em@groups)), 16)
  # ... and those names should match the dimnames of the BISM
  expect_setequal(names(res.em@groups), unlist(dimnames(res.em)))

  # the names of the strata should be used as names of the subprobs list
  expect_equal(names(findSubproblems(res.em)), letters[1:4])

  ### these next few tests are related to eM(), so I'm putting the test here,
  ### but it is implemented in fullmatch.R
  # the result of the fullmatch should use the original names
  fm <- fullmatch(res.em, data=d)
  expect_true(all(1:16 %in% names(fm)))

  # the prefixes shoudl be used in the levels of the factor
  expect_true(all(fm %in% apply(expand.grid(letters[1:4], 1:4), 1, function(r) { paste(r, collapse = ".") })))
})

test_that("t() maintains stratification", {
  Z <- rep(c(0,1), 8)
  B <- rep(letters[1:4], each = 4)

  em <- exactMatch(Z ~ B)
  em.t <- t(em)

  expect_equal(length(findSubproblems(em)), 4)
  expect_equal(length(findSubproblems(em.t)), 4)
})

test_that("Cbind/rbind an exact match", {
  n <- 16
  Z <- rep(c(0,1), each = n/2)
  my.names <- c(LETTERS[1:(n/2)], letters[(26 - n/2 + 1):26])
  names(Z) <- my.names

  W <- rnorm(16)
  B <- rep(c(0,1), n/2)
  test.data <- data.frame(Z, W, B)

  res <- exactMatch(B, treatment = Z) # factor, factor implementation

  mc <- matrix(c(rep(1, n/2), rep(2, n/2)), ncol = 2,
    dimnames = list(letters[(26 - n/2 + 1):26], c("new.1", "new.2")))

  res.cbind <- cbind(res, mc)
  expect_equal(dim(res.cbind), c(n/2, n/2 + 2))

  mr <- t(mc)
  colnames(mr) <- LETTERS[1:(n/2)]
  res.rbind <- rbind(res, mr)
  expect_equal(dim(res.rbind), c(n/2 + 2, n/2))

})

test_that("exactMatch objs can be update()'d", {
  Z <- rep(c(0,1), 8)
  B <- rep(letters[1:4], each = 4)

  simple <- exactMatch(Z ~ B)
  expect_equal(length(levels(simple@groups)), 4)

  B <- rep(letters[1:2], each = 8)
  updated <- update(simple)
  expect_equal(length(levels(updated@groups)), 2)
})


test_that("antiExactMatch", {
  x <- as.factor(c(1,1,2,2,3,3))
  z <- c(0,1,0,1,0,1)
  names(x) <- paste0("X", 1:6)

  ex <- matrix(c(Inf, 0, 0, 0, Inf, 0, 0, 0, Inf), nrow = 3, ncol = 3,
               dimnames =list(
                   treatment = c("X2", "X4", "X6"),
                   control = c("X1", "X3", "X5")))

  res <- antiExactMatch(x, z)

  expect_equal(as.matrix(res), ex)

})

test_that("#123: exactmatch accepts NA treatment", {
  data <- data.frame(z = rep(0:1, each = 5),
                     b = rep(0:1, times = 5))

  m <- match_on(z ~ b, data = data)
  e <- exactMatch(z ~ b, data = data)
  expect_equal(dim(m), dim(e))
  expect_equal(length(e@groups), sum(dim(m)))
  expect_equal(rownames(m), rownames(e))
  expect_equal(colnames(m), colnames(e))
  expect_setequal(names(e@groups), unlist(dimnames(m)))

  data$z[1] <- NA

  m <- match_on(z ~ b, data = data)
  e <- exactMatch(z ~ b, data = data)
  expect_equal(dim(m), dim(e))
  expect_equal(length(e@groups), sum(dim(m)))
  expect_equal(rownames(m), rownames(e))
  expect_equal(colnames(m), colnames(e))
  expect_setequal(names(e@groups), unlist(dimnames(m)))

  data$z[c(2,4,6,7)] <- NA

  m <- match_on(z ~ b, data = data)
  e <- exactMatch(z ~ b, data = data)
  expect_equal(dim(m), dim(e))
  expect_equal(length(e@groups), sum(dim(m)))  
  expect_equal(rownames(m), rownames(e))
  expect_equal(colnames(m), colnames(e))
  expect_setequal(names(e@groups), unlist(dimnames(m)))
})

test_that("#149: exactMatch fails on unique RHS values", {
  t <- rep(0:1, each = 3)
  x <- rnorm(6)
  names(t) <- names(x) <- letters[1:6]

  expect_error(exactMatch(x, t), "no overlap")
  # if x is factor, let it go
  expect_silent(exactMatch(as.factor(x), t))
  
  x <- c(1, 1, 2, 3, 4, 4)
  names(t) <- names(x) <- letters[1:6]
  expect_error(exactMatch(x, t), "no overlap")
  # if x is factor, again let it go
  expect_silent(exactMatch(as.factor(x), t))
})

test_that("#206: maintain dimension if x has NAs", {
  data <- data.frame(z = rep(0:1, each = 5),
                     b = rep(0:1, times = 5))

  m <- match_on(z ~ b, data = data)
  e <- exactMatch(z ~ b, data = data)
  expect_equal(dim(m), dim(e))
  expect_equal(length(e@groups), sum(dim(m)))  
  expect_equal(rownames(m), rownames(e))
  expect_equal(colnames(m), colnames(e))
  expect_setequal(names(e@groups), unlist(dimnames(m)))

  a <- antiExactMatch(setNames(data$b, rownames(data)), data$z)
  expect_equal(dim(m), dim(a))
  expect_equal(rownames(m), rownames(a))
  expect_equal(colnames(m), colnames(a))
  

  data$b[1] <- NA

  e <- exactMatch(z ~ b, data = data)
  expect_equal(dim(m), dim(e))
  expect_equal(length(e@groups), sum(dim(m)))  
  expect_equal(rownames(m), rownames(e))
  expect_equal(colnames(m), colnames(e))
  expect_setequal(names(e@groups), unlist(dimnames(m)))

  a <- antiExactMatch(setNames(data$b, rownames(data)), data$z)
  expect_equal(dim(m), dim(a))
  expect_equal(rownames(m), rownames(a))
  expect_equal(colnames(m), colnames(a))

  data$b[c(2,4,6,7)] <- NA

  e <- exactMatch(z ~ b, data = data)
  expect_equal(dim(m), dim(e))
  expect_equal(length(e@groups), sum(dim(m)))  
  expect_equal(rownames(m), rownames(e))
  expect_equal(colnames(m), colnames(e))
  expect_setequal(names(e@groups), unlist(dimnames(m)))

  a <- antiExactMatch(setNames(data$b, rownames(data)), data$z)
  expect_equal(dim(m), dim(a))
  expect_equal(rownames(m), rownames(a))
  expect_equal(colnames(m), colnames(a))  
})
