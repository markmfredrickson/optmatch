################################################################################
### Tests for the distUnion function
###############################################################################

context("distUnion tests")

test_that("distUnion Basics", {

  mab <- matrix(1:4, nrow = 2, ncol = 2, dimnames = list(c("A", "B"), c("Y", "Z")))
  mac <- matrix(5:6, nrow = 2, ncol = 1, dimnames = list(c("A", "C"), c("X")))
  mbc <- matrix(7:8, nrow = 2, ncol = 1, dimnames = list(c("B", "C"), c("W")))
  expected <- matrix(c(Inf, 7, 8, 5, Inf, 6, 1, 2, Inf, 3, 4, Inf),
                     nrow = 3, ncol = 4,
                     dimnames = list(treated = c("A", "B", "C"), control = c("W", "X", "Y", "Z")))

  iab <- optmatch:::as.InfinitySparseMatrix(mab)
  iac <- optmatch:::as.InfinitySparseMatrix(mac)
  ibc <- optmatch:::as.InfinitySparseMatrix(mbc)
  iex <- optmatch:::as.InfinitySparseMatrix(expected)

  res.m <- distUnion(mab, mac, mbc)
  res.i <- distUnion(iab, iac, ibc)

  expect_equal(res.m, res.i)

  standardize <- function(i) {
    tmp <- as.matrix(i)
    tmp <- tmp[c("A", "B", "C"), ]
    tmp <- tmp[, c("W", "X", "Y", "Z")]
    return(tmp)
  }

  expect_equal(standardize(res.m), expected)
  expect_equal(standardize(res.i), expected)

  # uh oh, duplicates!

  mdup <- matrix(99, nrow = 1, ncol = 1, dimnames = list(c("B"), c("Z")))

  res.dup <- distUnion(mab, mac, mbc, mdup)

  expect_equal(length(res.dup), 8)

  expect_equal(standardize(res.dup), expected)
})
