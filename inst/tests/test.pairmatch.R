################################################################################
# Pairmatch tests
################################################################################

library(testthat)
library(optmatch)

context("pairmatch function")

test_that("No cross strata matches", {
  # test data
  Z <- rep(c(0,1), 4)
  B <- rep(c(0,1), each = 4)
  distances <- 1 + exactMatch(Z ~ B)

  res <- pairmatch(distances)
  expect_false(any(is.na(res)))
  expect_false(any(res[1:4] %in% res[5:8]))

})

test_that("Remove unmatchables", {
  A <- matrix(c(1,1,Inf,1,1,Inf,1,1,Inf,1,1,Inf), nrow = 3)
  dimnames(A) <- list(1:3, 4:7)

  expect_true(all(is.na(pairmatch(A, remove.unmatchables = F))))
  expect_true(!all(is.na(pairmatch(A, remove.unmatchables = T))))

  Ai <- as.InfinitySparseMatrix(A)
  expect_true(all(is.na(pairmatch(Ai, remove.unmatchables = F))))
  expect_true(!all(is.na(pairmatch(Ai, remove.unmatchables = T))))


  
})

