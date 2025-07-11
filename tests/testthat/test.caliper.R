################################################################################
# Caliper Tests
################################################################################

context("Caliper")

test_that("Caliper return values", {
  m <- matrix(c(1,Inf, 2, 3), nrow = 2, ncol = 2,
              dimnames = list(treatment = c("A", "B"),
                              control = c("C", "D")))
  A <- as.InfinitySparseMatrix(m)

  # use the Mahalanobis distance match_on method
  result <- caliper(A, 2)
  expect_true(validDistanceSpecification(result))

  expect_equal(result@.Data, c(0,0))

  # make sure that matrix input does same thing
  expect_equivalent(caliper(A, 2), caliper(m, 2))

  # values option returns the value from the object directly
  expect_equivalent(caliper(A,2) + A, caliper(A, 2, values = T))
})

test_that("Caliper exclusion", {
  m <- matrix(c(3,Inf, 1, 3), nrow = 2, ncol = 2,
              dimnames = list(treatment = c("A", "B"),
                              control = c("C", "D")))
  A <- as.InfinitySparseMatrix(m)

  # use the Mahalanobis distance match_on method
  result <- caliper(A, 2, exclude = c("B"))


  m2 <- matrix(c(Inf,Inf, 0, 0), nrow = 2, ncol = 2,
              dimnames = list(treatment = c("A", "B"),
                              control = c("C", "D")))

  expect_equal(as.matrix(result), m2)

})

test_that("Caliper respects groups", {

  # set up the exact match
  Z <- rep(c(T,F), each = 4)
  names(Z) <- c(LETTERS[1:4], letters[23:26])
  B <- c(T,T,F,F,T,T,F,F)
  em <- exactMatch(Z ~ B)

  expect_equal(length(findSubproblems(em)), 2)

  expect_equal(length(findSubproblems(caliper(em, 2))), 2)
  # f <- function(d) { fullmatch(d, min.controls = 1, omit.fraction = 0.75)}

  # expect_equal(sum(is.na(f(m))), 8) # should fail entirely
  # expect_equal(sum(is.na(f(em))), 0) # everything should work within strata

  # here is the real test, can we combine the two to firewall the failure?
  # expect_equal(sum(is.na(f(m + em))), 4)
})


test_that("update() caliper objects", {
  Z <- rep(c(0,1), each = 10)
  S <- rep(1:10 * 2, 2)
  names(Z) <- names(S) <- letters[1:20]

  basic <- caliper(match_on(S, z = Z), 2)
  expect_equal(length(basic), 28)

  S <- rep(1:10 * 3, 2)
  names(S) <- letters[1:20]
  updated <- update(basic)
  expect_equal(length(updated), 10)
})
