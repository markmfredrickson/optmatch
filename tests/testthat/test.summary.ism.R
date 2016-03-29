context("summary method for ISM and related")

test_that("summary for ISM", {
  set.seed(1)
  d <- data.frame(z=rep(0:1, each=5),
                  x=rnorm(10))
  rownames(d) <- letters[1:10]
  m1 <- match_on(z ~ x, data=d)
  sm1 <- summary(m1)

  m2 <- m1 + caliper(m1, width=1)
  sm2 <- summary(m2)
  expect_true(is(sm2, "summary.InfinitySparseMatrix"))
  expect_true(is.list(sm2))
  expect_equal(attr(sm2, "ismname"), "m2")
  expect_equal(sm2$total$treatment, 5)
  expect_equal(sm2$total$control, 5)
  expect_equal(sm2$total$matchable, 12)
  expect_equal(sm2$total$unmatchable, 25-12)
  expect_equal(length(sm2$matchable$treatment), 5)
  expect_equal(length(sm2$matchable$control), 4)
  expect_equal(sm2$unmatchable$treatment, character(0))
  expect_equal(sm2$unmatchable$control, "d")
  expect_true(is(sm2$distances, "summaryDefault"))

  m3 <- m2
  m3[1:2] <- Inf
  sm3 <- summary(m3)

  expect_equal(sm3$total$matchable, 10)
  expect_equal(sm3$total$unmatchable, 25-10)
  expect_equal(length(sm3$matchable$treatment), 4)
  expect_equal(length(sm3$matchable$control), 4)
  expect_equal(sm3$unmatchable$treatment, "f")
  expect_equal(sm3$unmatchable$control, "d")
  expect_true(is(sm3$distances, "summaryDefault"))
  expect_true(all(is.finite(sm3$distances)))

  m4 <- m1 + caliper(m1, width=.0001)
  sm4 <- summary(m4)
  expect_equal(sm4$matchable$treatment, character(0))
  expect_equal(sm4$matchable$control, character(0))
  expect_true(is.null(sm4$distances))

  m5 <- m3
  m5@.Data <- rep(Inf, length(m5))
  sm5 <- summary(m5)
  expect_equal(sm5$matchable$treatment, character(0))
  expect_equal(sm5$matchable$control, character(0))
  expect_true(is.null(sm5$distances))


})

test_that("summary for BlockedISM", {
  set.seed(1)
  d <- data.frame(z=rep(0:1, each=5),
                  x=rnorm(10),
                  q=rep(c("a", "d"), times=5))
  rownames(d) <- letters[1:10]
  m1 <- match_on(z ~ x + strata(q), data=d, caliper=1)
  sm1 <- summary(m1)

  expect_true(is(sm1, "summary.BlockedInfinitySparseMatrix"))
  expect_true(is.list(sm1))
  expect_equal(length(sm1), 3)
  expect_equal(names(sm1), c("a", "d", "overall"))
  expect_equal(attr(sm1, "ismname"), "m1")
  expect_equal(attr(sm1, "blocknames"), c("a", "d"))
  expect_equal(attr(sm1, "printAllBlocks"), FALSE)
  expect_equal(attr(sm1, "blockStructure"), TRUE)
  expect_true(is(sm1[["a"]], "summary.InfinitySparseMatrix"))
  expect_true(is(sm1[["d"]], "summary.InfinitySparseMatrix"))

  sm2 <- summary(m1, printAllBlocks=TRUE, blockStructure=FALSE)
  expect_equal(attr(sm2, "printAllBlocks"), TRUE)
  expect_equal(attr(sm2, "blockStructure"), FALSE)

  expect_true(all.equal(c(5,5,4,21),
                        unlist(sm1$overall$total),
                        check.attributes=FALSE))

  # Alternate ways of calling blocks
  suma1 <- sm1[['a']]
  suma2 <- sm1$`a`
  suma3 <- sm1[[1]]
  expect_identical(suma1, suma2)
  expect_identical(suma1, suma3)


  expect_equal(attr(sm1$`a`, "blockname"), "a")
  expect_equal(attr(sm1$`d`, "blockname"), "d")

})


test_that("summary for DenseMatrix", {
  set.seed(1)
  d <- data.frame(z=rep(0:1, each=5),
                  x=rnorm(10))
  rownames(d) <- letters[1:10]
  m1 <- match_on(z ~ x, data=d)
  sm1 <- summary(m1)
  expect_true(is(sm1, "summary.DenseMatrix"))
  expect_true(is.list(sm1))
  expect_equal(attr(sm1, "ismname"), "m1")
  expect_equal(sm1$total$treatment, 5)
  expect_equal(sm1$total$control, 5)
  expect_equal(sm1$total$matchable, 25)
  expect_equal(sm1$total$unmatchable, 0)
  expect_equal(length(sm1$matchable$treatment), 5)
  expect_equal(length(sm1$matchable$control), 5)
  expect_equal(sm1$unmatchable$treatment, character(0))
  expect_equal(sm1$unmatchable$control, character(0))
  expect_true(is(sm1$distances, "summaryDefault"))

  m2 <- m1
  m2[1,] <- Inf
  sm2 <- summary(m2)
  expect_true(is(sm2, "summary.DenseMatrix"))
  expect_true(is.list(sm2))
  expect_equal(sm2$total$treatment, 5)
  expect_equal(sm2$total$control, 5)
  expect_equal(sm2$total$matchable, 20)
  expect_equal(sm2$total$unmatchable, 25-20)
  expect_equal(length(sm2$matchable$treatment), 4)
  expect_equal(length(sm2$matchable$control), 5)
  expect_equal(sm2$unmatchable$treatment, "f")
  expect_equal(sm2$unmatchable$control, character(0))
  expect_true(is(sm2$distances, "summaryDefault"))
})

test_that("distanceSummary suppresses distance", {
  set.seed(1)
  d <- data.frame(z=rep(0:1, each=5),
                  x=rnorm(10),
                  q=rep(c("a", "d"), times=5))
  rownames(d) <- letters[1:10]
  m1 <- match_on(z ~ x, data=d)

  expect_true(!is.null(summary(m1)$distances))
  expect_true(!is.null(summary(m1, distanceSummary=TRUE)$distances))
  expect_true(is.null(summary(m1, distanceSummary=FALSE)$distances))

  m2 <- match_on(z ~ x, data=d, caliper=1)
  expect_true(!is.null(summary(m2)$distances))
  expect_true(!is.null(summary(m2, distanceSummary=TRUE)$distances))
  expect_true(is.null(summary(m2, distanceSummary=FALSE)$distances))

  m3 <- match_on(z ~ x + strata(q), data=d, caliper=1)
  sm3.1 <- summary(m3)
  expect_true(!is.null(sm3.1$overall$distances))
  expect_true(!is.null(sm3.1$a$distances))
  expect_true(!is.null(sm3.1$d$distances))

  sm3.2 <- summary(m3, distanceSummary=TRUE)
  expect_true(!is.null(sm3.2$overall$distances))
  expect_true(!is.null(sm3.2$a$distances))
  expect_true(!is.null(sm3.2$d$distances))

  sm3.3 <- summary(m3, distanceSummary=FALSE)
  expect_true(is.null(sm3.3$overall$distances))
  expect_true(is.null(sm3.3$a$distances))
  expect_true(is.null(sm3.3$d$distances))
})
