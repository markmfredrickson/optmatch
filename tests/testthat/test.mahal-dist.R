
context("Mahalanobis distances")

test_that("", {
  data(nuclear, package="boot")
  m1 <- match_on(pr ~ cap, data = nuclear)
  m2 <- match_on(pr ~ date + cum.n, data = nuclear)
  m3 <- match_on(pr ~ date + cum.n, data = nuclear, within = exactMatch(pr ~ pt, data = nuclear))
  expect_true("DenseMatrix" %in% class(m1))
  expect_true("DenseMatrix" %in% class(m2))
  expect_true("BlockedInfinitySparseMatrix" %in% class(m3))

  if ( (require(splines)) ) {
    m4 <- match_on(pr ~ ns(date,df=3) + cum.n, data = nuclear)
    m5 <- match_on(pr ~ ns(date,df=3) + cum.n, data = nuclear, exactMatch(pr ~ pt, data = nuclear))
    expect_true("DenseMatrix" %in% class(m4))
    expect_true("BlockedInfinitySparseMatrix" %in% class(m5))
  }

  cum.n.q <- cut(nuclear$cum.n, quantile(nuclear$cum.n), include.lowest=TRUE)
  m6 <- match_on(pr ~ date + cum.n.q, data = nuclear)
  m7 <- match_on(pr ~ date + cum.n.q, data = nuclear, within = exactMatch(pr~pt, data = nuclear))
  expect_true("DenseMatrix" %in% class(m6))
  expect_true("BlockedInfinitySparseMatrix" %in% class(m7))

  ### should give error, incorrect mode
  #expect_error(match_on(as.factor(pr) ~ cap, data = nuclear))
})
