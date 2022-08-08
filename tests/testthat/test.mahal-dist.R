
context("Mahalanobis distances")

test_that("", {
  data(nuclearplants)
  m1 <- match_on(pr ~ cap, data = nuclearplants)
  m2 <- match_on(pr ~ date + cum.n, data = nuclearplants)
  m3 <- match_on(pr ~ date + cum.n, data = nuclearplants, within = exactMatch(pr ~ pt, data = nuclearplants))
  expect_true("DenseMatrix" %in% class(m1))
  expect_true("DenseMatrix" %in% class(m2))
  expect_true("BlockedInfinitySparseMatrix" %in% class(m3))

  if ( (require(splines)) ) {
    m4 <- match_on(pr ~ ns(date,df=3) + cum.n, data = nuclearplants)
    m5 <- match_on(pr ~ ns(date,df=3) + cum.n, data = nuclearplants, exactMatch(pr ~ pt, data = nuclearplants))
    expect_true("DenseMatrix" %in% class(m4))
    expect_true("BlockedInfinitySparseMatrix" %in% class(m5))
  }

  cum.n.q <- cut(nuclearplants$cum.n, quantile(nuclearplants$cum.n), include.lowest=TRUE)
  m6 <- match_on(pr ~ date + cum.n.q, data = nuclearplants)
  m7 <- match_on(pr ~ date + cum.n.q, data = nuclearplants, within = exactMatch(pr~pt, data = nuclearplants))
  expect_true("DenseMatrix" %in% class(m6))
  expect_true("BlockedInfinitySparseMatrix" %in% class(m7))

  ### should give error, incorrect mode
  expect_error(match_on(as.factor(pr) ~ cap, data = nuclearplants))
})
