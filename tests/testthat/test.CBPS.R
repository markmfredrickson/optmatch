context("CBPS integration")

test_that("internal predict.CBPS function", {
  library(CBPS)
  cpsm <- CBPS(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n, ATT = 0, data = nuclearplants)
  expect_true(is(fullmatch(cpsm, data = nuclearplants), "optmatch"))
})
