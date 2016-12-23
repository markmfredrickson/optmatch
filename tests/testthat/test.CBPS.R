context("CBPS integration")

test_that("internal predict.CBPS function", {
    if (require('CBPS'))
        {
  cpsm <- CBPS(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n, ATT = 0, data = nuclearplants)
  expect_true(is(fullmatch(cpsm, data = nuclearplants), "optmatch"))

  # See if scores is properly imputing.

  nuclearplants$t1[c(2,5,10)] <- NA
  nuclearplants$ne[c(6,25,23)] <- NA

  cpsm <- CBPS::CBPS(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n, ATT = 0, data = nuclearplants)
  expect_true(length(fullmatch(cpsm, data = nuclearplants)) == nrow(nuclearplants))
}
})
