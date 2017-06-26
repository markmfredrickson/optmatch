context("CBPS integration")

test_that("internal predict.CBPS function", {
    if (require('CBPS'))
        {
  cpsm <- CBPS(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n, ATT = 0, data = nuclearplants)
  expect_true(is(fullmatch(cpsm, data = nuclearplants), "optmatch"))

  # Check recognition of type argument.

  mo0 <- match_on(cpsm, data=nuclearplants)
  mo1 <- match_on(cpsm, data=nuclearplants, type="link")
  expect_equivalent(mo0,mo1)
  expect_equivalent(fullmatch(mo0, data=nuclearplants), fullmatch(cpsm, data=nuclearplants))

  mo2 <- match_on(cpsm, data=nuclearplants, type="response")
  expect_false(isTRUE(all.equal(mo1, mo2, check.attributes =FALSE)))
  expect_equivalent(fullmatch(mo2, data=nuclearplants), fullmatch(cpsm, data=nuclearplants, type="response"))
  
  # See if scores is properly imputing.

  nuclearplants$t1[c(2,5,10)] <- NA
  nuclearplants$ne[c(6,25,23)] <- NA

  cpsm <- CBPS::CBPS(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n, ATT = 0, data = nuclearplants)
  expect_true(length(fullmatch(cpsm, data = nuclearplants)) == nrow(nuclearplants))
}
})
