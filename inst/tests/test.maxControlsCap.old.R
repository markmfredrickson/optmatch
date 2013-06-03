library(testthat)

context("maxControlsCap function old")

test_that("maxControlsCap", {
  data(nuclearplants)
  mhd2a <- match_on(pr ~ date + cum.n, data = nuclearplants,
                    within = exactMatch(pr ~ pt, data = nuclearplants))
  mhd2a <- t(mhd2a)

  mhd2a.caliper <- mhd2a + caliper(mhd2a, 3)
  s1 <- stratumStructure(fullmatch(mhd2a.caliper, data=nuclearplants)) # Works OK:
  expect_true(attr(s1, 'comparable.num.matched.pairs') > 0)

  mx1 <- maxControlsCap(mhd2a.caliper)              # no unmatchable Tx
  s2 <- stratumStructure(fullmatch(mhd2a.caliper, max=1, data=nuclearplants))
  s3 <- stratumStructure(fullmatch(mhd2a.caliper, max=1/2, data=nuclearplants))
  s4 <- stratumStructure(fullmatch(mhd2a + caliper(mhd2a, 2), data=nuclearplants)) # Problem in version <= .5-9:
  expect_true(attr(s2, 'comparable.num.matched.pairs') > 0)
  expect_true(attr(s3, 'comparable.num.matched.pairs') == 0)
  expect_true(attr(s4, 'comparable.num.matched.pairs') > 0)

  mx2 <- maxControlsCap(mhd2a + caliper(mhd2a, 2))     # caused by unmatchable Tx
  s5 <- stratumStructure(fullmatch(mhd2a + caliper(mhd2a, 2), max=mx2$strictest, data=nuclearplants))
  s6 <- stratumStructure(fullmatch(mhd2a + caliper(mhd2a, 2), max=1/2, data=nuclearplants))
  expect_true(attr(s5, 'comparable.num.matched.pairs') > 0)
  expect_true(attr(s6, 'comparable.num.matched.pairs') == 0)
})
