
context("maxControlsCap function old")

test_that("maxControlsCap", {
  data(nuclearplants)
  mhd2a <- match_on(pr ~ date + cum.n, data = nuclearplants,
                    within = exactMatch(pr ~ pt, data = nuclearplants))
  mhd2a <- t(mhd2a)

  mhd2a.caliper <- mhd2a + caliper(mhd2a, 3)
  s1 <- stratumStructure(fullmatch(mhd2a.caliper, data=nuclearplants)) # Works OK:
  expect_equal(names(s1), c("5:1", "4:1", "2:1", "1:1", "1:2"))
  expect_equal(as.vector(s1), c(1,2,3,2,1))
  expect_equal(attr(s1, "comparable.num.matched.pairs"), 12.2)

  mx1 <- maxControlsCap(mhd2a.caliper)              # no unmatchable Tx
  expect_true(all.equal(unlist(mx1),c(0, 0, .5, 1), check.attributes=FALSE))
  s2 <- stratumStructure(fullmatch(mhd2a.caliper, max=1, data=nuclearplants))
  s3 <- stratumStructure(fullmatch(mhd2a.caliper, max=1/2, data=nuclearplants))
  s4 <- stratumStructure(fullmatch(mhd2a + caliper(mhd2a, 2), data=nuclearplants)) # Problem in version <= .5-9:
  expect_equal(names(s2), c("5:1", "4:1", "2:1", "1:1"))
  expect_equal(as.vector(s2), c(1,2,2,5))
  expect_equal(attr(s2, "comparable.num.matched.pairs"), 12.53333333)
  expect_equal(names(s3), c("1:0", "0:1"))
  expect_equal(as.vector(s3), c(22,10))
  expect_equal(attr(s3, "comparable.num.matched.pairs"), 0)
  expect_equal(names(s4), c("5:1", "4:1", "2:1", "1:1", "1:2"))
  expect_equal(as.vector(s4), c(1,2,3,2,1))
  expect_equal(attr(s4, "comparable.num.matched.pairs"), 12.2)

  mx2 <- maxControlsCap(mhd2a + caliper(mhd2a, 2))     # caused by unmatchable Tx
  expect_true(all.equal(unlist(mx2),c(0, 0, .5, 1), check.attributes=FALSE))
  s5 <- stratumStructure(fullmatch(mhd2a + caliper(mhd2a, 2), max=mx2$strictest, data=nuclearplants))
  s6 <- stratumStructure(fullmatch(mhd2a + caliper(mhd2a, 2), max=1/2, data=nuclearplants))
  expect_equal(names(s5), c("1:0", "1:1", "0:1"))
  expect_equal(as.vector(s5), c(19,3,7))
  expect_equal(attr(s5, "comparable.num.matched.pairs"), 3)
  expect_equal(names(s6), c("1:0", "0:1"))
  expect_equal(as.vector(s6), c(22,10))
  expect_equal(attr(s6, "comparable.num.matched.pairs"), 0)
})
