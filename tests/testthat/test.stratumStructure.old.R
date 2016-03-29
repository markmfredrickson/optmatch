
context('stratumStructure tests old')

test_that("stratumStructure old", {

  data(plantdist)

  expect_warning(plantsfm <- fullmatch(plantdist))
  expect_warning(plantsfm1 <- fullmatch(plantdist,min.controls=2, max.controls=3))

  s1 <- stratumStructure(plantsfm)
  s2 <- stratumStructure(plantsfm1)

  expect_equal(names(s1), c("1:1", "1:2", "1:3", "1:4", "1:6"))
  expect_equal(as.vector(s1), c(2,2,1,1,1))
  expect_true(all.equal(attr(s1, "comparable.num.matched.pairs"), 9.480952, tol=1e-6))

  expect_equal(names(s2), c("1:2", "1:3"))
  expect_equal(as.vector(s2), c(2,5))
  expect_true(all.equal(attr(s2, "comparable.num.matched.pairs"), 10.1666666))

  data(nuclearplants)
  psd <- match_on(glm(pr~.-(pr+cost), family=binomial(),
                      data=nuclearplants),
                  standardization.scale = sd) # sd was standard in < 0.7-2
  s3 <- stratumStructure(fullmatch(psd, data=nuclearplants))
  expect_equal(names(s3), c("5:1", "3:1", "1:3", "1:17"))
  expect_equal(as.vector(s3), rep(1,4))
  expect_true(all.equal(attr(s3, "comparable.num.matched.pairs"), 6.5555555))

  s4 <- stratumStructure(fm.psd.cal <- fullmatch(psd/(psd<.25), data=nuclearplants))
  expect_equal(names(s4), c("1:0", "3:1", "2:1", "1:2", "1:7", "0:1"))
  expect_equal(as.vector(s4), c(3,1,1,1,1,11))
  expect_true(all.equal(attr(s4, "comparable.num.matched.pairs"), 5.9166666))

  psd2 <- match_on(glm(pr~.-(pr+cost), family=binomial(),
                       data=nuclearplants),
                   standardization.scale = sd,
                   within = exactMatch(pr ~ pt, nuclearplants))
  s5 <- stratumStructure(fullmatch(psd2,min.controls=1,max.controls=1, data=nuclearplants))
  s6 <- stratumStructure(fullmatch(psd2,min.controls=3,max.controls=3, data=nuclearplants))
  expect_true(all(sapply(strsplit(names(s5), ":"), "[", 2) <= 1))
  expect_true(all(sapply(strsplit(names(s6), ":"), "[", 2) <= 3))


  ### Tests of min.controls, max.controls
  s7 <- stratumStructure(fm.psd.cal, min.controls=.5)
  s8 <- stratumStructure(fm.psd.cal, max.controls=3)
  s9 <- stratumStructure(fm.psd.cal, min.controls=.5, max.controls=3)

  expect_true(all(substr(lapply(strsplit(names(s7), ":"), "[", 1),1,1) <= 2))
  expect_true(all(substr(lapply(strsplit(names(s8), ":"), "[", 2),1,1) <= 3))
  expect_true(all(substr(lapply(strsplit(names(s9), ":"), "[", 1),1,1) <= 2))
  expect_true(all(substr(lapply(strsplit(names(s9), ":"), "[", 2),1,1) <= 3))


  ### tests of stratumStructure on non-optmatch objects
  #fac <- as.factor(fm.psd.cal)
  #tvec <- attr(fm.psd.cal, "contrast.group")
  #stratumStructure(fac, tvec)
  #stratumStructure(as.integer(fac),tvec)
})
