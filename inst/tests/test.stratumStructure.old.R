library(testthat)

context('stratumStructure tests old')

test_that("stratumStructure old", {

  data(plantdist)

  expect_warning(plantsfm <- fullmatch(plantdist))
  expect_warning(plantsfm1 <- fullmatch(plantdist,min.controls=2, max.controls=3))

  stratumStructure(plantsfm)
  stratumStructure(plantsfm1)

  data(nuclearplants)
  psd <- match_on(glm(pr~.-(pr+cost), family=binomial(),
                      data=nuclearplants),
                  standardization.scale = sd) # sd was standard in < 0.7-2
  stratumStructure(fullmatch(psd, data=nuclearplants))
  stratumStructure(fm.psd.cal <- fullmatch(psd/(psd<.25), data=nuclearplants))

  psd2 <- match_on(glm(pr~.-(pr+cost), family=binomial(),
                       data=nuclearplants),
                   standardization.scale = sd,
                   within = exactMatch(pr ~ pt, nuclearplants))
  s1 <- stratumStructure(fullmatch(psd2,min.controls=1,max.controls=1, data=nuclearplants))
  s2 <- stratumStructure(fullmatch(psd2,min.controls=3,max.controls=3, data=nuclearplants))
  expect_true(all(sapply(strsplit(names(s1), ":"), "[", 2) <= 1))
  expect_true(all(sapply(strsplit(names(s2), ":"), "[", 2) <= 3))


  ### Tests of min.controls, max.controls
  s1 <- stratumStructure(fm.psd.cal, min.controls=.5)
  s2 <- stratumStructure(fm.psd.cal, max.controls=3)
  s3 <- stratumStructure(fm.psd.cal, min.controls=.5, max.controls=3)

  expect_true(all(substr(lapply(strsplit(names(s1), ":"), "[", 1),1,1) <= 2))
  expect_true(all(substr(lapply(strsplit(names(s2), ":"), "[", 2),1,1) <= 3))
  expect_true(all(substr(lapply(strsplit(names(s3), ":"), "[", 1),1,1) <= 2))
  expect_true(all(substr(lapply(strsplit(names(s3), ":"), "[", 2),1,1) <= 3))


  ### tests of stratumStructure on non-optmatch objects
  #fac <- as.factor(fm.psd.cal)
  #tvec <- attr(fm.psd.cal, "contrast.group")
  #stratumStructure(fac, tvec)
  #stratumStructure(as.integer(fac),tvec)
})
