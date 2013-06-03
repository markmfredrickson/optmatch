library(testthat)

context("more optmatch methods")

test_that("", {
  data(plantdist)
  expect_warning(plantsfm <- fullmatch(plantdist))
  plantsfm[1:10]
  attributes(plantsfm[1:10])
  attributes(plantsfm[5:10,drop=TRUE])
  plantsfm[1:26 <11]
  attributes(plantsfm[1:26 <6])
  attributes(plantsfm[1:26 <6,drop=TRUE])
  plantsfm[names(plantsfm)[1:10] ]
  attributes(plantsfm[names(plantsfm)[1:10] ])
  attributes(plantsfm[names(plantsfm)[5:10],drop=TRUE])


  plantsfm[5] <- "1.4"
  plantsfm[1:5]
  attributes(plantsfm)

  expect_warning(plantsfm <- fullmatch(plantdist))
  plantsfm[26:1]
  attributes(plantsfm[26:1])

  ### arises in lme4:::lmerFactorList , which is called in lme4::lmer
  ### at following line:
  ###
  ###   fl <- lapply(bars, function(x) eval(substitute(as.factor(fac)[,
  ###        drop = TRUE], list(fac = x[[3]])), mf))
  ###
  ### (caused [.optmatch to die in optmatch version 0.4-0 on R-2.6.0 +)
  p1 <- as.factor(plantsfm)[, drop = TRUE]
  p2 <- plantsfm[, drop = TRUE]
  p3 <- plantsfm[, drop = FALSE]
  expect_true(all(p1==p2))
  expect_true(all(p1==p3))
})
