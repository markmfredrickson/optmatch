library(testthat)

context('summary optmatch')

test_that("summary.optmatch", {
  data(plantdist)
  expect_warning(summary(fullmatch(1 * (plantdist < 10)))) # a zero-1 matrix
  expect_error(summary(pairmatch(plantdist + caliper(plantdist, 1)))) # Matching fails everywhere

  data(nuclearplants)
  psm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
  psd <- match_on(psm, standardization.scale = sd) # backwards compatible to 0.7-2
  psfm <- fullmatch(psd + caliper(psd, 0.25), data = nuclearplants)
  summary(psfm)
  ## not run as it causes an error in a subclass
  expect_error(pspm <- pairmatch(caliper(match_on(psm, standarization.scale = sd, within =
                     exactMatch(pr ~ pt, data = nuclearplants)), width=2))) # Fails in subclass '1'
  expect_error(summary(pspm), "object 'pspm' not found")
  psd[1,] <- psd[1,] + rep(100,22)

  # due to slight differences in the match on different platforms, just check that the
  # total.distances remain the same
  summary(pairmatch(psd, controls=2, data = nuclearplants))$total.distance

  summary(psfm, propensity.model=psm)
  require('RItools')
  summary(psfm, propensity.model='foo')
  summary(psfm, propensity.model=psm)
  summary(psfm, psm)
  psm2 <- glm(pr~ cut(date, c(67, 69.5, 72)) +
              t1 + t2 + cap + ne + ct + bw + cum.n + pt,
              family=binomial, data=nuclearplants)
  psd2 <- match_on(psm2, standardization.scale = sd)
  psd2summary <- summary(pairmatch(psd2, data = nuclearplants), propensity.model=psm2)

  # due to slight differences in the match on different platforms, just check that the
  # total.distances are the same and that the chi-squared value is 9.5 +- 0.5

  psd2summary$total.distance
  chisquared.value <- psd2summary$balance$overall$chisquare
  expect_true(abs(9.5 - chisquared.value) < 0.5)
})
